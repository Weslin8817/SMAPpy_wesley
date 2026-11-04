"""Precision-weighted centerlines, independent of GUI and localization storage.

An open, longitudinal curve describes the center of a filled distribution;
there is no shell model, fitted radius, or precision correction to distances.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline, PPoly
from scipy.integrate import cumulative_trapezoid


RADIAL_PERCENTILES = (30, 50, 70, 90, 95, 100)


def radial_percentiles(distances):
    values = np.percentile(distances, RADIAL_PERCENTILES) if len(distances) else [None]*len(RADIAL_PERCENTILES)
    return {f"R{p}_nm": None if value is None else float(value)
            for p, value in zip(RADIAL_PERCENTILES, values)}


def result_summary(entry):
    """Shared table/export values, including results accumulated with old R keys."""
    values = {key: value for key, value in entry.summary.items()
              if key not in {f"R{p}_nm" for p in (*RADIAL_PERCENTILES, 80)}}
    values.update(radial_percentiles(entry.distances))
    return values


@dataclass
class CenterlineResult:
    points: np.ndarray
    indices: np.ndarray
    curve: np.ndarray
    distances: np.ndarray
    arc: np.ndarray
    nearest: np.ndarray
    summary: dict


def manual_distance_mask(distances, maximum_distance_nm):
    """Select an inclusive radius about the initial centerline, in nm."""
    if not np.isfinite(maximum_distance_nm) or maximum_distance_nm <= 0:
        raise ValueError("Maximum distance must be a positive finite radius in nm")
    return np.asarray(distances) <= maximum_distance_nm


def conservative_distance_mask(distances):
    """Independent automatic cutoff: median + four scaled MADs.

    A 1 nm scale floor avoids a vanishing threshold for an exactly straight
    or quantized synthetic cloud. No assumed structure radius is used.
    """
    distances = np.asarray(distances, dtype=float)
    median = np.median(distances)
    scale = max(1., 1.4826 * np.median(np.abs(distances-median)))
    cutoff = float(median + 4.*scale)
    return distances <= cutoff, cutoff


def fit_selected_centerline(points, precision=None, *, bins=16, smoothing=1.,
                            rejection="None", maximum_distance_nm=50.):
    """One initial fit, one optional selection, at most one final fit.

    Modes are mutually exclusive. The mask is frozen before the final fit;
    distances are then recomputed, never clipped or iteratively screened.
    Returns (result, effective cutoff or None, excluded finite point count).
    """
    if rejection not in ("None", "Manual maximum distance", "Auto conservative"):
        raise ValueError("Unknown noise rejection mode")
    if rejection == "Manual maximum distance":
        manual_distance_mask([], maximum_distance_nm)  # validate before fitting
    initial = fit_centerline(points, precision, bins=bins, smoothing=smoothing)
    if rejection == "None":
        return initial, None, 0
    if rejection == "Manual maximum distance":
        cutoff = float(maximum_distance_nm)
        keep = manual_distance_mask(initial.distances, cutoff)
    else:
        keep, cutoff = conservative_distance_mask(initial.distances)
    excluded = int((~keep).sum())
    if not excluded:
        return initial, cutoff, 0
    if keep.sum() < 48:
        raise ValueError(f"{rejection} retained only {int(keep.sum())} localizations; "
                         "at least 48 are required. Increase maximum distance or use None.")
    selected_indices = initial.indices[keep]
    # Resolve missing uncertainty once from the input, so point selection
    # cannot silently change its fallback precision between the two fits.
    errors, _ = sanitize_precision(precision, len(points), initial.points.shape[1])
    result = fit_centerline(initial.points[keep], errors[selected_indices],
                            bins=bins, smoothing=smoothing)
    result.indices = selected_indices
    result.summary["excluded_nonfinite"] = initial.summary["excluded_nonfinite"]
    result.summary["precision_replacements"] = initial.summary["precision_replacements"]
    # A rejection boundary was defined about the initial curve. Use a
    # conservative inward margin when judging a radius about the final curve.
    result.summary["rejection_refit_shift_nm"] = float(np.max(
        np.linalg.norm(result.nearest-initial.nearest[keep], axis=1)))
    return result, cutoff, excluded


def sanitize_precision(precision, n, dim):
    """Invalid errors use axis median, then 20 nm XY / 60 nm Z.

    Errors outside [0.1, 10000] nm are invalid. A 1 nm floor prevents a
    nearly zero uncertainty dominating an entire bin. These are fitting
    safeguards, not measurement corrections.
    """
    p = np.full((n, dim), np.nan) if precision is None else np.array(precision, float, copy=True)
    if p.shape != (n, dim):
        raise ValueError("precision must have one uncertainty per point and axis")
    replaced = 0
    for j in range(dim):
        good = np.isfinite(p[:, j]) & (p[:, j] >= .1) & (p[:, j] <= 10000)
        fallback = np.median(p[good, j]) if good.any() else (60. if j == 2 else 20.)
        replaced += int((~good).sum())
        p[~good, j] = fallback
    return np.maximum(p, 1.), replaced


def _spatial_median(points):
    center = points.mean(axis=0)
    for _ in range(60):
        w = 1 / np.maximum(np.linalg.norm(points - center, axis=1), 1e-9)
        new = np.average(points, axis=0, weights=w)
        if np.linalg.norm(new - center) < 1e-8:
            break
        center = new
    return center


def project_curve(points, spline):
    """Global Euclidean projection onto every continuous cubic span.

    On each span, stationary squared distances are roots of a quintic.
    Test all real roots within that span plus endpoints. Unlike nearest
    sample distances or a single local optimizer, this includes every
    candidate on the fitted curve. Arc integrals use a dense speed grid.
    """
    points = np.asarray(points, dtype=float)
    polys = [PPoly.from_spline((spline.t, spline.c[:, j], 3))
             for j in range(points.shape[1])]
    best = np.full(len(points), np.inf)
    parameter = np.zeros(len(points))
    nearest = np.zeros_like(points)
    for k, (left, right) in enumerate(zip(polys[0].x[:-1], polys[0].x[1:])):
        if right <= left or left < 0 or right > 1:
            continue
        coeff = np.stack([p.c[:, k] for p in polys], axis=1)
        derivative = coeff[:-1] * np.array([3, 2, 1])[:, None]
        for i, point in enumerate(points):
            residual = coeff.copy()
            residual[-1] -= point
            stationary = sum((np.convolve(residual[:, j], derivative[:, j])
                              for j in range(points.shape[1])))
            roots = np.roots(stationary)
            real = roots.real[np.abs(roots.imag) < 1e-9]
            u = np.r_[0., real[(real > 0) & (real < right-left)], right-left]
            candidates = np.stack([np.polyval(coeff[:, j], u)
                                   for j in range(points.shape[1])], axis=1)
            dist2 = np.sum((candidates-point)**2, axis=1)
            pick = np.argmin(dist2)
            if dist2[pick] < best[i]:
                best[i] = dist2[pick]
                parameter[i] = left + u[pick]
                nearest[i] = candidates[pick]
    grid = np.linspace(0, 1, 4097)
    arc = cumulative_trapezoid(np.linalg.norm(spline(grid, 1), axis=1), grid, initial=0)
    return np.sqrt(best), np.interp(parameter, grid, arc), nearest, arc[-1]


def fit_centerline(points, precision=None, *, bins=16, smoothing=1.0):
    """Robust PCA, local anisotropic centers, penalized low-rank B-spline.

    At most 12 coefficients per coordinate and at least 12 points per bin.
    The PCA ordering assumes an unbranched curve without a longitudinal
    fold. Endpoints are the 0.5/99.5 percentiles of its robust axial support.
    Outliers remain in the returned distribution, even if downweighted in
    estimation. Nonfinite coordinates alone are excluded.
    """
    points = np.asarray(points, float)
    if points.ndim != 2 or points.shape[1] not in (2, 3):
        raise ValueError("coordinates must be an N x 2 or N x 3 array in nm")
    if not 4 <= bins <= 64 or int(bins) != bins or not np.isfinite(smoothing) or not .01 <= smoothing <= 100:
        raise ValueError("bins must be 4–64; smoothing must be 0.01–100")
    valid = np.isfinite(points).all(axis=1)
    indices = np.flatnonzero(valid)
    p, replacements = sanitize_precision(precision, len(points), points.shape[1])
    p = p[valid]
    points = points[valid].copy()
    if len(points) < 48:
        raise ValueError("at least 48 finite localizations are required")
    origin = _spatial_median(points)
    centered = points-origin
    radius = np.linalg.norm(centered, axis=1)
    if np.median(radius) < 1e-9:
        raise ValueError("point cloud has insufficient longitudinal extension")
    robust = np.minimum(1., (2.5 * np.median(radius) / np.maximum(radius, 1e-9))**4)
    w = robust / np.mean(p*p, axis=1)
    eig, vectors = np.linalg.eigh((centered*w[:, None]).T @ centered / w.sum())
    if eig[-1] < 1. or eig[-1] < 2. * eig[-2]:
        raise ValueError("point cloud has insufficient longitudinal extension (PCA ratio < 2)")
    axis = vectors[:, -1]
    # A reproducible sign; reversal changes only the arc origin, not distances.
    if axis[np.argmax(np.abs(axis))] < 0:
        axis = -axis
    longitudinal = centered @ axis
    support = longitudinal[robust > .5]
    low, high = np.quantile(support, [.005, .995])
    if high-low < 1.:
        raise ValueError("longitudinal extent is less than 1 nm")
    t = (longitudinal-low)/(high-low)
    count = min(int(bins), len(points)//12)
    # Quantize only bin membership, not coordinates. Roundoff at a bin
    # boundary must not change a fit after translation or rotation.
    scaled = np.round(t*count, 10)
    membership = np.floor(scaled).astype(int)
    membership[scaled == count] = count-1
    centers, positions, information = [], [], []
    for bin_index in range(count):
        mask = (membership == bin_index) & (scaled >= 0) & (scaled <= count)
        if mask.sum() < 6:
            continue
        q, sig = centered[mask], p[mask]
        center = _spatial_median(q)
        scale = max(np.median(np.linalg.norm(q-center, axis=1)), 1.)
        # Shared robust residual factor preserves rotational equivariance
        # for isotropic errors; coordinate inverse variances retain anisotropy.
        for _ in range(30):
            r = np.linalg.norm(q-center, axis=1)
            weights = np.minimum(1., 1.5*scale/np.maximum(r, 1e-9))[:, None]/sig**2
            new = np.sum(weights*q, axis=0)/weights.sum(axis=0)
            if np.linalg.norm(new-center) < 1e-8:
                center = new
                break
            center = new
        centers.append(center)
        positions.append((center @ axis-low)/(high-low))
        information.append(weights.sum(axis=0))
    if len(centers) < 4:
        raise ValueError("too few populated longitudinal bins; reduce bins or redraw ROI")
    centers, positions, information = map(np.asarray, (centers, positions, information))
    ncoeff = min(12, max(4, len(centers)//2+1))
    knots = np.r_[np.zeros(4), np.linspace(0, 1, ncoeff-2)[1:-1], np.ones(4)]
    basis = BSpline(knots, np.eye(ncoeff), 3)(positions)
    # Integrated squared curvature penalizes shape, leaving straight lines
    # unpenalized. The common weight scale keeps anisotropic information.
    grid = np.linspace(0, 1, 257)
    d2 = BSpline(knots, np.eye(ncoeff), 3)(grid, 2)
    penalty = d2.T @ d2 / len(grid)
    information /= np.median(information)
    coeff = np.empty((ncoeff, points.shape[1]))
    try:
        for j in range(points.shape[1]):
            weights = information[:, j]
            coeff[:, j] = np.linalg.solve(
                basis.T @ (weights[:, None]*basis) + smoothing*.0001*penalty,
                basis.T @ (weights*centers[:, j]))
        spline = BSpline(knots, coeff, 3)
        distances, arc, nearest, length = project_curve(centered, spline)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as e:
        raise ValueError(f"centerline fitting/projection failed: {e}") from e
    curve = spline(np.linspace(0, 1, 1001)) + origin
    if not np.isfinite(distances).all() or not np.isfinite(curve).all() or not np.isfinite(length):
        raise ValueError("centerline fitting produced nonfinite results")
    q25, median, q75 = np.percentile(distances, [25, 50, 75])
    summary = dict(dimension=f"{points.shape[1]}D", N=len(points),
                   centerline_length_nm=float(length), mean_radial_distance_nm=float(distances.mean()),
                   median_radial_distance_nm=float(median), rms_radial_distance_nm=float(np.sqrt(np.mean(distances**2))),
                   sd_nm=float(distances.std()), iqr_nm=float(q75-q25), **radial_percentiles(distances),
                   bins=int(bins), populated_bins=len(centers),
                   coefficients=ncoeff, smoothing=float(smoothing), excluded_nonfinite=int((~valid).sum()),
                   precision_replacements=replacements, algorithm="robust_binned_pspline_v1")
    return CenterlineResult(points, indices, curve, distances, arc, nearest+origin, summary)


class ROIResults:
    """Plugin-owned snapshots. Identifiers never get reused after deletion."""
    def __init__(self):
        self.entries = {}
        self._next = 1

    def add(self, result, **provenance):
        identifier = f"ROI_{self._next:03d}"
        self._next += 1
        result.summary = dict(ROI=identifier, **result.summary, **provenance)
        self.entries[identifier] = result
        return identifier

    def delete(self, identifiers):
        for identifier in identifiers:
            self.entries.pop(identifier, None)

    def clear(self):
        self.entries.clear()

    def export(self, directory):
        entries = list(self.entries.values())
        if not entries:
            raise ValueError("no accumulated ROI results to export")
        directory = Path(directory)
        summary_path = directory / "centerline_radial_summary.csv"
        distance_path = directory / "centerline_radial_all_distances.csv"
        from .centerline_export import FIELDS, summary_row, write_dictionary
        summaries = [summary_row(result_summary(entry)) for entry in entries]
        with summary_path.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(summaries)
        with distance_path.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.writer(handle)
            writer.writerow(["ROI", "distance_nm", "localization_index", "dimension",
                             "centerline_arc_coordinate_nm", "nearest_centerline_x_nm",
                             "nearest_centerline_y_nm", "nearest_centerline_z_nm", "in_density_domain"])
            for entry in entries:
                for i, distance in enumerate(entry.distances):
                    xyz = list(entry.nearest[i])
                    if len(xyz) == 2:
                        xyz.append("")
                    writer.writerow([entry.summary["ROI"], distance, entry.indices[i],
                                     entry.summary["dimension"], entry.arc[i], *xyz,
                                     bool(entry.density_support[i]) if hasattr(entry, "density_support") else ""])
        segment_rows = [dict(ROI=entry.summary["ROI"],
                            density_dimension=entry.summary.get("density_dimension", ""),
                            density_units=entry.summary.get("density_units", ""), **row) for entry in entries
                        for row in getattr(entry, "density_segments", [])]
        segment_path = directory / "centerline_density_segments.csv"
        paths = [summary_path, distance_path]
        # Re-export after deleting density results must not leave old segment
        # rows beside a new summary. The GUI confirms overwriting all three.
        if segment_rows or segment_path.exists():
            with segment_path.open("w", newline="", encoding="utf-8-sig") as handle:
                fields = list(dict.fromkeys(k for row in segment_rows for k in row)) or ["ROI", "segment"]
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(segment_rows)
            paths.append(segment_path)
        for dimension in (2, 3):
            density_path = directory / f"centerline_density_{dimension}d_distances.csv"
            density_entries = [entry for entry in entries if hasattr(entry, "density_entry")
                               and entry.density_entry.points.shape[1] == dimension]
            if density_entries or density_path.exists():
                with density_path.open("w", newline="", encoding="utf-8-sig") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(["ROI", "grouped_table_index", f"distance_{dimension}d_nm", f"arc_{dimension}d_nm",
                                     "nearest_x_nm", "nearest_y_nm", "nearest_z_nm", "in_density_domain"])
                    for entry in density_entries:
                        fitted = entry.density_entry
                        for i in range(len(fitted.points)):
                            nearest = list(fitted.nearest[i]) + ([""] if dimension == 2 else [])
                            writer.writerow([entry.summary["ROI"], fitted.indices[i], fitted.distances[i],
                                             fitted.arc[i], *nearest, bool(fitted.density_support[i])])
                paths.append(density_path)
        dictionary_path = directory / "centerline_summary_fields.csv"
        write_dictionary(dictionary_path)
        paths.append(dictionary_path)
        return tuple(paths)

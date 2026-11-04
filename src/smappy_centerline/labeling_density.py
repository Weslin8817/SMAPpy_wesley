"""Point-supported geometry, with explicit quality gates, for sampling estimates.

Alpha is a maximum circumradius (a length, not its reciprocal). Summing the
retained Delaunay simplices measures their union, including disconnected parts
and excluding holes. It does not infer a biological boundary from empty space.
"""
from __future__ import annotations

import math
import numpy as np


def density_nyquist(count, measure_um, dimension):
    """Return count/um**d and the sampling length in nm."""
    if dimension not in (2, 3):
        raise ValueError("dimension must be 2 or 3")
    if not np.isfinite(measure_um) or measure_um <= 0:
        raise ValueError("area/volume must be finite and positive")
    if not np.isfinite(count) or count <= 0:
        raise ValueError("count must be finite and positive")
    density = count / measure_um
    return float(density), float(2000 / density ** (1 / dimension))


def reconstruct(points_nm, method="Alpha Shape", radius_nm=0.0,
                auto_factor=2.0, sensitivity_limit=0.25, max_points=100000,
                min_coverage=0.95):
    """Return diagnostics; ``measure_um`` is None when quality gates fail.

    Coordinates are deduplicated only for geometry, never for event counts.
    No jitter, imputation or silent convex-hull fallback. ``support_mask``
    identifies input events supported by the accepted shape, so callers must
    use that same subset in density numerators, including duplicate events.
    The +/-20% radius sweep reuses one tessellation. Quality thresholds are
    conservative heuristics, not confidence intervals or proof of validity.
    """
    from scipy.spatial import ConvexHull, Delaunay, QhullError, cKDTree

    p = np.asarray(points_nm, dtype=float)
    if p.ndim != 2 or p.shape[1] not in (2, 3):
        raise ValueError("points must have shape (n, 2) or (n, 3)")
    if method not in ("Alpha Shape", "Convex Hull"):
        raise ValueError("unknown geometry method")
    if (not np.isfinite([radius_nm, auto_factor, sensitivity_limit, min_coverage]).all()
            or radius_nm < 0 or auto_factor <= 0 or not 0 < sensitivity_limit <= 1
            or not 0 < min_coverage <= 1):
        raise ValueError("invalid alpha parameters")
    d = p.shape[1]
    out = dict(measure_um=None, dimension=d, method=method, status="unavailable",
               n_unique=0, radius_nm=None, coverage=None,
               sensitivity=None, sweep_low_um=None, sweep_high_um=None,
               event_coverage=None, min_coverage=min_coverage,
               nn_median_nm=None, nn_max_nm=None, nn_max_median_ratio=None,
               n_nn_flagged=0, auto_seed_nm=None, auto_max_nm=None,
               radius_scan=[], candidate_measure_um=None,
               components=None, alpha_hull_ratio=None,
               support_mask=np.zeros(len(p), dtype=bool),
               warnings=[])

    def fail(message):
        out["warnings"].append(message)
        return out

    if not np.isfinite(p).all():
        return fail("Non-finite coordinates: fix the filter/data; no rows were silently dropped.")
    p, inverse, multiplicity = np.unique(p, axis=0, return_inverse=True,
                                       return_counts=True)
    out["n_unique"] = len(p)
    if len(p) < d + 1:
        return fail("Too few distinct points for full-dimensional geometry.")
    if len(p) > max_points:
        return fail(f"More than {max_points} unique points; narrow the ROI/filter (no subsampling).")
    # Translation and isotropic scaling improve conditioning without changing
    # distances or using a bounding box as the estimated geometry.
    p = p - p.mean(axis=0)
    scale = float(np.max(np.abs(p)))
    if scale == 0 or np.linalg.matrix_rank(p / scale) < d:
        return fail("Degenerate geometry: collinear/coplanar points.")
    p /= scale
    singular = np.linalg.svd(p, compute_uv=False)
    if singular[-1] / singular[0] < 1e-6:
        return fail("Near-degenerate geometry: use a lower-dimensional model.")
    if len(p) < 10 * d:
        return fail(f"Sparse geometry: require at least {10*d} distinct points; this is a heuristic minimum.")
    try:
        nn = cKDTree(p).query(p, k=2)[0][:, 1]
        median_nn = float(np.median(nn))
        flagged = int(np.count_nonzero(nn > 10 * median_nn))
        out.update(nn_median_nm=median_nn * scale, nn_max_nm=float(np.max(nn)) * scale,
                   nn_max_median_ratio=float(np.max(nn) / median_nn), n_nn_flagged=flagged)
        if flagged:
            out["warnings"].append(
                f"Uneven spacing: {flagged} unique points have nearest-neighbour distance "
                ">10 times the median. Repeated blinking, sparse regions or outliers "
                "can cause this; no points were removed by this diagnostic.")
        if method == "Convex Hull":
            measure = float(ConvexHull(p).volume * (scale / 1000) ** d)
            if not np.isfinite(measure) or measure <= 0:
                return fail("Convex hull has no finite positive area/volume.")
            out.update(measure_um=measure, coverage=1.0, event_coverage=1.0,
                       support_mask=np.ones(len(inverse), dtype=bool), status="convex envelope comparison")
            out["warnings"].append("Convex hull fills holes and concavities and bridges components; not occupied geometry.")
            if flagged:
                out["warnings"].append("Convex envelope includes all selected points and may be inflated by isolated points; inspect the ROI before interpreting density.")
            return out
        tri = Delaunay(p)
        simplex = p[tri.simplices]
        edges = simplex[:, 1:] - simplex[:, :1]
        det = np.abs(np.linalg.det(edges))
        good = det > np.finfo(float).eps * 100
        radii = np.full(len(edges), np.inf)
        # Relative circumcentres avoid cancellation from large stage offsets.
        centres = np.linalg.solve(2 * edges[good],
                                  np.sum(edges[good] ** 2, axis=2)[..., None])[..., 0]
        radii[good] = np.linalg.norm(centres, axis=1)
        measures = det / math.factorial(d) * (scale / 1000) ** d
        # Each vertex first acquires d-dimensional support at the smallest
        # circumradius of its incident simplices. This is not a molecular
        # identity estimate. Unlike median NN alone it also sees the sparsely
        # sampled portion of a cloud containing close repeat localizations.
        support_radius = np.full(len(p), np.inf)
        for j in range(d + 1):
            np.minimum.at(support_radius, tri.simplices[:, j], radii)
        if radius_nm:
            candidates = [radius_nm / scale]
        else:
            needed = int(np.ceil(min_coverage * len(p))) - 1
            coverage_radius = float(np.partition(support_radius, needed)[needed])
            if not np.isfinite(coverage_radius):
                return fail("Numerical tessellation cannot support the requested coverage; inspect geometry or lower minimum coverage.")
            seed = max(auto_factor * median_nn, coverage_radius * 1.05)
            candidates = seed * np.array([1, 1.25, 1.5, 2, 3, 4])
            out.update(auto_seed_nm=seed * scale, auto_max_nm=float(candidates[-1]) * scale)
        # Do not increase a manual radius, and never search without a bound.
        # A stable result may still bridge holes; report hull ratio/components.
        for radius in candidates:
            keep = radii <= radius
            support = support_radius <= radius
            measure = float(measures[keep].sum())
            low = float(measures[radii <= radius * .8].sum())
            high = float(measures[radii <= radius * 1.2].sum())
            coverage = float(np.mean(support))
            event_coverage = float(multiplicity[support].sum() / len(inverse))
            sensitivity = (high - low) / measure if measure > 0 else None
            accepted = (np.isfinite(measure) and measure > 0
                        and coverage >= min_coverage and event_coverage >= min_coverage
                        and sensitivity is not None and sensitivity <= sensitivity_limit)
            out["radius_scan"].append(dict(radius_nm=float(radius * scale),
                measure_um=measure, coverage=coverage, event_coverage=event_coverage,
                sensitivity=sensitivity, accepted=bool(accepted)))
            out.update(radius_nm=float(radius * scale), coverage=coverage,
                       event_coverage=event_coverage, sensitivity=sensitivity,
                       sweep_low_um=low, sweep_high_um=high, candidate_measure_um=measure)
            if accepted:
                break
        if not accepted:
            if not np.isfinite(measure) or measure <= 0:
                out["warnings"].append("Alpha Shape contains no finite positive area/volume.")
            if coverage < min_coverage or event_coverage < min_coverage:
                out["warnings"].append("Alpha support coverage is below the requested minimum for unique points or localization events.")
            if sensitivity is not None and sensitivity > sensitivity_limit:
                out["warnings"].append("Alpha Shape is radius-sensitive: +/-20% sweep exceeds the allowed relative change.")
            return fail("No accepted radius. Inspect radius_scan; choose a physically justified manual radius or Convex Hull comparison. Diagnostic candidate measures are not accepted density denominators.")
        out.update(measure_um=measure, support_mask=support[inverse],
                   status="exploratory shape" if coverage == 1 else "exploratory partial-support shape")
        if coverage < 1:
            out["warnings"].append("Partial support: density uses only events represented in retained simplices. ROI totals are reported separately; this is not a whole-ROI density.")
        if len(out["radius_scan"]) > 1:
            out["warnings"].append("Automatic radius increased within the bounded scan; inspect the scan for filling of gaps/holes.")
        # Components joined through any vertex; no extra triangulation needed.
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        vertices = tri.simplices[keep]
        rows = np.repeat(vertices[:, 0], d)
        cols = vertices[:, 1:].ravel()
        graph = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(p), len(p))).tocsr()
        out["components"] = int(connected_components(graph[support][:, support], directed=False, return_labels=False))
        hull_measure = float(measures.sum())
        out["alpha_hull_ratio"] = measure / hull_measure if hull_measure > 0 else None
        if out["components"] > 1:
            out["warnings"].append(f"Shape has {out['components']} disconnected components; blinking clusters or separate structures can both cause this.")
        if np.all(keep[good]):
            out["warnings"].append("All resolvable simplices retained: this radius gives a convex envelope and may fill holes.")
        out["warnings"].append("Radius stability does not establish a true boundary. Repeated blinking, noise, holes and nonuniform sampling can bias geometry.")
        return out
    except (QhullError, np.linalg.LinAlgError) as exc:
        return fail(f"Geometry reconstruction failed ({type(exc).__name__}); no fallback used.")

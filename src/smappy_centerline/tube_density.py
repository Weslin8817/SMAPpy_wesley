"""Exploratory filled-strip (2D) or filled-tube (3D) centerline estimates.

Local second moments assume uniform intervals or circular cross-sections.
No localization-precision fallback is treated as a physical correction.
"""
import numpy as np

from .labeling_density import density_nyquist


def tube_density(entry, group_ids=None, *, segment_length_nm=200., min_events=30,
                 cutoff_nm=None, sensitivity_limit=.25, radial_policy="Model warning", n_rejected=None):
    """Return summary, segment rows and a mask matching the density numerator.

    Input rows must already be grouped (one centroid per event group).
    Flat end planes delimit the sampled length. Endpoint projections are
    excluded from both the count and radius estimation. Interior holes are
    filled by the explicitly assumed solid-tube model. Radii are apparent
    (precision-broadened), so all densities/Nyquist lengths are exploratory.
    """
    if not np.isfinite(segment_length_nm) or segment_length_nm <= 0 or min_events < 12:
        raise ValueError("Segment length must be positive; minimum events must be at least 12")
    if not np.isfinite(sensitivity_limit) or not 0 < sensitivity_limit <= 1:
        raise ValueError("Sensitivity tolerance must be in (0, 1]")
    if radial_policy not in ("Model warning", "Strict"):
        raise ValueError("Unknown radial QC policy")
    if n_rejected is not None and (not np.isfinite(n_rejected) or n_rejected < 0 or int(n_rejected) != n_rejected):
        raise ValueError("Rejected count must be a nonnegative integer or None")
    dimension = entry.points.shape[1]
    if dimension not in (2, 3):
        raise ValueError("Density requires 2D or 3D coordinates")
    measure_key = "Area_um2" if dimension == 2 else "Volume_um3"
    n = len(entry.points)
    ids = np.arange(1, n+1) if group_ids is None else np.asarray(group_ids)
    if ids is not None and (ids.shape != (n,) or not np.isfinite(ids).all()
                           or np.any(ids <= 0) or np.any(ids != np.floor(ids)) or len(np.unique(ids)) != n):
        raise ValueError("Group IDs must be unique positive integers aligned with retained events")
    record = dict(N_grouped=n, density_dimension=f"{dimension}D",
                  N_independent_estimated=None, Volume_um3=None, Area_um2=None,
                  grouped_event_density=None, grouped_event_nyquist_nm=None,
                  density_status="unavailable", density_segment_length_nm=segment_length_nm,
                  density_min_events=min_events, density_sensitivity_limit=sensitivity_limit,
                  density_radial_policy=radial_policy, density_point_influence_limit=sensitivity_limit,
                  density_effective_cutoff_nm=cutoff_nm,
                  density_rejected_count=n_rejected,
                  density_pre_rejection_groups=None if n_rejected is None else n+int(n_rejected),
                  density_model="uniform filled circular tube; local apparent moment radii",
                  density_scope="post-rejection grouped centroids with interior centerline projections",
                  density_units=f"groups/um^{dimension}", nyquist_units="nm",
                  independent_reason="No calibrated blinking/identity model; groups are not independent labels.",
                  density_interpretation="Exploratory sampling estimates, not image resolution. "
                  "Repeated blinking and localization uncertainty bias counts and apparent radii.")
    support = np.zeros(n, bool)
    if entry.nearest.shape != entry.points.shape or entry.curve.shape[1] != dimension:
        record["density_warnings"] = "Inconsistent centerline coordinate dimensions."
        return record, [], support
    if dimension == 2:
        record["density_model"] = "uniform filled planar strip; local apparent moment half-widths"
    length = entry.summary["centerline_length_nm"]
    support = (entry.arc > 1e-6) & (entry.arc < length-1e-6)
    record.update(N_endpoint_groups_excluded=int((~support).sum()),
                  N_grouped_in_geometry=int(support.sum()))
    residual = entry.points-entry.nearest
    # Curve samples are used only for local frames and geometric diagnostics;
    # the input distances and arc coordinates use continuous spline projection.
    steps = np.linalg.norm(np.diff(entry.curve, axis=0), axis=1)
    curve_arc = np.r_[0., np.cumsum(steps)]
    curve_arc *= length/curve_arc[-1]
    tangents = np.gradient(entry.curve, axis=0)
    tangents /= np.linalg.norm(tangents, axis=1)[:, None]

    def segments(count):
        edges = np.linspace(0., length, count+1)
        rows = []
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            mask = support & (entry.arc >= lo) & (entry.arc < hi)
            k = int(mask.sum())
            row = dict(segment=i+1, arc_start_nm=lo, arc_end_nm=hi,
                       length_nm=hi-lo,
                       N_grouped=k,
                       radius_nm=None, Volume_um3=None, Area_um2=None, status="unavailable", warning="")
            if k < min_events or len(np.unique(entry.points[mask], axis=0)) < 12:
                row["warning"] = "Too few events or unique positions in segment"
            else:
                r2 = entry.distances[mask]**2
                radius = float(np.sqrt((3 if dimension == 2 else 2)*np.mean(r2)))
                row["radius_nm"] = radius
                row[measure_key] = (2*radius*(hi-lo)/1e6 if dimension == 2
                                    else np.pi*radius**2*(hi-lo)/1e9)
                # For a uniformly filled disk, r^2 is uniform. A grossly
                # inconsistent moment ratio flags shells, residual outliers,
                # and strongly nonuniform sampling without claiming a fit test.
                moment_ratio = np.mean(r2*r2)/np.mean(r2)**2 if radius > 0 else np.inf
                tangent = np.array([np.interp((lo+hi)/2, curve_arc, tangents[:, a]) for a in range(dimension)])
                tangent /= np.linalg.norm(tangent)
                _, _, basis = np.linalg.svd(tangent.reshape(1, dimension))
                transverse = residual[mask] @ basis[1:].T
                eig = np.linalg.eigvalsh(transverse.T @ transverse/k)
                ratio = eig[-1]/max(eig[0], 1e-12)
                row.update(radial_moment_ratio=float(moment_ratio), transverse_axis_ratio=float(np.sqrt(ratio)))
                # Describe the influence of the farthest centroid without
                # deleting it or treating this diagnostic as a fitted model.
                largest_fraction = float(np.max(r2)/np.sum(r2)) if radius > 0 else None
                area_change = ((k*largest_fraction-1)/(k-1)
                               if largest_fraction is not None and k > 1 else None)
                row.update(radial_moment_lower_limit=1.12, radial_moment_upper_limit=2.5,
                           largest_r2_fraction=largest_fraction,
                           area_change_without_farthest=area_change)
                reasons, advisories = [], []
                if radius <= 0 or not np.isfinite(radius):
                    reasons.append("Degenerate radius")
                if moment_ratio < 1.12 or moment_ratio > 2.5:
                    direction = "below" if moment_ratio < 1.12 else "above"
                    limit = 1.12 if direction == "below" else 2.5
                    detail = ("radial distances are narrowly distributed; inspect centering and cross-section sampling"
                              if direction == "below" else
                              "broad radial tail; inspect localization blur, residual outliers and within-segment radius variation")
                    influence = f"; farthest-point area influence={area_change:.1%}" if area_change is not None else ""
                    message = (f"Radial-moment heuristic: Q={moment_ratio:.4f} {direction} {limit:g}, "
                                   f"N_grouped={k}{influence}; {detail}. "
                                   "This threshold is not a statistical test or proof of a hollow structure")
                    if radial_policy == "Strict" or moment_ratio < 1.12 or not np.isfinite(moment_ratio):
                        reasons.append(message)
                    else:
                        advisories.append(message + "; apparent filled-tube result is model-dependent")
                if area_change is not None and area_change > sensitivity_limit:
                    reasons.append(f"Single-centroid area influence {area_change:.1%} exceeds {sensitivity_limit:.1%}")
                if ratio > 4:
                    reasons.append("Strongly anisotropic cross-section; circular model uncertain")
                if cutoff_nm is not None and radius >= .8*cutoff_nm:
                    details = (f"Radius near rejection cutoff: radius={radius:.2f} nm, "
                               f"effective cutoff={cutoff_nm:.2f} nm, trigger={.8*cutoff_nm:.2f} nm")
                    if n_rejected == 0:
                        advisories.append(details + "; 0 grouped points removed by this rejection step; "
                                          "cutoff proximity alone does not demonstrate censoring")
                    else:
                        count_text = ("rejected count unknown" if n_rejected is None else
                                      f"3D rejection removed {n_rejected}/{n+int(n_rejected)} grouped points "
                                      f"({n_rejected/(n+n_rejected):.1%}) across the ROI")
                        reasons.append(details + "; " + count_text +
                                       "; full cross-section may be censored. Review the cutoff against "
                                       "structure width and localization uncertainty")
                local = (curve_arc >= lo) & (curve_arc <= hi)
                if local.sum() > 2:
                    curvature = np.linalg.norm(np.gradient(tangents[local], curve_arc[local], axis=0), axis=1)
                    if np.max(curvature)*radius >= .5:
                        reasons.append("Bend too tight for local tube approximation")
                row["warning"] = "; ".join(reasons+advisories)
                row["status"] = "unavailable" if reasons else "exploratory_warning" if advisories else "exploratory"
            rows.append(row)
        return rows

    count = max(1, int(np.ceil(length/segment_length_nm)))
    # Avoid unbounded allocation for accidentally tiny lengths.
    if count > 1000:
        record["density_warnings"] = "More than 1000 segments requested; increase segment length."
        return record, [], support
    rows = segments(count)
    failures = [f"Segment {r['segment']}: {r['warning']}" for r in rows if r["status"] == "unavailable"]
    advisories = [f"Segment {r['segment']}: {r['warning']}" for r in rows
                  if r["status"] == "exploratory_warning"]
    alternate = segments(max(1, count//2))
    advisories += [f"Coarse segment {r['segment']}: {r['warning']}" for r in alternate
                   if r["status"] == "exploratory_warning"]
    if not failures and all(r["status"].startswith("exploratory") for r in alternate):
        volume = sum(r[measure_key] for r in rows)
        coarse = sum(r[measure_key] for r in alternate)
        sensitivity = abs(coarse-volume)/volume
        record.update({("density_area_coarse_um2" if dimension == 2 else "density_volume_coarse_um3"): coarse, "density_segment_sensitivity": sensitivity})
        if count == 1:
            failures.append("Only one segment; local-radius sensitivity cannot be assessed")
        elif sensitivity > sensitivity_limit:
            failures.append("Area/volume sensitive to segment length")
        # Detect nonlocal overlap; neighbouring cross-sections necessarily touch.
        sample = np.linspace(0, length, 201)
        xyz = np.column_stack([np.interp(sample, curve_arc, entry.curve[:, a]) for a in range(dimension)])
        centers = [(r["arc_start_nm"]+r["arc_end_nm"])/2 for r in rows]
        radii = np.interp(sample, centers, [r["radius_nm"] for r in rows])
        separation = np.linalg.norm(xyz[:, None]-xyz[None, :], axis=2)
        radius_sum = radii[:, None]+radii[None, :]
        nonlocal_pair = np.abs(sample[:, None]-sample[None, :]) > 2*radius_sum
        if np.any(nonlocal_pair & (separation < radius_sum)):
            failures.append("Nonlocal tube overlap; summed volume would double count")
        if not failures:
            record.update({measure_key: volume, "density_status": "exploratory_warning" if advisories else "exploratory"})
            for basis, number in (("grouped_event", record["N_grouped_in_geometry"]),):
                if number:
                    record[basis+"_density"], record[basis+"_nyquist_nm"] = density_nyquist(number, volume, dimension)
    elif not failures:
        failures.append("Coarser segmentation failed quality checks; sensitivity unavailable")
    for row in rows:
        for basis, key in (("grouped_event", "N_grouped"),):
            value = (None, None)
            if record["density_status"].startswith("exploratory") and row[key]:
                value = density_nyquist(row[key], row[measure_key], dimension)
            row[basis+"_density"], row[basis+"_nyquist_nm"] = value
    record["density_warnings"] = "; ".join(failures+advisories) or "Apparent radii include localization blur; verify uniform filling and full transverse ROI/slab coverage."
    return record, rows, support

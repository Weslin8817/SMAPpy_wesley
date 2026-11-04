"""Current Render ROI centerline and continuous radial distance measurement."""
from copy import copy, deepcopy
from dataclasses import asdict, dataclass
import json

import numpy as np

from smappy.plugins import Context, FormState, Plot, Plugin, Result, Selection, param
from .centerline import CenterlineResult, ROIResults, fit_selected_centerline, result_summary
from smappy.locs import Localizations
from .tube_density import tube_density


@dataclass
class Settings:
    rejection: str = param("None", label="Noise rejection",
                           choices=(("None", "None"), ("Manual maximum distance", "Manual maximum distance"),
                                    ("Auto conservative", "3D auto conservative")),
                           help="Choose before analysis. Manual and Auto are independent alternatives; "
                                "Manual and Auto select in XYZ (requires z), before either analysis. "
                                "No rejection iteration.")
    maximum_distance_nm: float = param(50., label="Maximum distance", unit="nm", min=.1,
                                       help="Manual mode only: radius from the initial centerline, "
                                            "not diameter. A 100 nm diameter corresponds to 50 nm. "
                                            "Final refitting can move distances slightly across this cutoff.")
    rejection_bins: int = param(16, label="Denoising bins", min=4, max=64,
                                help="Maximum centerline bins for the shared XYZ noise selection.")
    rejection_smoothing: float = param(1., label="Denoising smoothing", min=.01, max=100.,
                                      help="Curvature penalty for the shared XYZ noise selection.")
    radial: bool = param(True, label="Enable Radial Distance",
                         help="Calculate radial distances and percentiles in the dimension selected below.")
    dimension: str = param("2D", label="Dimension", choices=("2D", "3D"),
                           help="Radial-distance analysis: XY or XYZ, after shared 3D rejection.")
    bins: int = param(16, label="Centerline bins", min=4, max=64,
                      help="Maximum longitudinal bins; automatically reduced for small ROIs.")
    smoothing: float = param(1., label="Smoothing", min=.01, max=100.,
                             help="Curvature penalty: larger values give a smoother centerline.")
    density: bool = param(True, label="Enable Labeling Density",
                          help="Grouped-event density and Nyquist after shared 3D noise rejection.")
    density_dimension: str = param("3D", label="Dimension", choices=("2D", "3D"),
                                  help="2D: apparent filled-strip area (groups/um²). 3D: filled-tube volume (groups/um³).")
    density_bins: int = param(16, label="Centerline bins", min=4, max=64,
                              help="Maximum longitudinal bins for the density centerline, independent of radial analysis.")
    density_smoothing: float = param(1., label="Smoothing", min=.01, max=100.,
                                    help="Curvature penalty for the density centerline.")
    segment_length_nm: float = param(200., label="Segment length", unit="nm", min=1.,
                                    help="Maximum segment length. Each segment has its own apparent strip half-width (2D) or tube radius (3D); no precision correction.")
    min_segment_events: int = param(30, label="Min. groups", min=12,
                                    help="Minimum grouped centroids per segment. Sparse segments invalidate total area/volume.")
    density_sensitivity: float = param(.25, label="Sensitivity", min=.01, max=1.,
                                       help="Maximum relative area/volume change across segment scales or squared transverse-size change when omitting the farthest centroid (diagnostic only).")
    radial_policy: str = param("Model warning", label="Radial QC", choices=("Model warning", "Strict"),
                               help="Model warning: a broad radial tail warns if influence/geometry checks pass. Strict: Q above 2.5 also blocks density. Neither proves physical model validity.")


def roi_coordinates(roi, axes=None, locs=None):
    """Bounding-box center of the captured Render ROI, in physical XY nm.

    Non-spatial Render axes cannot supply a physical XY ROI center.
    """
    x0, y0, x1, y1 = roi.bounds
    center = ((x0+x1)/2., (y0+y1)/2.)
    names = ("x_nm", "y_nm") if axes is None else axes.names(locs)
    scales = (1., 1.) if axes is None else (axes.x_scale, axes.y_scale)
    values = {name: float(value*scale) for name, value, scale in zip(names, center, scales)}
    return {"roi_center_x_nm": values.get("x_nm", ""),
            "roi_center_y_nm": values.get("y_nm", "")}


def display_statistics(entry):
    """R labels denote radial-distance percentiles in nm."""
    return result_summary(entry)


def plots(entry):
    def fit(figure):
        dim = entry.points.shape[1]
        pairs = [(0, 1)] if dim == 2 else [(0, 1), (0, 2), (1, 2)]
        axes = figure.subplots(1, len(pairs), squeeze=False)[0]
        for ax, (a, b) in zip(axes, pairs):
            ax.scatter(entry.points[:, a], entry.points[:, b], s=2, alpha=.3)
            ax.plot(entry.curve[:, a], entry.curve[:, b], color="crimson", lw=2)
            ax.set(xlabel=f"{'xyz'[a]} (nm)", ylabel=f"{'xyz'[b]} (nm)")
            ax.set_aspect("equal", adjustable="datalim")
    def distribution(ax):
        if not len(entry.distances):
            ax.set_title("Distance analysis unavailable; see result text")
            return
        values = display_statistics(entry)
        ax.hist(entry.distances, bins="auto", color="steelblue", alpha=.8)
        for name, color in [("R30", "green"), ("R50", "black"),
                            ("R70", "orange"), ("R90", "crimson"),
                            ("R95", "purple"), ("R100", "gray")]:
            ax.axvline(values[f"{name}_nm"], label=name, color=color)
        ax.set(xlabel="Radial distance (nm)", ylabel="Localizations")
        ax.legend()
    # Plot with panels=1 expects an axis; the fitting plot wants a figure.
    def xy(ax):
        ax.scatter(entry.points[:, 0], entry.points[:, 1], s=2, alpha=.3)
        ax.plot(entry.curve[:, 0], entry.curve[:, 1], color="crimson", lw=2)
        ax.set(xlabel="x (nm)", ylabel="y (nm)")
        ax.set_aspect("equal", adjustable="datalim")
    result = {"Centerline": Plot(xy if entry.points.shape[1] == 2 else fit,
                                 panels=1 if entry.points.shape[1] == 2 else 3),
              "Radial distribution": distribution}
    if not entry.summary.get("radial_enabled", True):
        result = {}
    if getattr(entry, "density_segments", []):
        def radii(ax):
            rows = entry.density_segments
            for status, color, label in (("exploratory", "steelblue", "Passed segment checks"),
                                         ("exploratory_warning", "darkorange", "Model warning"),
                                         ("unavailable", "crimson", "Failed segment checks")):
                selected = [r for r in rows if r["radius_nm"] is not None and
                            r["status"] == status]
                if selected:
                    ax.errorbar([(r["arc_start_nm"]+r["arc_end_nm"])/2 for r in selected],
                                [r["radius_nm"] for r in selected],
                                xerr=[r["length_nm"]/2 for r in selected],
                                fmt="o", color=color, label=label)
            ax.set(xlabel="Centerline arc (nm)", ylabel="Apparent filled-tube radius (nm)",
                   title=f"{entry.summary['density_dimension']} grouped density: {entry.summary['density_status']}")
            ax.set_ylabel("Apparent strip half-width (nm)" if entry.summary["density_dimension"] == "2D"
                          else "Apparent filled-tube radius (nm)")
            if ax.has_data():
                ax.legend()
        result["Segment radii"] = radii
    return result


class CenterlineRadial(Plugin):
    description = "3D noise rejection → Radial Distance and Labeling Density. Each analysis has its own dimension."
    favorite = False
    version = "9"
    changes_locs = False
    # Group on a scratch session: never attach grouping columns to the input.
    grouping = "ungrouped"
    path = "Analysis/Measure/Radial distance and density"
    name = "Radial distance and density"
    Settings = Settings

    def __init__(self):
        self.results = ROIResults()

    def active(self, settings):
        return {"maximum_distance_nm": settings.rejection == "Manual maximum distance",
                "rejection_bins": settings.rejection != "None",
                "rejection_smoothing": settings.rejection != "None",
                **{key: settings.radial for key in ("dimension", "bins", "smoothing")},
                **{key: settings.density for key in
                   ("density_dimension", "density_bins", "density_smoothing", "segment_length_nm", "min_segment_events", "density_sensitivity", "radial_policy")}}

    def prepare_context(self, ctx, settings):
        if ctx.session is None:
            return ctx
        source = ctx.session
        index = source.render_layer_index
        if source.roi is None:
            raise ValueError("Draw a current ROI in Render before analysis")
        if not 0 <= index < len(source.layers) or source.layers[index].is_image:
            raise ValueError("Select a localization layer in Render")
        session = source.scratch()
        locs, selection = session.table(index, grouping="grouped")
        # Copy only selected rows on the GUI thread: filter, grouping, ROI and
        # original table row indices cannot change underneath the worker.
        selected = Localizations({k: np.array(v[selection.mask], copy=True)
                                  for k, v in locs.columns.items()}, deepcopy(locs.metadata))
        snapshot = Context(locs=selected, selection=Selection(
            np.ones(len(selected), bool), layer=index, roi=deepcopy(session.roi), name=selection.name))
        snapshot.centerline_indices = selection.indices.copy()
        layer = session.layers[index]
        snapshot.centerline_group_ids = None
        if settings.density and "group_id" in selected:
            ids = np.asarray(selected["group_id"])
            if np.isfinite(ids).all() and np.all(ids > 0) and np.all(ids == np.floor(ids)):
                snapshot.centerline_group_ids = ids.copy()
        snapshot.centerline_provenance = dict(source_roi=repr(session.roi), layer=index,
            layer_name=session.layers[index].name, grouped=True,
            selection=selection.name, source_path=str(session.path or ""),
            filter_ranges=json.dumps(layer.ranges, default=str),
            slab=repr(session.slab) if session.selects_slab else "")
        if settings.density:
            snapshot.centerline_provenance.update(grouped=True,
                render_grouped=layer.grouped, grouping_parameters=json.dumps(asdict(session.group_settings)),
                grouping_reason="One centroid per grouped-table row; grouped coordinates and precision used for all calculations.")
        snapshot.centerline_provenance.update(roi_coordinates(
            snapshot.selection.roi, session.layers[index].state.settings.axes, locs))
        return snapshot

    def run(self, ctx, settings):
        ctx = self.prepare_context(ctx, settings)
        if ctx.selection.roi is None:
            raise ValueError("Draw a current ROI in Render before analysis")
        if not (settings.radial or settings.density):
            raise ValueError("Enable at least one analysis module")
        for enabled, dimension in ((settings.radial, settings.dimension),
                                   (settings.density, settings.density_dimension)):
            if enabled and dimension not in ("2D", "3D"):
                raise ValueError("Dimension must be 2D or 3D")
        locs = ctx.selection.apply(ctx.locs)
        rejection = ("Auto conservative" if settings.rejection == "3D auto conservative"
                     else settings.rejection)
        rejection_label = "3D auto conservative" if rejection == "Auto conservative" else rejection

        def coordinates(dimension):
            if dimension not in ("2D", "3D"):
                raise ValueError("Dimension must be 2D or 3D")
            axes = "xy" if dimension == "2D" else "xyz"
            for axis in axes:
                if f"{axis}_nm" not in locs:
                    raise ValueError(f"{dimension} requires {axis}_nm coordinates in nm")
            points = np.column_stack([locs[f"{a}_nm"] for a in axes])
            if dimension == "3D" and not np.isfinite(points[:, 2]).any():
                raise ValueError("3D requires valid finite z_nm data")
            precision = np.column_stack([locs.columns.get(f"{a}_err_nm",
                locs.columns.get("xy_err_nm", np.full(len(locs), np.nan))
                if a in "xy" else np.full(len(locs), np.nan)) for a in axes])
            return points, precision

        # One XYZ selection precedes both analyses; projection never reruns rejection.
        retained = np.arange(len(locs))
        cutoff, excluded = None, 0
        rejection_shift = 0.
        noise_nonfinite = None
        fits = {}
        if rejection != "None":
            xyz, errors = coordinates("3D")
            denoised, cutoff, excluded = fit_selected_centerline(
                xyz, errors, bins=settings.rejection_bins,
                smoothing=settings.rejection_smoothing, rejection=rejection,
                maximum_distance_nm=settings.maximum_distance_nm)
            retained = denoised.indices
            noise_nonfinite = denoised.summary["excluded_nonfinite"]
            rejection_shift = denoised.summary.get("rejection_refit_shift_nm", 0.)
            fits[("3D", settings.rejection_bins, settings.rejection_smoothing)] = denoised

        def analyze(dimension, bins, smoothing):
            key = (dimension, bins, smoothing)
            if key not in fits:
                points, precision = coordinates(dimension)
                fitted, _, _ = fit_selected_centerline(
                    points[retained], precision[retained], bins=bins,
                    smoothing=smoothing, rejection="None")
                fitted.indices = retained[fitted.indices]
                fitted.summary["excluded_nonfinite"] = int(
                    (~np.isfinite(points).all(axis=1)).sum())
                fits[key] = fitted
            # Identical settings can share geometry, but never mutable summaries.
            fitted = copy(fits[key])
            fitted.summary = dict(fitted.summary)
            return fitted

        def empty_entry(status, warning=""):
            return CenterlineResult(np.empty((0, 2)), np.empty(0, dtype=int),
                np.empty((0, 2)), np.empty(0), np.empty(0), np.empty((0, 2)),
                dict(dimension=settings.dimension if settings.radial else "",
                     N=0, distance_status=status, distance_warning=warning,
                     excluded_nonfinite=0, precision_replacements=None))

        entry = empty_entry("disabled")
        radial_error = None
        if settings.radial:
            try:
                entry = analyze(settings.dimension, settings.bins, settings.smoothing)
            except ValueError as error:
                if not settings.density:
                    raise
                radial_error = error
                entry = empty_entry("unavailable", str(error))
        entry.summary.update(noise_rejection=rejection_label, noise_cutoff_nm=cutoff,
                             N_before_rejection=len(locs), N_noise_excluded=excluded,
                             N_after_rejection=len(retained), noise_excluded_nonfinite=noise_nonfinite,
                             noise_bins=settings.rejection_bins, noise_smoothing=settings.rejection_smoothing,
                             bins=settings.bins, smoothing=settings.smoothing,
                             radial_enabled=settings.radial, density_enabled=settings.density,
                             plugin_version=self.version)
        if settings.density:
            try:
                density_entry = analyze(settings.density_dimension,
                                        settings.density_bins, settings.density_smoothing)
                density_cutoff, density_excluded = cutoff, excluded
                group_ids = getattr(ctx, "centerline_group_ids", None)
                record, segments, support = tube_density(density_entry,
                    None if group_ids is None else group_ids[density_entry.indices],
                    segment_length_nm=settings.segment_length_nm,
                    min_events=settings.min_segment_events,
                    cutoff_nm=None if density_cutoff is None else density_cutoff-rejection_shift,
                    sensitivity_limit=settings.density_sensitivity, radial_policy=settings.radial_policy,
                    n_rejected=density_excluded)
                record.update(density_length_nm=density_entry.summary["centerline_length_nm"],
                              density_noise_cutoff_nm=density_cutoff, density_noise_excluded=density_excluded,
                              density_excluded_nonfinite=density_entry.summary["excluded_nonfinite"])
                density_export = copy(density_entry)
                density_export.density_support = support
                density_export.indices = getattr(ctx, "centerline_indices", ctx.selection.indices)[density_entry.indices]
                entry.density_support = np.isin(entry.indices, density_entry.indices[support])
                entry.density_entry = density_export
            except ValueError as error:
                if radial_error is not None:
                    raise radial_error
                record = dict(density_dimension=settings.density_dimension, density_status="unavailable", density_warnings=str(error),
                              Volume_um3=None, Area_um2=None, N_grouped=None, N_grouped_in_geometry=None,
                              grouped_event_density=None, grouped_event_nyquist_nm=None)
                segments = []
            entry.summary.update(record)
            entry.summary.update(density_bins=settings.density_bins,
                                 density_smoothing=settings.density_smoothing,
                                 density_segment_length_nm=settings.segment_length_nm,
                                 density_min_events=settings.min_segment_events,
                                 density_sensitivity_limit=settings.density_sensitivity,
                                 density_point_influence_limit=settings.density_sensitivity,
                                 density_radial_policy=settings.radial_policy)
            entry.density_segments = segments
        source_indices = getattr(ctx, "centerline_indices", ctx.selection.indices)
        entry.indices = source_indices[entry.indices]
        provenance = getattr(ctx, "centerline_provenance", dict(source_roi=repr(ctx.selection.roi),
                             layer=ctx.selection.layer, selection=ctx.selection.name))
        if "roi_center_x_nm" not in provenance:
            provenance = dict(provenance, **roi_coordinates(ctx.selection.roi))
        identifier = self.results.add(entry, **provenance)
        selection_text = ("" if cutoff is None else
                          f" {rejection_label}: initial-distance cutoff {cutoff:.2f} nm; "
                          f"excluded {excluded} points (one pass).")
        median_text = f"{np.median(entry.distances):.2f} nm" if len(entry.distances) else "N/A"
        radial_text = (f"{identifier}: {len(entry.points)} localizations; "
                           f"median radial distance {median_text}. "
                           f"Excluded nonfinite: {entry.summary['excluded_nonfinite']}; "
                           f"precision fallbacks: {entry.summary['precision_replacements']}."
                       if settings.radial else f"{identifier}: Radial Distance disabled.")
        return Result(text=radial_text
                           + selection_text + ("\n" + "\n".join(
                               f"{k}: {'N/A' if v is None else v}" for k, v in entry.summary.items())
                               if settings.density else ""),
                      data={"centerline_roi": identifier, **entry.summary}, plots=plots(entry), settings=settings)

    def describe(self, ctx, settings):
        return FormState(active=self.active(settings),
                         shown={"maximum_distance_nm": settings.rejection == "Manual maximum distance"})

    def create_result_widget(self, panel):
        from .gui import CenterlineResultsWidget, configure_panel
        configure_panel(panel)
        return CenterlineResultsWidget(panel)

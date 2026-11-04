---
version: "9"
covers: [smappy_centerline.centerline.fit_centerline, smappy_centerline.centerline.fit_selected_centerline, smappy_centerline.centerline.manual_distance_mask, smappy_centerline.centerline.conservative_distance_mask, smappy_centerline.centerline.project_curve, smappy_centerline.centerline.sanitize_precision, smappy_centerline.centerline.ROIResults, smappy_centerline.tube_density.tube_density]
---

## What it measures

Draw one ROI in Render and select its localization layer. Open **Centerline &
Radial Distance** from Plugins / Analysis / Measure (or add it to a tab).
Configure the three sections before **Analyze/Add current ROI**:
1. **3D noise rejection**: select None, Manual maximum distance, or 3D auto conservative, with independent denoising bins and smoothing.
2. **Radial Distance**: enabled by default; choose its own 2D/3D Dimension, bins and smoothing.
3. **Labeling Density**: enabled by default; choose its own 2D/3D Dimension, bins, smoothing and density parameters.

Disabled modules retain their settings. At least one module must be enabled.
The default dimensions remain 2D for radial distance and 3D for density. Each
successful run adds an independent snapshot with a new ROI identifier.
All calculations use the current layer's grouped localization table, its filter,
ROI and active slab restriction, regardless of the Render grouping display toggle.
Missing or stale groups are rebuilt on an isolated scratch session using the
current Session grouping settings, as required by the newer shared-table API.
The original localization table, filters and Render state remain unchanged.
Each grouped row contributes once; its coordinates and combined precision are used.
Image layers are refused. Coordinates must be available in nm.

The table, histogram and summary CSV show **R30, R50, R70, R90, R95 and R100**,
the corresponding percentiles of radial distance in nm. R100 is the maximum
retained radial distance, not a fitted structural boundary. These remain
radial distances: no values are multiplied by two. Mean, median, RMS,
SD, IQR and manual cutoff retain their original radial-distance values.
Summary CSV percentiles use the same calculation as the table; R80 is replaced
by the current set above. Per-localization distance CSV values are unchanged.

Results appear in three vertically stacked tables: ROI information and shared
noise rejection, Radial Distance, and Labeling Density. Every table retains
the same ROI identifier and row order, including rows for disabled modules.
Selecting a row in any table selects the same ROI in the other two tables;
delete and clear operate on the shared ROI record. Sorting is disabled to
preserve correspondence. Multi-row selection is supported.

The ROI table includes bounding-box center X/Y in nm, input groups, invalid
XYZ coordinates, noise removals, retained groups and the actual cutoff.
With rejection set to None, invalid XYZ is shown as a dash: coordinate validity
is handled separately by each analysis dimension. Retained groups then means
the unchanged input count, before analysis-specific finite-coordinate checks.

Radial results show status, dimension, count, length, mean and R30/R50/R70/R90/
R95/R100. R50 already represents the median. Density results show status,
dimension, retained and geometry-domain groups, area or volume with its units,
group density with its units, and Nyquist. Disabled modules keep their ROI row
with blank values. Warnings and unavailable results have explicit text and
colored status cells. Full statistics (including RMS, SD and IQR), diagnostic
messages and provenance are available in the selected ROI details beneath
the tables. Display columns do not restrict the summary CSV: it still exports
all summary fields for each ROI in one row.

Select a result row to inspect the fit and radial histogram. 3D fits show
XY, XZ and YZ projections. Delete selected and Clear all affect only these
measurements. They never change localizations, ROI geometry or rendering.
Results live in this plugin instance until it is discarded; they are not
automatically saved into localization HDF5 or GUI settings files.

## Filled-tube density and Nyquist

Enable **Labeling Density** and select its independent **Dimension**. This replaces the separate
Labeling Density plugin. Density uses the selected Render layer's grouped
filter, ROI and slab. Noise rejection runs first; retained grouped localizations
supply both geometry and counts. Group IDs refer to the same grouped rows.
Only grouped-event density and Nyquist are reported, not raw-event density.
Group counts are not independent molecule counts.

Equal arc-length segments (maximum length **200 nm**, minimum **30 events**
and 12 unique positions each) independently estimate apparent radii. For a
uniformly filled circular disk, radius squared is twice mean squared radial
distance. Segment volume is pi times radius squared times arc length, converted
from nm cubed to micrometres cubed. This fills internal unsampled spaces by
model assumption. Radii include localization blur; precision fallback values
are never used for physical correction. This is not a shell model or a
LocMoFit likelihood fit. LocMoFit's curve-plus-tube design is described in its
[model library](https://locmofit.readthedocs.io/en/stable/LocMoFit.modelLibrary.html).

Flat end planes delimit the sampled centerline length. Events projecting onto
endpoints are excluded from density and counted separately; their distances
contain axial offsets. No missing physical ends are extrapolated. Per-point
CSV records membership in this interior domain. Total and segment group
counts use one centroid per group. Each centroid belongs to one segment,
so segment group counts sum to the total count in the volume domain.

Approximately doubling segment length tests volume sensitivity (default
tolerance **0.25**). Sparse segments, narrow radial distributions, strong
anisotropy, tight bends, suspected overlap, or radii near a rejection cutoff
invalidate the total. One segment cannot establish local-radius sensitivity.
Candidate segment measurements remain exportable with quality notes. These
diagnostics cannot prove uniform filling or detect every clipped structure.
Ensure the ROI/slab includes full cross-sections. Branches, folds, hollow
fibers and anisotropic localization blur require particular care.

Diagnostic limits: radial fourth/second-moment ratio outside 1.12–2.5,
transverse second-moment eigenvalue ratio above 4, curvature times radius at
least 0.5, and radius at least 80% of the rejection cutoff after subtracting
maximum centerline refit displacement. These heuristics are not confidence
intervals. The Segment radii plot displays local apparent radii and failed
segments; its title reports total density status.

Radial-moment failures report Q, whether it is below/above the limit, group
count, and the relative change in estimated cross-section area if the
farthest centroid were omitted. That last value is diagnostic only: no point
is actually removed. Q is mean(r^4)/mean(r^2)^2; the ideal noiseless uniform
disk gives 4/3. A high Q indicates a broad radial tail, not a hollow structure.
Localization blur, residual outliers and radius changes within a segment can
all affect it. These fixed cutoffs do not account for sampling uncertainty;
a marginal exceedance is not proof that the biological model is invalid.

Only grouped-event density and Nyquist are calculated. Density is retained
interior groups per micrometre cubed; Nyquist
is 2000 divided by cube root of density, in nm. Both are exploratory sampling
estimates affected by blinking, not measured image resolution. Independent
labels remain unknown without a calibrated identity/blinking model;
effective labeling density/Nyquist are not calculated. Results retain rejection settings, counts,
segmentation, sensitivity and grouping provenance. Export additionally writes
**centerline_density_segments.csv** with each segment's radius, volume, counts,
density, Nyquist and quality notes when segment results exist.

## Optional noise rejection

Choose one independent mode before running:

- **None** preserves all finite input coordinates, as before.
- **Manual maximum distance** keeps points at or below the specified radius
  from an initial centerline. The default is **50 nm**, corresponding to a
  **100 nm diameter**. It is a selection boundary, not a fitted tube radius.
- **3D auto conservative** always uses XYZ distances, even with Dimension set
  to 2D. It requires valid Z coordinates and uses median(distance) + 4 × max(1 nm, 1.4826 × MAD)
  from the initial centerline distances. MAD is the median absolute deviation
  from their median. It does not use the manual radius or assume a 100 nm
  diameter; it is intended to exclude conspicuously distant points.

Both rejection modes use one initial precision-weighted centerline fit,
one selection, and one final fit on retained points if any were excluded.
There is no repeated rejection/refitting loop. The final distances are
recomputed in Euclidean space, without clipping. Because the final curve can
move, some final distances can exceed the initial selection cutoff. The log
reports that cutoff and the number excluded. Fewer than 48 retained points
is an error and does not add a result row.

Automatic rejection uses 3D centerline distances and Z precision for both
dimensions. With a 2D centerline, the retained XYZ point indices are then
fitted in XY without another rejection pass. Each module fits the retained set with its own bins and smoothing. Missing valid Z is an error, not a fallback to 2D.
Manual rejection also selects once in XYZ before both modules.
Automatic rejection is a heuristic, not a biological boundary,
and cannot guarantee separation when signal and noise overlap.

## Estimation

A robust, precision-weighted PCA defines longitudinal ordering. Equal-width
axial bins provide robust centers using coordinate-specific inverse variance
weights with a common Huber residual factor. A cubic B-spline with at most
12 coefficients per coordinate fits the centers, weighted by their summed
precision information. Integrated squared curvature is penalized; increasing
Smoothing suppresses curvature. Bins controls the maximum number of local
centers and is reduced to allow at least 12 input points per bin on average.
Bins with fewer than six points are omitted; at least four must survive.
Neither knots nor a fitted radius are assigned to individual localizations.

The fit uses x_err_nm and y_err_nm, falling back to xy_err_nm when an axis
column is absent. In 3D it also uses z_err_nm independently. Invalid errors
(nonfinite or outside 0.1–10000 nm) use the valid axis median. With none,
XY uses 20 nm and Z uses 60 nm. A 1 nm floor limits extreme precision weight.
The number of replaced uncertainty values is reported. These defaults are
assumptions when calibration is absent, not inferred measurement precisions.
Isotropic supplied uncertainties give rotationally equivariant estimation;
anisotropic uncertainties refer to the acquisition coordinate axes.

Distances are **Euclidean spatial distances**, with no precision subtraction,
deconvolution, FWHM conversion or hollow-tube assumption. Every retained finite
localization contributes to all statistics, including points downweighted
while fitting. With Noise rejection set to None this includes every finite
input point. Nonfinite coordinates are excluded and counted. SD uses the
population convention (ddof=0); IQR is P75−P25; R30/R50/R70/R90/R95/R100 are empirical
percentiles and RMS is sqrt(mean(distance²)).

## Continuous projection and endpoints

For each cubic polynomial span, all real roots of the derivative of squared
distance are evaluated alongside its endpoints. The best candidate across
all spans gives the nearest point, not merely the nearest sampled position.
The arc coordinate is integrated from spline speed on a dense grid; it starts
at the deterministically oriented first endpoint. Endpoints use the 0.5 and
99.5 percentiles of robust axial support, so points outside this finite curve
include an endcap contribution in their distances. The curve is sampled for
display and tube-quality diagnostics; distances use continuous projection.
Root solving scales with localization count and curve spans.

## Export

Export CSV writes summary and distance files, plus a segment file when density
segments exist, into a chosen directory. The summary contains
one row per retained ROI, all statistics, fitting settings, actual populated
bins/coefficient count, precision fallback count, excluded coordinate count,
and source ROI/layer/grouping provenance. The distances file starts with
ROI,distance_nm and contains **all retained ROIs and their selected points**, followed by the source
table row index, dimension, arc coordinate and nearest XYZ coordinates (Z
blank in 2D). For grouped data the index refers to the grouped table. IDs
are not reused after deletion or clearing. Export is the persistent record.

## Applicability

Designed for a single elongated, open, unbranched structure with a continuous
filled radial distribution. At least 48 finite points, 1 nm extent and PCA
variance ratio of 2 are required. PCA ordering cannot resolve loops, strong
foldbacks, intersections, branches or multiple nearby filaments. Nonuniform
labeling and asymmetric biological occupancy can shift the estimated axis;
precision weighting cannot distinguish those from physical asymmetry.
Validate smoothing, endpoint bias and density sensitivity against synthetic
controls and experimental references before biological interpretation.


### Independent dimensions (version 9)

Both active rejection modes select once in XYZ. Each enabled module fits that
retained set in its chosen dimension without repeating rejection. Changing
radial settings does not change density. With None selected, rejection is
skipped; 2D analysis can use XY-only data. Active rejection requires valid Z.

3D density uses the filled circular tube model above. 2D density uses a
uniform filled planar strip: apparent half-width is sqrt(3 * mean(r²)),
segment area is 2 * half-width * length, and density is groups/µm².
The 2D Nyquist estimate is 2000 / sqrt(density) nm. This is a planar model,
not a reconstruction of volume or a correction for projecting a 3D tube.
The same heuristic stability gates apply; transverse anisotropy has only
one component in 2D. Cross-section influence is conservatively evaluated
on the squared half-width. The ideal uniform strip gives Q = 9/5.

Summary and segment exports use Area_um2 for 2D and Volume_um3 for 3D;
density_units records the appropriate units. Density distances are written
to centerline_density_2d_distances.csv or centerline_density_3d_distances.csv,
with a blank nearest Z for 2D. Disabled radial analysis produces no radial
plots or distances. A failed analysis is reported as unavailable while the
other module can still succeed.
### Radial QC policy (version 7)

**Radial QC** defaults to **Model warning**. Q above 2.5 identifies a broad
radial tail; it no longer alone rejects an otherwise stable apparent volume.
Such results are explicitly **exploratory_warning**, with orange segments in
the radius plot and model warnings retained in CSV. These remain conditional
filled-tube estimates: passing stability checks does not validate uniform
filling or remove localization blur. All retained centroids still contribute.

**Strict** retains the high-Q rejection. In both modes, Q below 1.12 still
rejects the filled-disk estimate. Omitting the farthest centroid hypothetically
must change each segment's area by no more than **Sensitivity** (default 25%).
This new influence guard applies at both segment scales, alongside the
existing volume sensitivity, sparse-segment, clipping, anisotropy, curvature
and overlap checks. It is not a confidence interval and does not establish
robustness to multiple outliers. The omitted-point calculation is diagnostic;
no localization is actually removed. Policy and influence limit are exported.


### Rejection-cutoff diagnostics (version 8)

The cutoff is a radial distance, not a diameter. The 80% proximity rule uses
cutoff minus the maximum centerline refit displacement. This is conservative;
a displacement at one location reduces the effective cutoff for all segments.
Proximity is a heuristic, not proof that a particular segment lost points.

When 3D rejection removed zero finite points, cutoff proximity alone now gives
an advisory instead of blocking density: this rejection step did not truncate
the observed point set. This does not establish complete labeling, nor rule
out earlier filter/ROI/slab clipping. If any points were removed (or the count
is unknown), the existing blocking rule remains. Messages show radius,
effective cutoff, trigger and ROI-wide 3D removed count/fraction. Counts and
cutoff are saved in the summary. Radial QC Model warning does not override
this protection. Clean simulated structures can lose valid peripheral points
to a manual cutoff; localization blur can also extend coordinates beyond a
physical radius. Remaining points alone cannot identify a fully censored outer
boundary without additional assumptions or ground truth.

## CSV export format 2

The summary remains one row per ROI, in ROI insertion order. Its header is fixed
across runs, module combinations and failure states. Columns are arranged as:
ROI/shared noise rejection; radial core results; density core results;
additional statistics/diagnostics; analysis parameters; warnings/provenance.
The UI continues to use its existing three tables.

Field names use roi_, noise_, radial_ and density_ prefixes (ROI remains the
shared identifier). For example, N becomes radial_groups, R50_nm becomes
radial_r50_nm, Area_um2 becomes density_area_um2 and Volume_um3 becomes
density_volume_um3. Existing CSV-processing scripts must use the new names.
export_format_version is 2. centerline_summary_fields.csv lists every summary
column in order, with its Chinese description, unit, applicability, empty-value
meaning and legacy field name where renamed.

Module statuses are disabled, unavailable, complete or warning. Density
complete/warning results remain exploratory: density_interpretation preserves
this qualification. Disabled modules have blank values. Failed radial fitting
does not export the placeholder count 0 as a measured result. Empty cells mean
not applicable, unavailable, disabled or not recorded; inspect status and
warning fields. Actual zeros remain zero. Numeric precision is not rounded to
the UI display precision.

Only R50 is exported for the median; duplicate median and shared rejection
aliases are omitted. Actual selection cutoff and effective diagnostic cutoff
remain distinct. Area and volume have separate numeric columns, only populated
for the corresponding dimension. Detailed settings now include independent
density bins and smoothing, even when an enabled density analysis fails.
Shared denoising parameters are blank when rejection is None. Unknown extension
metadata is retained in additional_metadata_json, without changing the header.

The point-distance CSV formats remain unchanged. Density segment CSV additionally
includes density_dimension and density_units; missing numeric values are blank.
The field dictionary is exported alongside the summary, and is included in
the overwrite confirmation.

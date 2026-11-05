# Nikon ND2 input

Completed Nikon ND2 files can be selected in **Localize > source > file** and
used directly for preview, Gaussian or spline fitting, and HDF5 output.
Batch folder discovery includes ND2, including mixed-case extensions; explicit
file inputs work too. Each file is a separate acquisition, even beside a TIFF
series or an NDTiff dataset. The default output is `<name>_locs.hdf5`.

## Installation

The optional reader is imported only when an ND2 file is opened. Other image
formats work without it; attempting ND2 without it gives an installation hint.
Install into the same environment that runs SMAPpy:

```powershell
.\.venv\python.exe -m pip install "nd2>=0.10" --index-url https://pypi.tuna.tsinghua.edu.cn/simple
```

Package installations can also request `smappy-smlm[nd2]`. This checkout requires
Python >=3.12, as declared by the upstream project.

## Reading and metadata

`src/smappy/io/nd2.py` ports the reader from `G:/SMAPpy_wes`, using the public
`nd2.ND2File.read_frame()` API. Opening reads dimensions and metadata, not the
whole image array. `frame(index)` reads one frame; `frames(chunk, start, stop)`
allocates one block at a time. Stop is exclusive and indices are zero-based.
The fitting pipeline's existing bounded read-ahead may retain additional blocks.
Pixels keep their original dtype, values, orientation and temporal order.

Each random read and each block iterator owns its own context-managed reader.
Preview and fitting threads never share an ND2 handle. Callers abandoning an
iterator early should close it (for example, with `contextlib.closing`). Read
errors propagate and close the handle; invalid frames are not silently cast.

Only YX/TYX, single-channel, single-Z, single-position grayscale data is accepted.
Non-spatial singleton dimensions are allowed; multiple channels, Z planes,
positions, other non-time dimensions and RGB are rejected, never flattened.

Dimensions, significant bit depth, channel and objective names are retained in
the acquisition summary and written to HDF5 by the existing pipeline. Pixel
sizes are used in micrometres only when **both X and Y are explicitly calibrated**
and finite and positive. Unequal X/Y values are preserved. Missing or uncalibrated
values stay unset. Missing descriptive names do not discard valid calibration.

Nikon gain fields are not interpreted as conversion, offset or EM gain. The
reader supplies only validated pixel sizes through `camera_values`; the newer
camera database, microscope resolution and user override precedence remain in
place. Choose a camera preset or enter camera parameters as usual. If pixel size
is unknown, the newer fitter's existing assumed-scale warning still applies.

## Limits

- Use completed files only. Live/growing ND2 is unsupported. The source's
  `watch()`, `open_growing_stack()` and `watch_stack()` reject ND2 immediately,
  including a path whose file has not appeared yet. The adapter does not certify
  acquisition completion; stop acquisition before opening it.
- No conversion to TIFF, rescaling, rotation, dimension merging, or timestamp
  interpolation is performed.
- Per-frame exposure, timestamps and stage coordinates are not mapped to frame
  tags, matching the old reader.
- This adds acquisition input, not ND2 support to the separate rendered-image
  or overlay loader.
- Spline fitting requires a suitable PSF calibration; loading ND2 does not
  provide camera calibration or a PSF model.

## Integration boundaries for future upstream updates

There is no image-source plugin registry in this version. Integration therefore
uses small changes to existing entry points, without runtime patching:

| File | Change / interface to review after updates |
| --- | --- |
| `src/smappy/io/nd2.py` | Independent `ImageSource` subclass; dimensions, metadata, reads, lifetime and Live policy |
| `src/smappy/io/tiff.py` | Dispatch before directory detectors; `is_image_stack()`; normalised `camera_values` in `resolve_camera()` |
| `src/smappy/io/watch.py` | ND2 rejection before polling in both public Live entry points |
| `src/smappy/plugins/fit.py` | File filter, help text, standalone default output naming |
| `src/smappy/gui/app.py` | Main window File > Open filter; drops use `is_image_stack()` |
| `src/smappy/batch.py` | Case-insensitive ND2 default discovery pattern |
| `pyproject.toml` | Optional `nd2` dependency group |

TIFF/NDTiff/DCIMG readers, camera/microscope databases, fitting engines, HDF5
writer and unrelated custom plugins are not replaced. When updating upstream,
check `ImageSource` fields and its `frame`/`frames` contracts, the above dispatch
and Live entry points, camera overlay semantics, and GUI `Plugin.describe()`.

## Validation in this checkout (2026-11-04)

Interpreter: `.venv/python.exe`, Python 3.12.16; `nd2==0.10.3`, NumPy 2.5.4.
The configured mirror's dependency resolver selected this compatible ND2
version. `pip check` passed and imports resolve to this target checkout.

The 36 ND2 tests passed, including the two opt-in real-file tests. Coverage
includes missing dependencies, all extension cases, singleton axes, rejected
dimensions/RGB, incomplete reads, handle release on errors/early close,
concurrent preview/iteration, metadata, manual overrides/presets, anisotropic
Gaussian/spline fits, Batch processing and HDF5 reload. Spline coverage uses a
synthetic bead calibration and synthetic emitter frames, testing pipeline
operation rather than a physical calibration for the real acquisition.

Real file (read only):
`G:/ND2 Test/HeLa R1-MK-DKO-Tt SNAP-R3-FL-BG-AF647-001.nd2`.

| Check | Result |
| --- | --- |
| Size | 1,488,359,424 bytes |
| Data | 11,000 frames, 256 x 256, uint16 |
| Metadata | 14 significant bits; calibrated X=Y=0.16 um; channel 647; Apo TIRF 100x Oil DIC N2 |
| Complete stream | 87 blocks of at most 127 frames |
| Largest pixel block | 16,646,144 bytes, excluding reader/pipeline overhead |
| Complete stream vs direct ND2 | Identical SHA-256 (below) |
| Random access | 13 frames across the acquisition, pixel-exact equality |
| Partial blocks | Frames 5500:5523, chunks of 7, pixel-exact equality |
| GUI | Offscreen Qt file-dialog callback, form hints and preview passed |
| Gaussian fit | First 32 frames; 3,987 candidates, 3,681 fitted ROIs, 3,546 accepted localizations; HDF5 written |
| Input integrity | File size and modification time unchanged |

Complete pixel-stream SHA-256:
`3bb4b80b57df472d287e92843a6f660ada896834d7eb013aa8903029d23dd7a8`.

The real Gaussian smoke test assumes conversion=1 e-/ADU, offset=100 ADU,
EM off, gain=1. These are **test settings, not measured camera calibration**.
All frames were read; only 32 were fitted. There was no manual desktop
click-through. ND2 0.10.3 emits NumPy 2.5 deprecation warnings about setting
array shape; pixel comparisons and fitting still pass.

The existing regression selection finished with **246 passed, 4 failed,
9 skipped**. It includes TIFF, NDTiff (including Live), single-image TIFF sets,
DCIMG, metadata, cameras/microscopes, Batch, fit plugins, previews, raw frames,
file drops and the public API. The four failures all came from the sandbox
denying `os.link()` when saving PSF calibration files (`WinError 5`). Each was
reproduced with the original HEAD versions of all four modified processing
modules loaded in memory, without changing the working tree:

- `test_microscope_db.py::test_a_saved_calibration_keeps_its_beads_microscope`
- `test_camera_frames.py::test_a_calibration_from_the_beads_fits_the_astigmatic_frames_in_z`
- `test_camera_frames.py::test_frames_drawn_with_a_measured_psf_fit_back_with_it`
- `test_drop_files.py::test_fit_this_gives_a_split_camera_drawn_with_a_dual_calibration_to_spline_3d_2c`

Those original tests and the calibration writers were left unchanged. The new
ND2 spline test saves its disposable calibration with the existing
`overwrite=True` option, avoiding hard links. Regression reports are retained
locally in `.venv/nd2-regression.xml` and `.venv/nd2-baseline.xml`. Skipped tests
belong to the upstream slow tier or require external reference datasets; the
entire repository suite was not run.

After adding the main-window file filter, the final ND2 + file-drop selection
reported **54 passing tests**, with the already confirmed hard-link failure
deselected (`.venv/nd2-final.xml`). Checking process status separately revealed
a Qt native crash *after* reporting results (Windows exit `-1073741819`,
`0xC0000005`). Therefore that combined invocation is not a clean process-level
pass. Isolated verification established:

- ND2: **36 passed, exit 0**, including both real-file tests
  (`.venv/nd2-isolated.xml`, `.venv/nd2-isolated.log`).
- Existing file-drop tests: **18 passed, 1 deselected**, followed by the same
  exit crash, even without running any ND2 tests (`.venv/drop-isolated.log`).
- Original HEAD versions of all five modified Python integration modules,
  including `gui/app.py`, loaded in memory: the same **18 passed, 1 deselected**
  and the same exit crash (`.venv/drop-baseline.log`). This is a pre-existing
  GUI test teardown issue in this environment; no GUI teardown code was changed.

Repeat the ND2 tests in PowerShell (all outputs go to the test directory):

```powershell
$env:SMAPPY_TEST_ND2 = 'G:\ND2 Test\HeLa R1-MK-DKO-Tt SNAP-R3-FL-BG-AF647-001.nd2'
.\.venv\python.exe -m pytest tests/test_nd2.py -q -s --basetemp=.venv/pytest-nd2
```

Use a disposable `--basetemp`: pytest clears it on reuse. Without
`SMAPPY_TEST_ND2`, the two real-file tests are skipped.

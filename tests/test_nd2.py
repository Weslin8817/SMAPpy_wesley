"""ND2 adapter contracts, plus opt-in verification against a real acquisition.

Set SMAPPY_TEST_ND2 to a completed single-channel ND2 to run the integration
test. It only reads that file; all output goes to pytest's temporary folder.
"""
import os
from contextlib import closing
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
import pytest

from smappy.io import nd2 as adapter
from smappy.io.tiff import open_stack, camera_metadata, is_image_stack, resolve_camera


@pytest.fixture
def fake(monkeypatch):
    class Reader:
        sizes = {"T": 7, "C": 1, "Y": 4, "X": 5}
        is_rgb = False
        dtype = np.dtype("uint16")
        attributes = NS(heightPx=4, widthPx=5, sequenceCount=7,
                        bitsPerComponentSignificant=14)
        calibrated = (True, True, False)
        xy = (0.16, 0.16)
        reads = []
        active = 0
        data = (np.arange(140).reshape(7, 4, 5) * 451).astype("uint16")

        def __init__(self, path):
            pass

        def __enter__(self):
            Reader.active += 1
            return self

        def __exit__(self, *args):
            Reader.active -= 1

        @property
        def metadata(self):
            return NS(channels=[NS(volume=NS(axesCalibrated=self.calibrated),
                                   channel=NS(name="647"),
                                   microscope=NS(objectiveName="100x"))])

        def voxel_size(self):
            return NS(x=self.xy[0], y=self.xy[1])

        def read_frame(self, i):
            self.reads.append(i)
            return self.data[i]

    monkeypatch.setattr(adapter, "_reader", Reader)
    return Reader


def test_lazy_dispatch_random_chunks_and_handle_lifetime(fake, tmp_path):
    # An unrelated TIFF in the same directory must not capture the ND2 path.
    import tifffile
    for i in range(2):
        tifffile.imwrite(tmp_path / f"img_{i}.tif", np.zeros((4, 5), "uint16"))
    source = open_stack(tmp_path / "test.ND2")
    from smappy.plugins.fit import default_output_path
    assert default_output_path(tmp_path / "test.ND2") == tmp_path / "test_locs.hdf5"
    assert isinstance(source, adapter.ND2Source)
    assert fake.reads == [] and fake.active == 0
    assert (len(source), source.shape, source.dtype) == (7, (4, 5), np.dtype("uint16"))
    np.testing.assert_array_equal(source.frame(6), fake.data[6])
    assert fake.reads == [6] and fake.active == 0
    it = source.frames(chunk=2, start=1, stop=6)
    first, block = next(it)
    assert first == 1 and fake.reads == [6, 1, 2]
    np.testing.assert_array_equal(block, fake.data[1:3])
    it.close()
    assert fake.active == 0
    blocks = list(source.frames(chunk=3, start=2, stop=99))
    assert [i for i, _ in blocks] == [2, 5]
    got = np.concatenate([b for _, b in blocks])
    assert got.dtype == fake.dtype
    np.testing.assert_array_equal(got, fake.data[2:])
    assert list(source.frames(start=7)) == []
    assert list(source.frames(start=3, stop=2)) == []
    for i in (-1, 7):
        with pytest.raises(IndexError):
            source.frame(i)
    for args in ({"chunk": 0}, {"start": -1}, {"stop": -1}):
        with pytest.raises(ValueError):
            list(source.frames(**args))


@pytest.mark.parametrize("axis", ["C", "Z", "P", "S"])
def test_unsupported_dimensions_are_not_flattened(fake, axis):
    fake.sizes = {**fake.sizes, axis: 2}
    with pytest.raises(ValueError, match="single-channel"):
        open_stack("test.nd2")
    assert fake.active == 0


def test_metadata_and_manual_camera_values(fake):
    source = open_stack("test.nd2")
    cam = camera_metadata(source, require=False)
    assert cam.pixelsize_um == 0.16
    assert cam.conversion is None and cam.offset is None
    assert cam.em_on is None and cam.emgain is None
    with pytest.raises(ValueError):
        camera_metadata(source)
    cam = camera_metadata(source, overrides={"conversion": 2, "offset": 100,
                          "em_on": True, "emgain": 50, "pixelsize_um": 0.2})
    assert (cam.conversion, cam.offset, cam.emgain, cam.pixelsize_um) == (2, 100, 50, 0.2)
    fake.xy = (0.16, 0.18)
    assert camera_metadata(open_stack("test.nd2"), require=False).pixelsize_um_xy == (0.16, 0.18)
    fake.calibrated = (False, False, False)
    assert camera_metadata(open_stack("test.nd2"), require=False).pixelsize_um is None


def test_missing_library_has_actionable_error(monkeypatch):
    import sys
    monkeypatch.setitem(sys.modules, "nd2", None)
    with pytest.raises(ImportError, match="pip install"):
        adapter._reader("test.nd2")


def test_missing_library_does_not_affect_recognition_or_tiff(monkeypatch, tmp_path):
    import sys
    import tifffile
    monkeypatch.setitem(sys.modules, "nd2", None)
    assert is_image_stack(tmp_path / "test.Nd2")
    path = tmp_path / "plain.tif"
    data = np.arange(20, dtype="uint16").reshape(4, 5)
    tifffile.imwrite(path, data)
    np.testing.assert_array_equal(open_stack(path).frame(0), data)
    with pytest.raises(ImportError, match="pip install"):
        open_stack(tmp_path / "test.nd2")


@pytest.mark.parametrize("suffix", ["nd2", "ND2", "nD2", "Nd2"])
def test_recognition_and_batch_case_variants(fake, tmp_path, suffix):
    from smappy.batch import Input, Job, Output, resolve_inputs
    from smappy.chain import ChainSpec
    from smappy.plugins.fit import default_output_path
    from smappy.gui.drops import kind_of, STACK
    path = tmp_path / f"test.{suffix}"
    path.touch()
    assert is_image_stack(path)
    assert kind_of(path) == STACK
    assert isinstance(open_stack(path), adapter.ND2Source)
    assert default_output_path(path) == tmp_path / "test_locs.hdf5"
    job = Job(chain=ChainSpec(), inputs=[Input(folder=str(tmp_path))],
              output=Output(folder=str(tmp_path / "out")))
    assert [item.path for item in resolve_inputs(job, "images")] == [path]
    job.inputs = [Input(file=str(path))]
    assert [item.path for item in resolve_inputs(job, "images")] == [path]


def test_nd2_next_to_ndtiff_is_not_a_directory_acquisition(fake, tmp_path):
    from smappy.plugins.fit import default_output_path
    (tmp_path / "NDTiff.index").touch()
    path = tmp_path / "own.nd2"
    assert isinstance(open_stack(path), adapter.ND2Source)
    assert default_output_path(path) == tmp_path / "own_locs.hdf5"


@pytest.mark.parametrize("calibrated", [(), (True,), (True, False), (False, True)])
def test_both_axes_must_be_explicitly_calibrated(fake, calibrated):
    fake.calibrated = calibrated
    assert camera_metadata(open_stack("test.nd2"), require=False).pixelsize_um is None


@pytest.mark.parametrize("xy", [(0, .16), (-1, .16), (.16, float("inf"))])
def test_invalid_pixel_sizes_stay_unset(fake, xy):
    fake.xy = xy
    assert camera_metadata(open_stack("test.nd2"), require=False).pixelsize_um is None


def test_missing_descriptive_metadata_does_not_hide_calibration(fake, monkeypatch):
    monkeypatch.setattr(fake, "metadata", property(lambda self: NS(channels=[
        NS(volume=NS(axesCalibrated=(True, True)))])))
    assert camera_metadata(open_stack("test.nd2"), require=False).pixelsize_um == .16
    monkeypatch.setattr(fake, "metadata", None)
    with pytest.warns(RuntimeWarning, match="optional metadata"):
        source = open_stack("test.nd2")
    assert camera_metadata(source, require=False).pixelsize_um is None
    assert fake.active == 0


def test_explicit_camera_preset_and_user_overrides_still_work(fake, monkeypatch):
    from smappy import camera_db
    from smappy.camera_db import Camera, CameraDatabase, Parameter
    db = CameraDatabase([Camera(name="ND2 test camera", parameters={
        key: Parameter(fixed=value) for key, value in dict(
            conversion=2., offset=100., em_on=True, emgain=50., pixelsize_um=.1).items()})])
    monkeypatch.setattr(camera_db, "database", lambda: db)
    source = open_stack("test.nd2")
    resolved = resolve_camera(source, camera="ND2 test camera")
    cam = resolved.camera_metadata()
    assert (cam.conversion, cam.offset, cam.em_on, cam.emgain) == (2., 100., True, 50.)
    assert cam.pixelsize_um == .16
    assert resolved.sources["pixelsize_um"].kind == "metadata"
    assert camera_metadata(source, camera="ND2 test camera",
                           overrides={"pixelsize_um": (.2, .3)}).pixelsize_um_xy == (.2, .3)


def test_rgb_and_empty_sequences_are_rejected(fake):
    fake.is_rgb = True
    with pytest.raises(ValueError, match="single-channel"):
        open_stack("test.nd2")
    fake.is_rgb = False
    fake.attributes.sequenceCount = 0
    with pytest.raises(ValueError, match="empty or inconsistent"):
        open_stack("test.nd2")
    assert fake.active == 0


@pytest.mark.parametrize("bad", [np.zeros((3, 5), "uint16"), np.zeros((4, 5), "uint8")])
def test_bad_frame_releases_reader(fake, monkeypatch, bad):
    source = open_stack("test.nd2")
    monkeypatch.setattr(fake, "read_frame", lambda self, i: bad)
    with pytest.raises(ValueError, match="inconsistent"):
        source.frame(0)
    with pytest.raises(ValueError, match="inconsistent"):
        list(source.frames(chunk=2))
    assert fake.active == 0


def test_read_failure_releases_reader(fake, monkeypatch):
    source = open_stack("test.nd2")
    def broken(self, i):
        raise OSError("incomplete ND2 frame")
    monkeypatch.setattr(fake, "read_frame", broken)
    with pytest.raises(OSError, match="incomplete"):
        list(source.frames())
    assert fake.active == 0


def test_preview_and_background_iteration_own_separate_readers(fake, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from threading import Event, get_ident
    ready, release = Event(), Event()
    handles = []
    read = fake.read_frame
    def tracked(self, i):
        owner = getattr(self, "owner", get_ident())
        assert owner == get_ident()
        self.owner = owner
        if not any(self is handle for handle in handles):
            handles.append(self)
        return read(self, i)
    monkeypatch.setattr(fake, "read_frame", tracked)
    source = open_stack("test.nd2")
    def background():
        with closing(source.frames(chunk=2)) as frames:
            first, block = next(frames)
            ready.set()
            assert release.wait(10)
            return block
    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(background)
        try:
            assert ready.wait(10)
            np.testing.assert_array_equal(source.frame(6), fake.data[6])
            assert fake.active == 1
        finally:
            release.set()
        np.testing.assert_array_equal(future.result(), fake.data[:2])
    assert len(handles) == 2 and fake.active == 0


def test_all_live_entry_points_reject_nd2_without_waiting(fake):
    from smappy.io.watch import open_growing_stack, watch_stack
    from smappy.plugins.fit import SourceSettings, _open
    for path in ("test.nd2", "not-created.ND2"):
        for call in (lambda: open_growing_stack(path), lambda: next(watch_stack(path)),
                     lambda: _open(SourceSettings(path=path, live=True), watch=True)):
            with pytest.raises(ValueError, match="Live ND2"):
                call()
    with pytest.raises(ValueError, match="Live ND2"):
        open_stack("test.nd2").watch()


def test_single_frame_singleton_spatial_axes(fake):
    fake.sizes = {"Y": 1, "X": 5}
    fake.attributes = NS(heightPx=1, widthPx=5, sequenceCount=1,
                         bitsPerComponentSignificant=16)
    fake.data = np.array([[0, 1, 32768, 65534, 65535]], dtype="uint16")
    source = open_stack("single.nd2")
    assert len(source) == 1 and source.shape == (1, 5)
    np.testing.assert_array_equal(source.frame(0), fake.data)
    assert next(source.frames())[1].shape == (1, 1, 5)


def test_bad_pixel_calibration_and_inconsistent_length(fake):
    fake.xy = (float("nan"), 0.16)
    assert camera_metadata(open_stack("test.nd2"), require=False).pixelsize_um is None
    fake.sizes = {**fake.sizes, "T": 8}
    with pytest.raises(ValueError, match="inconsistent"):
        open_stack("test.nd2")
    assert fake.active == 0


def test_gui_filter_batch_discovery_and_live_guard(fake, tmp_path):
    from smappy.plugins.fit import SourceSettings, TIFF_FILTER, _open
    from smappy.batch import Input, Job, Output, resolve_inputs
    from smappy.chain import ChainSpec
    from smappy.gui.app import open_filter
    assert "*.nd2" in TIFF_FILTER
    assert all(f"*.{ext}" in open_filter() and f"*.{ext}" in TIFF_FILTER
               for ext in ("nd2", "ND2", "nD2", "Nd2"))
    path = tmp_path / "test.nd2"
    path.touch()
    job = Job(chain=ChainSpec(), inputs=[Input(folder=str(tmp_path))],
              output=Output(folder=str(tmp_path / "out")))
    assert [item.path for item in resolve_inputs(job, "images")] == [path]
    with pytest.raises(ValueError, match="Live ND2"):
        _open(SourceSettings(path=str(path), live=True), watch=True)


@pytest.fixture
def emitter_nd2(fake, tmp_path):
    """Small deterministic photon-like frames, exercising the real fitting engine."""
    fake.sizes = {"T": 6, "Y": 48, "X": 48}
    fake.attributes = NS(heightPx=48, widthPx=48, sequenceCount=6,
                         bitsPerComponentSignificant=16)
    y, x = np.mgrid[:48, :48]
    clean = 100 + 3000 * np.exp(-((x - 22.3)**2 + (y - 24.2)**2) / (2 * 1.3**2))
    fake.data = np.random.default_rng(4).poisson(clean, (6, 48, 48)).astype("uint16")
    fake.xy = (.16, .18)
    path = tmp_path / "emitters.ND2"
    path.touch()
    return path


@pytest.mark.parametrize("spline", [False, True], ids=["gaussian", "spline"])
def test_fitting_preview_and_hdf5_with_anisotropic_pixels(emitter_nd2, fake, tmp_path, spline):
    from smappy.io.hdf5 import load_localizations
    from smappy.plugins import Context
    from smappy.plugins.fit import (CameraSettings, GaussianFit, GaussianFitSettings,
                                    OutputSettings, SourceSettings, SplineFit,
                                    SplineFitSettings, SplineModelSettings)
    common = dict(source=SourceSettings(path=str(emitter_nd2), chunk=2, start=1, stop=5),
                  camera=CameraSettings(conversion=1., offset=100., em_on=False, emgain=1.),
                  output=OutputSettings(path=str(tmp_path / "fitted.hdf5"), raw_frames=2))
    if spline:
        from smappy.calibrate.core import CalibrationSettings, build_calibration, collect_beads
        from smappy.calibrate.input import BeadStack
        from smappy.io.calibration import save_spline_calibration
        from smappy.simulate import bead_stacks
        stacks, z = bead_stacks(1, seed=0)
        beads = collect_beads([BeadStack(s.astype(np.float32), z) for s in stacks],
                              CalibrationSettings(roi_size=21))
        calibration = build_calibration(beads).calibration
        cal_path = tmp_path / "synthetic-calibration.h5"
        save_spline_calibration(cal_path, calibration, overwrite=True)
        settings = SplineFitSettings(**common, model=SplineModelSettings(
            calibration=str(cal_path), focus="calibrated"))
        plugin = SplineFit()
    else:
        settings, plugin = GaussianFitSettings(**common), GaussianFit()
    hints = plugin.describe(None, settings).hints
    assert hints["camera.pixelsize_um"] == .16
    assert hints["camera.pixelsize_y_um"] == .18
    preview = plugin.preview(Context(), settings, frame=2)
    assert preview.plot is not None and preview.data["error"] == ""
    result = plugin.run(Context(), settings)
    assert len(result.locs) == 4
    np.testing.assert_array_equal(result.locs["frame"], [1, 2, 3, 4])
    assert np.isfinite(result.locs["x_nm"]).all()
    if spline:
        assert np.isfinite(result.locs["z_nm"]).all()
    else:
        np.testing.assert_allclose(result.locs["x_nm"], 22.3 * 160, atol=30)
        np.testing.assert_allclose(result.locs["y_nm"], 24.2 * 180, atol=30)
    saved = load_localizations(tmp_path / "fitted.hdf5")
    assert len(saved) == len(result.locs)
    np.testing.assert_array_equal(saved["x_nm"], result.locs["x_nm"])
    assert saved.metadata["acquisition"]["summary"]["Format"] == "Nikon ND2"
    assert saved.metadata["acquisition"]["summary"]["ND2-PixelSizeUmXY"] == [.16, .18]
    assert fake.active == 0


def test_batch_runs_nd2_through_fitting_and_output(emitter_nd2, tmp_path, fake):
    from smappy import batch
    from smappy.chain import ChainSpec, Step
    from smappy.io.hdf5 import load_localizations
    job = batch.Job(chain=ChainSpec(steps=[Step(plugin="Localize/Gaussian 2D", label="fit",
        values={"camera.conversion": 1., "camera.offset": 100.,
                "camera.em_on": False, "camera.emgain": 1., "source.chunk": 2})]),
        inputs=[batch.Input(folder=str(tmp_path))],
        output=batch.Output(folder=str(tmp_path / "batch-output"), figures="none"))
    assert batch.validate(job) == []
    report = batch.run(job)
    assert [f["status"] for f in report["files"]] == ["done"], report
    fitted = load_localizations(tmp_path / "emitters_locs.hdf5")
    processed = load_localizations(tmp_path / "batch-output" / "emitters" / "emitters_proc.h5")
    assert len(fitted) == len(processed) == 6
    assert fake.active == 0


@pytest.mark.skipif(not os.environ.get("SMAPPY_TEST_ND2"), reason="set SMAPPY_TEST_ND2")
def test_real_complete_pixel_stream(tmp_path):
    import hashlib
    import nd2
    path = Path(os.environ["SMAPPY_TEST_ND2"])
    before = path.stat()
    source = open_stack(path)
    streamed, direct = hashlib.sha256(), hashlib.sha256()
    count = blocks = max_bytes = 0
    for first, block in source.frames(chunk=127):
        assert first == count and len(block) <= 127
        streamed.update(block.tobytes())
        count += len(block)
        blocks += 1
        max_bytes = max(max_bytes, block.nbytes)
    with nd2.ND2File(path) as reference:
        for i in range(len(source)):
            direct.update(reference.read_frame(i).tobytes())
    assert count == len(source)
    assert streamed.digest() == direct.digest()
    after = path.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    print(f"FULL ND2: {count} frames, {blocks} blocks, max block={max_bytes} bytes; "
          f"SHA256={streamed.hexdigest()}")


@pytest.mark.skipif(not os.environ.get("SMAPPY_TEST_ND2"), reason="set SMAPPY_TEST_ND2")
def test_real_nd2_reading_gui_and_fit(tmp_path, monkeypatch):
    import nd2
    from smappy.plugins import Context
    from smappy.plugins.fit import (CameraSettings, GaussianFit, GaussianFitSettings,
                                    OutputSettings, SourceSettings)
    path = Path(os.environ["SMAPPY_TEST_ND2"])
    before = path.stat()
    source = open_stack(path)
    with nd2.ND2File(path) as reference:
        assert len(source) == reference.sizes.get("T", 1)
        assert source.dtype == reference.dtype
        assert source.shape == (reference.attributes.heightPx, reference.attributes.widthPx)
        indices = sorted(set([0, len(source) - 1, len(source) // 2,
                              *np.random.default_rng(42).integers(0, len(source), 10)]))
        for i in indices:
            np.testing.assert_array_equal(source.frame(int(i)), reference.read_frame(int(i)))
        start, stop = len(source) // 2, min(len(source), len(source) // 2 + 23)
        for first, block in source.frames(chunk=7, start=start, stop=stop):
            assert len(block) <= 7 and block.dtype == source.dtype
            for j, frame in enumerate(block):
                np.testing.assert_array_equal(frame, reference.read_frame(first + j))
        assert camera_metadata(source, require=False).pixelsize_um_xy == (
            reference.voxel_size().x, reference.voxel_size().y)

    # Exercise the GUI form and the same plugin the Run button calls. Camera
    # values below are smoke-test assumptions, NOT measured camera calibration.
    from PySide6.QtWidgets import QApplication, QFileDialog
    from smappy.gui.plugin_panel import PluginPanel
    from smappy.session import Session
    app = QApplication.instance() or QApplication([])
    panel = PluginPanel(GaussianFit, Session())
    def choose_nd2(parent, title, current, name_filter):
        assert "*.nd2" in name_filter
        return str(path), name_filter
    monkeypatch.setattr(QFileDialog, "getOpenFileName", choose_nd2)
    panel.form.fields["source"].fields["path"]._browse()
    assert panel.form.value().source.path == str(path)
    settings = GaussianFitSettings(
        source=SourceSettings(path=str(path), stop=min(32, len(source)), chunk=7),
        camera=CameraSettings(conversion=1.0, offset=100.0, em_on=False, emgain=1.0),
        output=OutputSettings(path=str(tmp_path / "nd2_smoke.hdf5"), raw_frames=3))
    panel.form.set(settings)
    assert panel.form.value().source.path == str(path)
    plugin = GaussianFit()
    assert plugin.describe(None, settings).hints["camera.pixelsize_um"] == camera_metadata(
        source, require=False).pixelsize_um_xy[0]
    preview = plugin.preview(Context(), settings, frame=len(source) // 2)
    assert preview.plot is not None and preview.data["error"] == ""
    assert len(preview.data["locs"]) > 0
    result = plugin.run(Context(), settings)
    assert len(result.locs) > 0
    assert np.isfinite(result.locs["x_nm"]).all()
    assert np.isfinite(result.locs["y_nm"]).all()
    assert result.locs["frame"].min() >= 0
    assert result.locs["frame"].max() < 32
    assert (tmp_path / "nd2_smoke.hdf5").exists()
    print(f"REAL ND2: {len(source)} frames, {source.shape}, {source.dtype}; "
          f"random frames={indices}; smoke localizations={len(result.locs)}; "
          f"stats={result.data['stats']}")
    panel.close()
    panel.deleteLater()
    app.processEvents()
    after = path.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)

"""The migration boundary: discovery, unchanged algorithms, frozen input."""
import ast
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import threading
import time

import numpy as np
import pytest

from smappy import config, docs, plugins
from smappy.locs import Localizations
from smappy.plugins import Context, Plugin
from smappy.regions import Region
from smappy.session import Session
from smappy.view3d import Slab
from smappy_centerline import plugin as implementation
from smappy_centerline.plugin import CenterlineRadial, Settings

ROOT = Path(__file__).resolve().parents[2]
PATH = "Analysis/Measure/Radial distance and density"


def session():
    x = np.repeat(np.linspace(-500, 500, 80), 4)
    y = np.tile([-15., -5., 5., 15.], 80)
    locs = Localizations(dict(x_nm=x, y_nm=y, z_nm=y * 2,
                             frame=np.arange(len(x)) * 10,
                             photons=np.full(len(x), 1000.),
                             xy_err_nm=np.full(len(x), 3.),
                             z_err_nm=np.full(len(x), 8.)))
    s = Session(locs)
    s.roi = Region.rect(-450, -100, 450, 100)
    return s


def test_source_algorithms_defaults_and_analysis_flow_unchanged():
    manifest = json.loads(Path(__file__).with_name("source_manifest.json").read_text())
    folder = Path(implementation.__file__).parent
    for name, digest in manifest["normalized_sha256"].items():
        assert hashlib.sha256((folder / name).read_text(encoding="utf-8").encode()).hexdigest() == digest
    tree = ast.parse(Path(implementation.__file__).read_text(encoding="utf-8"))
    functions = {n.name: hashlib.sha256(ast.dump(n).encode()).hexdigest()
                 for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    for name, digest in manifest["plugin_function_ast_sha256"].items():
        assert functions[name] == digest
    assert asdict(Settings()) == manifest["settings_defaults"]


def test_official_discovery_load_help_and_execution(monkeypatch):
    with monkeypatch.context() as patch:
        patch.setattr(config, "plugin_roots", lambda: [ROOT / "user_plugins"])
        found = plugins.discover(force=True)
        ref = found[PATH]
        assert ref.origin == ROOT / "user_plugins/Analysis/Measure/centerline_radial.py"
        cls = ref.load()
        assert issubclass(cls, CenterlineRadial) and not cls.changes_locs
        assert docs.page_for(cls).path == ref.origin.with_suffix(".md")
        s = session()
        result = s.run(cls(), Settings(density=False))
        assert result.data["N"] > 0 and result.locs is None
        assert result.data["R100_nm"] > result.data["R50_nm"]
    plugins.discover(force=True)


def test_normal_startup_discovers_checkout_plugin_without_configuration(monkeypatch, tmp_path):
    with monkeypatch.context() as patch:
        patch.setattr(config, "plugin_roots", lambda: [])
        patch.setattr(config, "user_plugin_dir", lambda: tmp_path / "empty-user-plugins")
        patch.chdir(tmp_path)
        ref = plugins.discover(force=True)[PATH]
        assert ref.root == "checkout"
        assert ref.origin == ROOT / "user_plugins/Analysis/Measure/centerline_radial.py"
        assert issubclass(ref.load(), CenterlineRadial)
        assert docs.page_for(ref.load()).path == ref.origin.with_suffix(".md")
        assert not (tmp_path / "empty-user-plugins").exists()
    plugins.discover(force=True)


@pytest.mark.parametrize("cached", [False, True])
@pytest.mark.parametrize("slab_on", [False, True])
def test_snapshot_preserves_input_and_uses_grouped_filter_roi_slab(cached, slab_on):
    s = session()
    if cached:
        s.table(grouping="grouped")
        # A genuine stale cache: change the table's linking key.
        s.set_group_settings(replace(s.group_settings, dt=2))
    s.add_layer(like=0)
    s.render_layer_index = 1
    s.layers[1].set_bound("x_nm", -300., 300.)
    s.slab = Slab.from_bounds(-400, 400, -100, 100, -12, 12)
    s.select_in_slab = slab_on
    s.slab_shown = True
    before = {k: v.copy() for k, v in s.locs.columns.items()}
    state = s.layer_configs()
    roi, slab = repr(s.roi), repr(s.slab)
    cached_table = s.tables.peek_grouped()
    expected_session = s.scratch()
    table, selected = expected_session.table(1, grouping="grouped")
    plugin = CenterlineRadial()
    snapshot = plugin.prepare_context(s.context(grouping=plugin.grouping), Settings())
    assert snapshot.session is None
    assert snapshot.selection.layer == 1
    np.testing.assert_array_equal(snapshot.centerline_indices, selected.indices)
    for key, value in snapshot.locs.columns.items():
        np.testing.assert_array_equal(value, table[key][selected.mask])
    if slab_on:
        assert np.all(np.abs(snapshot.locs["z_nm"]) <= 12)
    assert snapshot.centerline_provenance["slab"] == (slab if slab_on else "")
    assert s.layer_configs() == state
    assert (repr(s.roi), repr(s.slab)) == (roi, slab)
    assert s.tables.peek_grouped() is cached_table
    assert set(before) == set(s.locs.columns)
    for key, value in before.items():
        np.testing.assert_array_equal(s.locs[key], value)
    frozen = snapshot.locs["x_nm"].copy()
    s.locs["x_nm"][:] = 9999
    s.roi = None
    np.testing.assert_array_equal(snapshot.locs["x_nm"], frozen)
    assert plugin.run(snapshot, Settings(density=False)).data["N"] == len(frozen)


def test_image_layer_refused():
    from smappy.images import ImageData
    from smappy.session import Layer
    s = session()
    s.layers.append(Layer.from_image(ImageData(np.ones((4, 4)), pixelsize=10)))
    s.render_layer_index = len(s.layers) - 1
    with pytest.raises(ValueError, match="localization layer"):
        CenterlineRadial().prepare_context(s.context(grouping="ungrouped"), Settings())


def test_default_plugin_hooks_are_noops():
    class Measurement(Plugin):
        changes_locs = False
    ctx = Context()
    assert Measurement().prepare_context(ctx, None) is ctx
    assert Measurement().create_result_widget(None) is None


@pytest.mark.parametrize("manual_cutoff", [1., 10000.])
def test_auto_xyz_rejection_is_shared_with_independent_fit_settings(monkeypatch, manual_cutoff):
    from smappy.plugins import Selection
    x = np.repeat(np.linspace(-500, 500, 40), 24)
    angles = np.tile(np.repeat(np.arange(8) * np.pi / 4, 3), 40)
    radii = np.tile([10., 20., 30.], 320)
    points = np.column_stack([x, radii * np.cos(angles), radii * np.sin(angles)])
    points = np.vstack([points, [-200., 0., 300.], [200., 0., -300.]])
    ctx = Context(locs=Localizations({f"{a}_nm": points[:, i] for i, a in enumerate("xyz")}),
                  selection=Selection(np.ones(len(points), bool), roi=Region.rect(-600, -100, 600, 100)))
    calls = []
    fit = implementation.fit_selected_centerline
    def tracked(p, precision, **kwargs):
        calls.append((p.copy(), kwargs))
        return fit(p, precision, **kwargs)
    monkeypatch.setattr(implementation, "fit_selected_centerline", tracked)
    plugin = CenterlineRadial()
    result = plugin.run(ctx, Settings(rejection="Auto conservative", maximum_distance_nm=manual_cutoff,
                                     dimension="2D", bins=8, smoothing=2.,
                                     density_dimension="3D", density_bins=12, density_smoothing=.5))
    entry = plugin.results.entries["ROI_001"]
    assert result.data["N_noise_excluded"] == result.data["density_noise_excluded"] == 2
    assert result.data["N_after_rejection"] == 960
    assert len(calls) == 3
    assert sum(k["rejection"] != "None" for _, k in calls) == 1
    assert calls[0][0].shape[1] == 3
    assert (calls[1][1]["bins"], calls[1][1]["smoothing"]) == (8, 2.)
    assert (calls[2][1]["bins"], calls[2][1]["smoothing"]) == (12, .5)
    np.testing.assert_array_equal(calls[1][0], points[:-2, :2])
    np.testing.assert_array_equal(calls[2][0], points[:-2])
    np.testing.assert_array_equal(entry.indices, entry.density_entry.indices)


def test_gui_snapshot_precedes_worker_and_preserves_data(monkeypatch):
    from PySide6.QtWidgets import QApplication
    from smappy.gui.plugin_panel import PluginPanel
    app = QApplication.instance() or QApplication([])
    s = session()
    original = {k: v.copy() for k, v in s.locs.columns.items()}
    panel = PluginPanel(CenterlineRadial, s)
    panel.form.set(Settings(density=False))
    monkeypatch.setattr(panel, "plot", lambda **kw: None)
    waiting, proceed = threading.Event(), threading.Event()
    run = panel.plugin.run
    owner = threading.get_ident()
    threads = []
    def delayed(ctx, settings):
        threads.append(threading.get_ident())
        assert ctx.session is None
        waiting.set()
        assert proceed.wait(10)
        return run(ctx, settings)
    monkeypatch.setattr(panel.plugin, "run", delayed)
    try:
        panel.run()
        assert waiting.wait(10)
        s.roi = None
        s.layers[0].set_bound("photons", 0, 1)
        proceed.set()
        deadline = time.monotonic() + 20
        while (panel.busy or not panel.run_button.isEnabled()) and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(.01)
        app.processEvents()
        assert not panel.busy
        assert threads == [threads[0]] and threads[0] != owner
        assert panel.result_widget.table.rowCount() == 1
        entry = panel.plugin.results.entries["ROI_001"]
        assert entry.summary["N"] > 200
        assert set(original) == set(s.locs.columns)
        for key, value in original.items():
            np.testing.assert_array_equal(value, s.locs[key])
    finally:
        proceed.set()
        if panel._thread is not None:
            panel._thread.quit()
            panel._thread.wait(20000)
        panel.release()
        panel.close()

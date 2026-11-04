"""Real plugin panel, worker snapshot and result navigation."""
import time

import numpy as np
import pytest

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication, QFileDialog

from smappy.gui.plugin_panel import PluginPanel
from smappy.gui.render_tab import RenderTab
from smappy.locs import Localizations
from smappy_centerline.plugin import CenterlineRadial
from smappy.regions import Region
from smappy.session import Session


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def wait_for_worker(panel, app):
    deadline = time.monotonic()+20
    while time.monotonic() < deadline:
        app.processEvents()
        if panel._thread is not None and not panel._thread.isRunning() and panel.run_button.isEnabled():
            return
        time.sleep(.01)
    pytest.fail("centerline worker did not finish")


def test_panel_accumulates_selects_deletes_clears_and_exports(app, tmp_path, monkeypatch):
    locs = Localizations({"x_nm": np.repeat(np.linspace(-500, 500, 100), 4),
                          "y_nm": np.tile([-15., -5., 5., 15.], 100),
                          "frame": np.arange(400)*10})
    session = Session(locs)
    session.set_roi(Region.rect(-600, -100, 600, 100))
    panel = PluginPanel(CenterlineRadial, session)
    widget = panel.result_widget
    try:
        for n in (1, 2):
            panel.run()
            wait_for_worker(panel, app)
            assert widget.table.rowCount() == n, panel.output.toPlainText()
        widget.table.selectRow(0)
        assert panel._window is not None
        assert len(panel.result.figures()) == 2
        monkeypatch.setattr(QFileDialog, "getExistingDirectory", lambda *a: str(tmp_path))
        widget.export()
        assert (tmp_path/"centerline_radial_all_distances.csv").exists()
        widget.delete_selected()
        assert widget.table.rowCount() == 1
        assert list(panel.plugin.results.entries) == ["ROI_002"]
        widget.table.selectRow(0)
        widget.clear()
        assert widget.table.rowCount() == 0
        assert not panel.plugin.results.entries
        assert not panel.result_button.isEnabled()
        session.set_roi(None)
        panel.run()
        assert "ROI" in panel.status.text()
        assert widget.table.rowCount() == 0
    finally:
        if panel._window is not None:
            panel._window.close()
        panel.release()
        panel.close()


def test_render_tracks_selected_layer(app):
    session = Session(Localizations({"x_nm": np.arange(50.), "y_nm": np.arange(50.)}))
    session.add_layer(like=0)
    render = RenderTab(session)
    render._bind_layer(1)
    assert session.render_layer_index == 1
    render._bind_layer(0)
    assert session.render_layer_index == 0
    render.close()


def test_three_tables_keep_roi_identity_and_full_summary_export(app, tmp_path, monkeypatch):
    import csv
    from smappy_centerline.centerline import CenterlineResult
    panel = PluginPanel(CenterlineRadial, Session(Localizations({"x_nm": np.arange(50.), "y_nm": np.arange(50.)})))
    monkeypatch.setattr(panel, "plot", lambda: None)
    widget = panel.result_widget
    try:
        for radial, density, dimension in ((True, False, "3D"), (False, True, "2D"), (True, True, "3D")):
            entry = CenterlineResult(np.empty((0, 2)), np.empty(0, int), np.empty((0, 2)),
                np.empty(0), np.empty(0), np.empty((0, 2)),
                dict(radial_enabled=radial, density_enabled=density, dimension="2D",
                     N=50 if radial else 0, density_dimension=dimension,
                     density_status="exploratory_warning", density_warnings="Inspect sampling",
                     Area_um2=.2 if dimension == "2D" else None,
                     Volume_um3=.01 if dimension == "3D" else None,
                     density_units="groups/um^2" if dimension == "2D" else "groups/um^3",
                     rms_radial_distance_nm=12.5, noise_rejection="None",
                     N_before_rejection=50, N_after_rejection=50, N_noise_excluded=0))
            widget.store.add(entry)
        widget.refresh("ROI_002")
        for table in widget.tables:
            assert [table.item(row, 0).text() for row in range(3)] == ["ROI_001", "ROI_002", "ROI_003"]
            assert table.selectionModel().selectedRows()[0].row() == 1
        assert widget.radial_table.item(1, 1).text() == "Not enabled"
        assert widget.radial_table.item(1, 3).text() == "—"
        assert widget.density_table.item(0, 1).text() == "Not enabled"
        assert widget.density_table.item(1, 5).text() == "0.2"
        assert widget.density_table.item(1, 6).text() == "µm²"
        assert widget.density_table.item(2, 6).text() == "µm³"
        assert "rms_radial_distance_nm: 12.5" in widget.details.toPlainText()
        assert widget.density_table.item(2, 1).toolTip() == "Inspect sampling"
        monkeypatch.setattr(QFileDialog, "getExistingDirectory", lambda *a: str(tmp_path))
        widget.export()
        with (tmp_path / "centerline_radial_summary.csv").open(encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 3
        assert rows[1]["ROI"] == "ROI_002"
        assert rows[1]["radial_status"] == "disabled"
        assert rows[1]["density_area_um2"] == "0.2"
        assert rows[0]["radial_rms_nm"] == "12.5"
        widget.density_table.selectRow(0)
        widget.delete_selected()
        for table in widget.tables:
            assert [table.item(row, 0).text() for row in range(2)] == ["ROI_002", "ROI_003"]
        widget.radial_table.selectRow(1)
        assert all(t.selectionModel().selectedRows()[0].row() == 1 for t in widget.tables)
        widget.clear()
        assert all(t.rowCount() == 0 for t in widget.tables)
        assert not widget.details.toPlainText()
    finally:
        panel.release()
        panel.close()


def test_shared_noise_counts_are_distinct_from_analysis_validity():
    from smappy.plugins import Context, Selection
    from smappy_centerline.plugin import Settings
    x = np.linspace(-500, 500, 120)
    z = np.zeros(120)
    z[0] = np.nan
    z[60] = 300.
    ctx = Context(locs=Localizations(dict(x_nm=x, y_nm=np.sin(x/100)*5, z_nm=z)),
                  selection=Selection(np.ones(120, bool), roi=Region.rect(-600, -100, 600, 100)))
    plugin = CenterlineRadial()
    rejected = plugin.run(ctx, Settings(density=False, rejection="Manual maximum distance"))
    values = rejected.data
    assert values["noise_excluded_nonfinite"] == 1
    assert values["N_noise_excluded"] >= 1
    assert values["N_before_rejection"] == (values["N_after_rejection"] +
        values["N_noise_excluded"] + values["noise_excluded_nonfinite"])
    untouched = plugin.run(ctx, Settings(density=False))
    assert untouched.data["noise_excluded_nonfinite"] is None
    assert untouched.data["N_after_rejection"] == 120


def test_analysis_sections_and_independent_controls(app):
    from PySide6.QtWidgets import QLabel
    panel = PluginPanel(CenterlineRadial, Session(Localizations({"x_nm": np.arange(50.), "y_nm": np.arange(50.)})))
    try:
        fields = panel.form.fields
        assert panel.form.value().radial and panel.form.value().density
        headings = [w.text() for w in panel.form.findChildren(QLabel)]
        assert all(title in headings for title in
                   ("1. 3D noise rejection", "2. Radial Distance", "3. Labeling Density"))
        fields["density_dimension"].set("2D")
        fields["dimension"].set("3D")
        fields["radial"].set(False)
        assert not fields["dimension"].isEnabled()
        assert fields["density_dimension"].isEnabled()
        assert panel.run_button.isEnabled()
        fields["density"].set(False)
        assert not panel.run_button.isEnabled()
        fields["radial"].set(True)
        assert panel.run_button.isEnabled()
        assert panel.form.value().dimension == "3D"
        assert panel.form.value().density_dimension == "2D"
        fields["rejection"].set("Manual maximum distance")
        assert not fields["maximum_distance_nm"].isHidden()
        fields["rejection"].set("None")
        assert fields["maximum_distance_nm"].isHidden()
    finally:
        panel.release()
        panel.close()


def test_density_panel_worker_and_segment_plot(app, tmp_path, monkeypatch):
    from test_tube_density import tube
    from smappy_centerline.plugin import Settings
    p = tube(variable=True).points
    session = Session(Localizations(dict(x_nm=p[:, 0], y_nm=p[:, 1], z_nm=p[:, 2],
                                         frame=np.arange(len(p))*10)))
    session.set_roi(Region.rect(-10, -100, 1010, 100))
    panel = PluginPanel(CenterlineRadial, session)
    try:
        panel.form.set(Settings(dimension="3D", density=True, rejection="Auto conservative"))
        panel.run()
        wait_for_worker(panel, app)
        widget = panel.result_widget
        assert widget.table.rowCount() == 1, panel.output.toPlainText()
        assert len(panel.result.figures()) == 3
        index = [key for key, _ in widget.density_columns].index("geometry_measure")
        assert float(widget.density_table.item(0, index).text()) > 0
        monkeypatch.setattr(QFileDialog, "getExistingDirectory", lambda *a: str(tmp_path))
        widget.export()
        assert (tmp_path/"centerline_density_segments.csv").is_file()
    finally:
        if panel._window is not None:
            panel._window.close()
        panel.release()
        panel.close()


@pytest.mark.parametrize("dimension", ["2D", "3D"])
def test_result_figures_draw_without_gui(dimension):
    from matplotlib.figure import Figure
    from smappy.plugins import Context, Selection
    from smappy_centerline.plugin import Settings
    columns = {"x_nm": np.linspace(-500, 500, 120), "y_nm": np.sin(np.linspace(0, 4, 120))*20}
    if dimension == "3D":
        columns["z_nm"] = np.cos(np.linspace(0, 4, 120))*30
    ctx = Context(locs=Localizations(columns),
                  selection=Selection(np.ones(120, bool), roi=Region.rect(-600, -100, 600, 100)))
    result = CenterlineRadial().run(ctx, Settings(dimension=dimension))
    for plot in result.figures():
        figure = Figure()
        plot.draw_into(figure)
        assert len(figure.axes) == plot.panels

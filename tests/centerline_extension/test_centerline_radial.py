"""Numerical and accumulation contract; no GUI dependencies."""
import csv

import numpy as np
import pytest
from scipy.interpolate import BSpline
from scipy.spatial.transform import Rotation

from smappy_centerline.centerline import ROIResults, fit_centerline, project_curve
from smappy.locs import Localizations
from smappy.plugins import Context, Selection
from smappy_centerline.plugin import CenterlineRadial, Settings
from smappy.regions import Region
from smappy.session import Session


def cloud(dim=2, curved=False):
    x = np.repeat(np.linspace(-500, 500, 120), 4)
    offsets = np.tile([-18., -6., 6., 18.], 120)
    y = (120*(x/500)**2 if curved else np.zeros_like(x)) + offsets
    points = np.column_stack([x, y])
    if dim == 3:
        z = (70*np.sin(x/400) if curved else np.zeros_like(x)) + np.tile([-8., 8., -8., 8.], 120)
        points = np.column_stack([points, z])
    return points


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("curved", [False, True])
def test_straight_and_curved(dim, curved):
    points = cloud(dim, curved)
    result = fit_centerline(points, np.full_like(points, 10.))
    assert result.summary["N"] == len(points)
    assert 970 < result.summary["centerline_length_nm"] < 1100
    assert result.summary["coefficients"] <= 12
    if curved:
        expected = 120*(result.curve[:, 0]/500)**2
        np.testing.assert_allclose(result.curve[:, 1], expected, atol=6.)
        if dim == 3:
            np.testing.assert_allclose(result.curve[:, 2], 70*np.sin(result.curve[:, 0]/400), atol=6.)
    else:
        np.testing.assert_allclose(result.curve[:, 1:], 0., atol=1e-6)
        middle = np.abs(points[:, 0]) < 450
        np.testing.assert_allclose(result.distances[middle], np.linalg.norm(points[middle, 1:], axis=1), atol=1e-6)
    np.testing.assert_allclose(result.distances, np.linalg.norm(result.points-result.nearest, axis=1), atol=1e-10)
    assert result.summary["rms_radial_distance_nm"] == pytest.approx(np.sqrt(np.mean(result.distances**2)))
    assert result.summary["sd_nm"] == pytest.approx(np.std(result.distances))
    assert result.summary["iqr_nm"] == pytest.approx(np.diff(np.percentile(result.distances, [25, 75]))[0])
    for percentile in (30, 50, 70, 90, 95, 100):
        assert result.summary[f"R{percentile}_nm"] == pytest.approx(np.percentile(result.distances, percentile))
    assert (result.arc >= 0).all()
    assert (result.arc <= result.summary["centerline_length_nm"]).all()


@pytest.mark.parametrize("dim", [2, 3])
def test_translation_rotation_and_repeatability(dim):
    points = cloud(dim, True)
    precision = np.full_like(points, 10.)
    original = fit_centerline(points, precision)
    repeated = fit_centerline(points, precision)
    np.testing.assert_array_equal(original.distances, repeated.distances)
    translated = fit_centerline(points+100000., precision)
    np.testing.assert_allclose(original.distances, translated.distances, atol=1e-7)
    rotation = Rotation.from_rotvec([.5, -.3, .7]).as_matrix() if dim == 3 else np.array([[.8, -.6], [.6, .8]])
    rotated = fit_centerline(points @ rotation.T, precision)
    np.testing.assert_allclose(original.distances, rotated.distances, atol=1e-6)
    assert original.summary["centerline_length_nm"] == pytest.approx(rotated.summary["centerline_length_nm"])


@pytest.mark.parametrize("axis,dim", [(1, 2), (2, 3)])
def test_precision_conflict(axis, dim):
    x = np.linspace(-500, 500, 240)
    first = np.zeros((len(x), dim))
    first[:, 0] = x
    second = first.copy()
    second[:, axis] = 50
    points = np.vstack([first, second])
    equal = np.full_like(points, 20.)
    weighted = equal.copy()
    weighted[:len(x), axis] = 2.
    weighted[len(x):, axis] = 80.
    ordinary = fit_centerline(points, equal)
    precise = fit_centerline(points, weighted)
    assert np.mean(ordinary.curve[:, axis]) == pytest.approx(25., abs=1.)
    assert abs(np.mean(precise.curve[:, axis])) < 2.
    # The imprecise group still contributes its full spatial distance.
    assert np.median(precise.distances[len(x):]) > 45.


def test_outliers_are_downweighted_but_all_reported():
    points = cloud()
    points = np.vstack([points, [[0, 10000], [20, -20000], [9000, 8000]]])
    result = fit_centerline(points)
    assert result.summary["N"] == len(points)
    assert np.median(np.abs(result.curve[:, 1])) < 3.
    assert result.distances[-2] > 19000


def test_invalid_precision_and_nonfinite_coordinates():
    points = cloud()
    points[0, 0] = np.nan
    points[1, 1] = np.inf
    errors = np.full_like(points, np.nan)
    errors[2:5] = [[0, -1], [np.inf, 1e30], [1e-20, 10]]
    result = fit_centerline(points, errors)
    assert len(result.distances) == len(points)-2
    assert result.summary["excluded_nonfinite"] == 2
    assert np.isfinite(result.distances).all()
    assert result.summary["precision_replacements"] > 0


@pytest.mark.parametrize("points", [np.zeros((47, 2)), np.zeros((100, 3)),
    np.full((60, 2), np.nan), np.random.default_rng(9).normal(size=(500, 2))])
def test_insufficient_or_nonlongitudinal(points):
    with pytest.raises(ValueError):
        fit_centerline(points)


def test_projection_is_continuous_and_includes_endpoints():
    spline = BSpline([0, 0, 0, 0, 1, 1, 1, 1], [[0, 0], [100/3, 0], [200/3, 0], [100, 0]], 3)
    points = np.array([[12.345678, 7], [-3, 4], [103, 4.]])
    distance, arc, nearest, length = project_curve(points, spline)
    np.testing.assert_allclose(distance, [7, 5, 5], atol=1e-10)
    np.testing.assert_allclose(arc, [12.345678, 0, 100], atol=1e-9)
    np.testing.assert_allclose(nearest, [[12.345678, 0], [0, 0], [100, 0]], atol=1e-9)
    assert length == pytest.approx(100)


def context(points=None):
    points = cloud() if points is None else points
    locs = Localizations({f"{a}_nm": points[:, i].copy() for i, a in enumerate("xyz"[:points.shape[1]])})
    return Context(locs=locs, selection=Selection(np.ones(len(points), bool), roi=Region.rect(-1000, -1000, 1000, 1000)))


def test_accumulation_delete_and_both_csvs(tmp_path):
    plugin = CenterlineRadial()
    ctx = context()
    ctx.locs["y_nm"][0] = np.nan
    plugin.run(ctx, Settings())
    plugin.run(context(cloud(3, True)), Settings(dimension="3D"))
    plugin.run(ctx, Settings())
    assert list(plugin.results.entries) == ["ROI_001", "ROI_002", "ROI_003"]
    paths = plugin.results.export(tmp_path)
    with paths[1].open(encoding="utf-8-sig") as handle:
        all_rows = list(csv.DictReader(handle))
    assert len(all_rows) == 479*2+480
    assert len([r for r in all_rows if r["dimension"] == "3D" and r["nearest_centerline_z_nm"]]) == 480
    plugin.results.delete(["ROI_002"])
    paths = plugin.results.export(tmp_path)
    with paths[0].open(encoding="utf-8-sig") as handle:
        summaries = list(csv.DictReader(handle))
    with paths[1].open(encoding="utf-8-sig") as handle:
        distances = list(csv.DictReader(handle))
    assert len(summaries) == 2
    assert list(distances[0])[:2] == ["ROI", "distance_nm"]
    assert {r["ROI"] for r in summaries} == {r["ROI"] for r in distances} == {"ROI_001", "ROI_003"}
    for row in summaries:
        values = [float(r["distance_nm"]) for r in distances if r["ROI"] == row["ROI"]]
        assert len(values) == int(row["radial_groups"])
        assert np.mean(values) == pytest.approx(float(row["radial_mean_nm"]))
        from smappy_centerline.plugin import display_statistics
        shown = display_statistics(plugin.results.entries[row["ROI"]])
        assert "R80_nm" not in row
        for percentile in (30, 50, 70, 90, 95, 100):
            key = f"R{percentile}_nm"
            assert float(row[f"radial_r{percentile}_nm"]) == pytest.approx(np.percentile(values, percentile))
            assert float(row[f"radial_r{percentile}_nm"]) == shown[key]
        assert float(row["radial_r100_nm"]) == max(values)
    np.testing.assert_array_equal(ctx.locs["x_nm"], cloud()[:, 0])
    plugin.results.clear()
    assert not plugin.results.entries
    plugin.run(ctx, Settings())
    assert list(plugin.results.entries) == ["ROI_004"]


def test_failures_never_add_rows(monkeypatch):
    plugin = CenterlineRadial()
    with pytest.raises(ValueError, match="ROI"):
        plugin.run(Context(locs=context().locs), Settings())
    with pytest.raises(ValueError, match="z_nm"):
        plugin.run(context(), Settings(dimension="3D"))
    points = cloud(3)
    points[:, 2] = np.nan
    with pytest.raises(ValueError, match="z_nm"):
        plugin.run(context(points), Settings(dimension="3D"))
    def fail(*args, **kwargs):
        raise np.linalg.LinAlgError("test failure")
    monkeypatch.setattr(np.linalg, "solve", fail)
    with pytest.raises(ValueError, match="fitting/projection failed"):
        plugin.run(context(), Settings())
    assert not plugin.results.entries


@pytest.mark.parametrize("grouped", [False, True])
def test_snapshot_current_layer_roi_filter_grouping(grouped):
    ctx = context()
    ctx.locs.columns.update(frame=np.repeat(np.arange(len(ctx.locs)//2)*10, 2)
                            + np.tile([0, 1], len(ctx.locs)//2),
                            photons=np.arange(len(ctx.locs), dtype=float)+1)
    session = Session(ctx.locs)
    session.add_layer(like=0)
    session.render_layer_index = 1
    session.set_roi(Region.rect(-400, -100, 400, 100))
    layer = session.layers[1]
    layer.show_grouped(grouped)
    layer.set_bound("photons", 80, 390)
    table, selection = session.table(1, grouping="grouped")
    assert len(table) < len(session.locs)  # repeated detections really merged
    expected = table["x_nm"][selection.mask].copy()
    state = session.layer_configs()
    plugin = CenterlineRadial()
    frozen = plugin.prepare_context(session.context(), Settings())
    np.testing.assert_array_equal(frozen.locs["x_nm"], expected)
    assert session.layer_configs() == state
    # Changes after snapshot cannot alter data used by the worker.
    layer.set_bound("photons", 0, 1)
    session.set_roi(None)
    plugin.run(frozen, Settings())
    entry = plugin.results.entries["ROI_001"]
    assert entry.summary["N"] == len(expected)
    assert entry.summary["grouped"] is True
    assert layer.grouped == grouped
    assert entry.summary["layer"] == 1
    np.testing.assert_array_equal(entry.indices, selection.indices)


def test_empty_export(tmp_path):
    with pytest.raises(ValueError, match="no accumulated"):
        ROIResults().export(tmp_path)


@pytest.mark.parametrize("fields", [("x_err_nm", "y_err_nm"), ("xy_err_nm",)])
def test_plugin_uses_actual_precision_columns(fields):
    x = np.linspace(-500, 500, 240)
    points = np.column_stack([np.r_[x, x], np.r_[np.zeros(240), np.full(240, 50.)]])
    ctx = context(points)
    for name in fields:
        ctx.locs.columns[name] = np.r_[np.full(240, 2.), np.full(240, 80.)]
    plugin = CenterlineRadial()
    plugin.run(ctx, Settings())
    assert abs(np.mean(plugin.results.entries["ROI_001"].curve[:, 1])) < 2.


def test_smoothing_reduces_curve_bending():
    points = cloud(curved=True)
    points[:, 1] += 10*np.sin(points[:, 0]/50)
    flexible = fit_centerline(points, smoothing=.01)
    smooth = fit_centerline(points, smoothing=100.)
    def bending(entry):
        return np.sum(np.diff(entry.curve, n=2, axis=0)**2)
    assert bending(smooth) < bending(flexible)

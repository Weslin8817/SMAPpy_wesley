"""Analytical solid-tube truth and post-rejection counting contracts."""
import csv
from dataclasses import replace

import numpy as np
import pytest

from smappy_centerline.centerline import CenterlineResult
from smappy_centerline.tube_density import tube_density
from smappy.locs import Localizations
from smappy_centerline.plugin import CenterlineRadial, Settings
from smappy.regions import Region
from smappy.session import Session


def tube(curved=False, variable=False):
    # Symmetric disk quadrature: mean(r^2)=R^2/2, without random test noise.
    arc = np.repeat(np.arange(5., 1000., 10.), 48)
    theta = np.tile(np.repeat(np.arange(8)*np.pi/4, 6), 100)
    radius = (20.+20*(arc >= 500)) if variable else np.full(len(arc), 30.)
    radial = radius*np.sqrt(np.tile((np.arange(6)+.5)/6, 800))
    if curved:
        a = arc/2000
        nearest = np.column_stack([2000*np.sin(a), 2000*(1-np.cos(a)), np.zeros(len(a))])
        normal = np.column_stack([-np.sin(a), np.cos(a), np.zeros(len(a))])
        s = np.linspace(0, .5, 1001)
        curve = np.column_stack([2000*np.sin(s), 2000*(1-np.cos(s)), np.zeros(len(s))])
    else:
        nearest = np.column_stack([arc, np.zeros((len(arc), 2))])
        normal = np.tile([0., 1., 0.], (len(arc), 1))
        curve = np.column_stack([np.linspace(0, 1000, 1001), np.zeros((1001, 2))])
    points = nearest+normal*(radial*np.cos(theta))[:, None]
    points[:, 2] += radial*np.sin(theta)
    return CenterlineResult(points, np.arange(len(points)), curve, radial, arc, nearest,
                            dict(centerline_length_nm=1000., dimension="3D"))


@pytest.mark.parametrize("curved", [False, True])
@pytest.mark.parametrize("variable", [False, True])
def test_segment_radii_volume_and_nyquist_against_truth(curved, variable):
    entry = tube(curved, variable)
    record, rows, mask = tube_density(entry, segment_length_nm=250.)
    assert record["density_status"] == "exploratory", record["density_warnings"]
    radii = [20, 20, 40, 40] if variable else [30]*4
    np.testing.assert_allclose([r["radius_nm"] for r in rows], radii)
    volume = np.pi*np.sum(np.square(radii))*250/1e9
    assert record["Volume_um3"] == pytest.approx(volume)
    assert record["grouped_event_density"] == pytest.approx(mask.sum()/volume)
    assert record["grouped_event_nyquist_nm"] == pytest.approx(2000/(mask.sum()/volume)**(1/3))
    assert record["N_independent_estimated"] is None


@pytest.mark.parametrize("kind", ["sparse", "shell", "cutoff", "one_segment", "2d"])
def test_invalid_models_do_not_produce_total_density(kind):
    entry = tube()
    kwargs = {}
    if kind == "sparse":
        kwargs["min_events"] = 5000
    elif kind == "shell":
        entry.distances[:] = 30
    elif kind == "cutoff":
        kwargs["cutoff_nm"] = 32
    elif kind == "one_segment":
        kwargs["segment_length_nm"] = 2000
    else:
        entry.points = entry.points[:, :2]
    record, _, _ = tube_density(entry, **kwargs)
    assert record["Volume_um3"] is None
    assert record["grouped_event_nyquist_nm"] is None
    assert record["density_warnings"]


def test_groups_are_unique_globally_and_endpoint_counts_are_excluded():
    entry = tube()
    entry.arc[:48] = 0
    ids = np.arange(1, len(entry.points)+1)
    record, rows, support = tube_density(entry, ids)
    assert record["N_grouped_in_geometry"] == len(entry.points)-48
    assert sum(r["N_grouped"] for r in rows) == len(entry.points)-48
    assert record["N_endpoint_groups_excluded"] == 48
    assert not support[:48].any()


def test_rigid_transform_does_not_change_volume():
    entry = tube(curved=True, variable=True)
    original, _, _ = tube_density(entry)
    q, _ = np.linalg.qr(np.random.default_rng(12).normal(size=(3, 3)))
    for name in ("points", "nearest", "curve"):
        setattr(entry, name, getattr(entry, name) @ q + 1e6)
    transformed, _, _ = tube_density(entry)
    assert transformed["Volume_um3"] == pytest.approx(original["Volume_um3"])


def test_invalid_group_ids_never_become_group_count():
    entry = tube()
    with pytest.raises(ValueError, match="Group IDs"):
        tube_density(entry, np.zeros(len(entry.points)))


def test_planar_strip_area_and_nyquist():
    x = np.repeat(np.arange(5., 1000., 10.), 60)
    y = np.tile(np.linspace(-30, 30, 60), 100)
    nearest = np.column_stack([x, np.zeros_like(x)])
    points = np.column_stack([x, y])
    curve = np.column_stack([np.linspace(0, 1000, 501), np.zeros(501)])
    entry = CenterlineResult(points, np.arange(len(x)), curve, abs(y), x, nearest,
                             dict(centerline_length_nm=1000., dimension="2D"))
    record, rows, support = tube_density(entry, segment_length_nm=250.)
    area = 2*np.sqrt(3*np.mean(y*y))*1000/1e6
    assert record["density_status"] == "exploratory", record["density_warnings"]
    assert record["Area_um2"] == pytest.approx(area)
    assert record["Volume_um3"] is None
    assert record["density_units"] == "groups/um^2"
    assert record["grouped_event_density"] == pytest.approx(support.sum()/area)
    assert record["grouped_event_nyquist_nm"] == pytest.approx(2000/np.sqrt(support.sum()/area))


@pytest.mark.parametrize("radial_dim,density_dim", [("2D", "3D"), ("3D", "2D")])
def test_independent_modules_share_xyz_rejection(radial_dim, density_dim, tmp_path):
    from smappy.plugins import Context, Selection
    p = np.vstack([tube().points, [500., 0., 300.]])
    ctx = Context(locs=Localizations({f"{a}_nm": p[:, i] for i, a in enumerate("xyz")}),
                  selection=Selection(np.ones(len(p), bool), roi=Region.rect(-10, -100, 1010, 100)))
    plugin = CenterlineRadial()
    settings = Settings(dimension=radial_dim, density_dimension=density_dim,
                        rejection="Manual maximum distance", maximum_distance_nm=80.)
    result = plugin.run(ctx, settings)
    entry = plugin.results.entries["ROI_001"]
    assert result.data["dimension"] == radial_dim
    assert result.data["density_dimension"] == density_dim
    assert result.data["N_noise_excluded"] == result.data["density_noise_excluded"] == 1
    np.testing.assert_array_equal(entry.indices, entry.density_entry.indices)
    settings.radial = False
    settings.dimension = "invalid unused setting"
    only_density = plugin.run(ctx, settings)
    assert only_density.data["distance_status"] == "disabled"
    assert len(only_density.figures()) == 1
    paths = plugin.results.export(tmp_path)
    path = next(p for p in paths if p.name == f"centerline_density_{density_dim.lower()}_distances.csv")
    with path.open(encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert all(None not in row for row in rows)
    if density_dim == "2D":
        assert all(row["nearest_z_nm"] == "" for row in rows)
    settings.radial, settings.density, settings.dimension = True, False, radial_dim
    assert len(plugin.run(ctx, settings).figures()) == 2
    settings.radial = False
    with pytest.raises(ValueError, match="at least one"):
        plugin.run(ctx, settings)


@pytest.mark.parametrize("tail", [False, True])
def test_radial_gate_reports_value_direction_and_influence_without_removing_points(tail):
    entry = tube()
    if tail:
        entry.distances[0] = 1000.
    else:
        entry.distances[:] = 30.
    before = entry.distances.copy()
    record, rows, _ = tube_density(entry)
    row = rows[0]
    assert record["Volume_um3"] is None
    assert f"Q={row['radial_moment_ratio']:.4f}" in row["warning"]
    assert ("above 2.5" if tail else "below 1.12") in row["warning"]
    assert f"N_grouped={row['N_grouped']}" in row["warning"]
    assert "not a statistical test" in row["warning"]
    d2 = before[entry.arc < row["arc_end_nm"]]**2
    expected = (d2.mean() - np.sort(d2)[:-1].mean())/d2.mean()
    assert row["area_change_without_farthest"] == pytest.approx(expected)
    np.testing.assert_array_equal(entry.distances, before)


def test_broad_tail_can_warn_but_strict_policy_and_influence_guard_still_block():
    entry = tube()
    entry.distances[0] = 200.
    before = entry.distances.copy()
    record, rows, _ = tube_density(entry)
    assert record["density_status"] == "exploratory_warning"
    assert rows[0]["radial_moment_ratio"] > 2.5
    assert rows[0]["area_change_without_farthest"] < .25
    assert record["grouped_event_density"] > 0
    assert "model-dependent" in record["density_warnings"]
    from matplotlib.figure import Figure
    from smappy_centerline.plugin import plots
    entry.summary.update(record)
    entry.density_segments = rows
    ax = Figure().subplots()
    plots(entry)["Segment radii"](ax)
    assert "Model warning" in ax.get_legend_handles_labels()[1]
    strict, _, _ = tube_density(entry, radial_policy="Strict")
    assert strict["Volume_um3"] is None
    unstable, _, _ = tube_density(entry, sensitivity_limit=.01)
    assert unstable["Volume_um3"] is None
    assert "influence" in unstable["density_warnings"]
    np.testing.assert_array_equal(entry.distances, before)


def test_radial_policy_rejects_unknown_option():
    with pytest.raises(ValueError, match="policy"):
        tube_density(tube(), radial_policy="Ignore everything")


@pytest.mark.parametrize("removed", [None, 0, 32])
def test_cutoff_proximity_uses_actual_rejection_count(removed):
    entry = tube()
    record, rows, _ = tube_density(entry, cutoff_nm=35., n_rejected=removed)
    # A clean 30 nm-radius disk is close to 80% of 35 nm, even though
    # every sampled point lies inside the cutoff. Proximity is not truncation.
    assert entry.distances.max() < 35.
    if removed == 0:
        assert record["density_status"] == "exploratory_warning"
        assert record["Volume_um3"] == pytest.approx(np.pi*30**2*1000/1e9)
        assert "0 grouped points removed" in rows[0]["warning"]
    else:
        assert record["Volume_um3"] is None
        assert "may be censored" in rows[0]["warning"]
    assert "radius=30.00 nm" in rows[0]["warning"]
    assert "effective cutoff=35.00 nm" in rows[0]["warning"]
    assert "trigger=28.00 nm" in rows[0]["warning"]
    if removed:
        assert "32/4832" in rows[0]["warning"]


@pytest.mark.parametrize("removed", [-1, .5, np.nan])
def test_invalid_rejection_count_is_refused(removed):
    with pytest.raises(ValueError, match="Rejected count"):
        tube_density(tube(), n_rejected=removed)


def test_density_is_always_grouped_3d_independent_of_distance_dimension(tmp_path):
    from smappy.plugins import Context, Selection
    p = np.vstack([tube().points, [500., 0., 200.]])
    locs = Localizations({f"{a}_nm": p[:, i] for i, a in enumerate("xyz")})
    ctx = Context(locs=locs, selection=Selection(np.ones(len(p), bool), roi=Region.rect(-10, -100, 1010, 100)))
    results = []
    for dim in ("2D", "3D"):
        plugin = CenterlineRadial()
        result = plugin.run(ctx, Settings(dimension=dim, density=True,
            rejection="Manual maximum distance", maximum_distance_nm=80.))
        results.append(result.data)
        assert result.data["dimension"] == dim
        assert result.data["density_dimension"] == "3D"
        assert result.data["density_noise_excluded"] == 1
        assert "localization_density" not in result.data
        assert "effective_labeling_density" not in result.data
        paths = plugin.results.export(tmp_path)
        with (tmp_path/"centerline_density_3d_distances.csv").open(encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == result.data["N_grouped"]
        assert sum(r["in_density_domain"] == "True" for r in rows) == result.data["N_grouped_in_geometry"]
    assert results[0]["N_noise_excluded"] == 1
    for key in ("Volume_um3", "N_grouped", "N_grouped_in_geometry", "grouped_event_density", "grouped_event_nyquist_nm"):
        assert results[0][key] == pytest.approx(results[1][key])


def test_missing_z_does_not_invent_volume_or_destroy_2d_distances():
    from smappy.plugins import Context, Selection
    p = tube().points
    ctx = Context(locs=Localizations(dict(x_nm=p[:, 0], y_nm=p[:, 1])),
                  selection=Selection(np.ones(len(p), bool), roi=Region.rect(-10, -100, 1010, 100)))
    result = CenterlineRadial().run(ctx, Settings(density=True, dimension="2D"))
    assert result.data["dimension"] == "2D"
    assert result.data["density_status"] == "unavailable"
    assert result.data["grouped_event_nyquist_nm"] is None
    assert "z_nm" in result.data["density_warnings"]


def test_failed_xy_centerline_does_not_block_density_along_z(tmp_path):
    from smappy.plugins import Context, Selection
    from matplotlib.figure import Figure
    p = tube().points[:, [1, 2, 0]]
    ctx = Context(locs=Localizations({f"{a}_nm": p[:, i] for i, a in enumerate("xyz")}),
                  selection=Selection(np.ones(len(p), bool), roi=Region.rect(-100, -100, 100, 100)))
    plugin = CenterlineRadial()
    result = plugin.run(ctx, Settings(density=True, dimension="2D"))
    assert result.data["distance_status"] == "unavailable"
    assert result.data["density_status"] == "exploratory", result.text
    assert result.data["grouped_event_density"] > 0
    assert result.data["N"] == 0
    for plot in result.figures():
        plot.draw_into(Figure())
    paths = plugin.results.export(tmp_path)
    with (tmp_path/"centerline_density_3d_distances.csv").open(encoding="utf-8-sig") as handle:
        assert len(list(csv.DictReader(handle))) == result.data["N_grouped"]


def test_plugin_snapshot_rejection_counts_and_csv(tmp_path):
    entry = tube(variable=True)
    points = np.vstack([entry.points, [[100, 200, 0], [900, -200, 0]]]) + [5000, 5000, 0]
    locs = Localizations(dict(x_nm=points[:, 0], y_nm=points[:, 1], z_nm=points[:, 2],
                             frame=np.arange(len(points))*10, photons=np.full(len(points), 1000.)))
    session = Session(locs)
    session.show_grouped(0, True)
    session.roi = Region.rect(4900, 4700, 6100, 5300)
    settings = Settings(dimension="3D", density=True, rejection="Manual maximum distance",
                        maximum_distance_nm=80., segment_length_nm=250.)
    plugin = CenterlineRadial()
    before = {k: v.copy() for k, v in session.locs.columns.items()}
    snapshot = plugin.prepare_context(session.context(), settings)
    expected = snapshot.centerline_group_ids.copy()
    # Worker must use frozen ROI, group IDs and coordinates.
    session.roi = None
    session.set_group_settings(replace(session.group_settings, dt=2))
    result = plugin.run(snapshot, settings)
    stored = plugin.results.entries["ROI_001"]
    assert result.data["N_noise_excluded"] == 2
    assert result.data["N_grouped"] == len(entry.points)
    assert result.data["N_grouped"] == len(np.unique(expected[stored.indices]))
    assert result.data["Volume_um3"] is not None, result.text
    assert result.data["N_grouped_in_geometry"] == len(np.unique(expected[stored.indices][stored.density_support]))
    paths = plugin.results.export(tmp_path)
    assert len(paths) == 5
    with paths[2].open(encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) > 1
    assert float(rows[-1]["radius_nm"]) > 1.5*float(rows[0]["radius_nm"])
    from matplotlib.figure import Figure
    assert len(result.figures()) == 3
    for plot in result.figures():
        figure = Figure()
        plot.draw_into(figure)
        assert len(figure.axes) == plot.panels
    with paths[0].open(encoding="utf-8-sig") as handle:
        assert next(csv.DictReader(handle))["density_independent_labels"] == ""
    for key, value in before.items():
        np.testing.assert_array_equal(session.locs[key], value)
    # Exporting an entry without density into the same folder clears the old
    # segment rows rather than leaving a contradictory result beside it.
    del stored.density_segments
    plugin.results.export(tmp_path)
    with paths[2].open(encoding="utf-8-sig") as handle:
        assert list(csv.DictReader(handle)) == []


def test_stale_groups_rebuilt_and_grouped_filter_applies():
    entry = tube()
    p = entry.points
    locs = Localizations(dict(x_nm=p[:, 0], y_nm=p[:, 1], z_nm=p[:, 2],
                             photons=np.tile([100., 1000.], len(p)//2), frame=np.arange(len(p))*10))
    s = Session(locs)
    s.show_grouped(0, True)
    s.roi = Region.rect(-1, -100, 1001, 100)
    s.layers[0].set_bound("photons", 500., None)
    s.set_group_settings(replace(s.group_settings, dt=2))
    plugin = CenterlineRadial()
    snapshot = plugin.prepare_context(s.context(), Settings(density=True))
    assert snapshot.centerline_group_ids is not None
    assert snapshot.centerline_provenance["grouped"] is True
    assert np.all(snapshot.locs["photons"] == 1000.)
    assert len(snapshot.locs) > 0  # rebuilt on the scratch session

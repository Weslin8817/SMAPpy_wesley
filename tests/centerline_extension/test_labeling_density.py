"""Regression tests for reusable geometry and Nyquist numerical helpers."""
import numpy as np
import pytest

from smappy_centerline.labeling_density import density_nyquist, reconstruct


def grid(d=2, n=6):
    return np.stack(np.meshgrid(*([np.linspace(0, 1000, n)] * d)), axis=-1).reshape(-1, d)


@pytest.mark.parametrize("d", [2, 3])
@pytest.mark.parametrize("method", ["Alpha Shape", "Convex Hull"])
def test_unit_square_and_cube_have_unit_measure(d, method):
    result = reconstruct(grid(d), method)
    assert result["measure_um"] == pytest.approx(1)
    translated = reconstruct(grid(d) + 1e8, method)
    assert translated["measure_um"] == pytest.approx(1)
    scaled = reconstruct(grid(d) * 2, method)
    assert scaled["measure_um"] == pytest.approx(2 ** d)


@pytest.mark.parametrize("d,count,expected", [(2, 100, 200), (3, 1000, 200)])
def test_nyquist_converts_micrometres_to_nanometres(d, count, expected):
    density, nyquist = density_nyquist(count, 1, d)
    assert density == count
    assert nyquist == pytest.approx(expected)


@pytest.mark.parametrize("points", [np.empty((0, 2)), np.zeros((100, 2)),
    np.column_stack([np.arange(40), np.arange(40)]),
    np.column_stack([grid(), np.zeros(36)]), np.full((40, 3), np.nan), grid()[:4]])
def test_invalid_geometry_has_no_measure(points):
    result = reconstruct(points)
    assert result["measure_um"] is None
    assert result["warnings"]


def test_duplicates_do_not_change_geometry():
    assert reconstruct(np.repeat(grid(), 5, axis=0))["measure_um"] == pytest.approx(1)


def test_radius_sensitivity_rejects_successful_but_unstable_shape():
    result = reconstruct(grid(), radius_nm=150)
    assert result["measure_um"] is None
    assert result["sensitivity"] > .25
    assert any("radius-sensitive" in message for message in result["warnings"])
    assert len(result["radius_scan"]) == 1  # manual radii are never increased
    assert not result["support_mask"].any()


def test_tiny_radius_and_outliers_are_not_silently_filled_or_removed():
    assert reconstruct(grid(), radius_nm=1)["measure_um"] is None
    result = reconstruct(np.vstack([grid(), [1e6, 1e6]]))
    assert result["measure_um"] == pytest.approx(1)
    assert result["coverage"] == pytest.approx(36 / 37)
    np.testing.assert_array_equal(result["support_mask"], [True] * 36 + [False])
    assert result["n_nn_flagged"] == 1
    assert any("Partial support" in message for message in result["warnings"])


def test_alpha_retains_holes_and_nonconvexity():
    p = grid(n=21)
    ring = p[(np.abs(p[:, 0] - 500) >= 250) | (np.abs(p[:, 1] - 500) >= 250)]
    a = reconstruct(ring, radius_nm=60)
    assert .7 < a["measure_um"] < .9
    assert reconstruct(ring, "Convex Hull")["measure_um"] == pytest.approx(1)
    corner = p[(p[:, 0] <= 400) | (p[:, 1] <= 400)]
    assert reconstruct(corner, radius_nm=60)["measure_um"] < .7










def test_3d_alpha_preserves_a_cavity_and_disconnected_components():
    p = grid(3, 11)
    shell = p[np.any(np.abs(p - 500) >= 300, axis=1)]
    result = reconstruct(shell, radius_nm=120)
    assert .8 < result["measure_um"] < .95
    separated = np.vstack([grid(3), grid(3) + [3000, 0, 0]])
    assert reconstruct(separated)["measure_um"] == pytest.approx(2)
    assert reconstruct(separated, "Convex Hull")["measure_um"] == pytest.approx(4)


def test_qhull_failure_is_reported_without_a_fallback(monkeypatch):
    from scipy import spatial
    def fail(*args, **kwargs):
        raise spatial.QhullError("test failure")
    monkeypatch.setattr(spatial, "Delaunay", fail)
    result = reconstruct(grid())
    assert result["measure_um"] is None
    assert "no fallback" in result["warnings"][-1]








def repeat_lattice():
    # 259 spatial sites, 129 near-repeat events: the original hard NN gate
    # rejects this legitimate density variation even in Convex Hull mode.
    sites = np.stack(np.meshgrid(np.linspace(0, 300, 7), np.linspace(0, 700, 37)),
                     axis=-1).reshape(-1, 2)
    return np.vstack([sites, sites[:129] + [1, .2]])


@pytest.mark.parametrize("method", ["Convex Hull", "Alpha Shape"])
def test_close_repeat_events_do_not_veto_a_well_sampled_structure(method):
    p = repeat_lattice()
    before = p.copy()
    result = reconstruct(p, method)
    assert result["nn_max_median_ratio"] > 10
    assert result["n_nn_flagged"] > 0
    assert result["measure_um"] == pytest.approx(.21, rel=.01)
    assert result["support_mask"].all()
    assert result["event_coverage"] == 1
    np.testing.assert_array_equal(p, before)
    if method == "Alpha Shape":
        assert result["auto_seed_nm"] > 2 * result["nn_median_nm"]
        assert len(result["radius_scan"]) > 1
        assert result["radius_scan"][-1]["accepted"]
        assert result["radius_nm"] <= result["auto_max_nm"]
        assert result["components"] == 1


@pytest.mark.parametrize("d", [2, 3])
def test_convex_hull_reports_outlier_inflation_without_removing_points(d):
    from scipy.spatial import ConvexHull
    p = np.vstack([grid(d), np.full((1, d), 1e6)])
    result = reconstruct(p, "Convex Hull")
    assert result["measure_um"] == pytest.approx(ConvexHull(p).volume / 1000 ** d)
    assert result["support_mask"].all()
    assert result["n_nn_flagged"] == 1
    assert any("inflated" in message for message in result["warnings"])


@pytest.mark.parametrize("d", [2, 3])
def test_partial_support_preserves_original_order_and_duplicate_multiplicity(d):
    p = np.vstack([np.full((1, d), 1e6), grid(d), grid(d)[-1:]])
    result = reconstruct(p)
    expected = np.ones(len(p), dtype=bool)
    expected[0] = False
    np.testing.assert_array_equal(result["support_mask"], expected)
    assert result["event_coverage"] == pytest.approx((len(p) - 1) / len(p))
    assert result["measure_um"] == pytest.approx(1)




def test_unique_point_coverage_alone_cannot_hide_many_unsupported_events():
    p = np.vstack([grid(), np.tile([1e6, 1e6], (10, 1))])
    result = reconstruct(p, radius_nm=400)
    assert result["coverage"] > .95
    assert result["event_coverage"] < .95
    assert result["measure_um"] is None
    assert result["candidate_measure_um"] == pytest.approx(1)
    assert not result["support_mask"].any()


def test_strict_coverage_remains_available_without_silently_increasing_manual_radius():
    p = np.vstack([grid(), [1e6, 1e6]])
    result = reconstruct(p, radius_nm=400, min_coverage=1)
    assert result["measure_um"] is None
    assert result["radius_nm"] == pytest.approx(400)
    assert len(result["radius_scan"]) == 1


def test_automatic_search_is_bounded_and_does_not_relax_sensitivity():
    # No radius in the bounded search passes this deliberately tight tolerance.
    p = np.random.default_rng(6).uniform(0, 1000, (100, 2))
    result = reconstruct(p, sensitivity_limit=.00001)
    assert result["measure_um"] is None
    assert len(result["radius_scan"]) == 6
    assert all(not trial["accepted"] for trial in result["radius_scan"])
    assert result["radius_nm"] == pytest.approx(result["auto_seed_nm"] * 4)


@pytest.mark.parametrize("coverage", [0, 1.1, np.nan])
def test_invalid_coverage_is_rejected(coverage):
    with pytest.raises(ValueError, match="invalid alpha"):
        reconstruct(grid(), min_coverage=coverage)

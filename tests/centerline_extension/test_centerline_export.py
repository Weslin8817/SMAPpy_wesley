import csv
import json

from smappy_centerline.centerline_export import FIELDS, summary_row, write_dictionary


def test_fixed_schema_status_empty_values_and_precision(tmp_path):
    disabled = summary_row(dict(ROI="ROI_001", radial_enabled=False, density_enabled=False,
                                N=0, noise_rejection="None", noise_bins=16))
    failed = summary_row(dict(ROI="ROI_002", N=0, distance_status="unavailable",
                             dimension="2D", distance_warning="fit failed", density_enabled=True,
                             density_status="unavailable", density_dimension="3D",
                             density_bins=24, density_smoothing=2.))
    success = summary_row(dict(ROI="ROI_003", N=100, R50_nm=1.234567890123456,
                              median_radial_distance_nm=1.234567890123456,
                              density_enabled=True, density_status="exploratory_warning",
                              density_dimension="2D", Area_um2=.5, Volume_um3=123.,
                              noise_rejection="Manual maximum distance", N_noise_excluded=0,
                              future_diagnostic=[1, 2]))
    assert tuple(disabled) == tuple(failed) == tuple(success) == FIELDS
    assert len(FIELDS) == len(set(FIELDS))
    assert disabled["radial_status"] == disabled["density_status"] == "disabled"
    assert disabled["radial_groups"] is None
    assert disabled["noise_bins"] is None
    assert failed["radial_groups"] is None
    assert failed["radial_warning"] == "fit failed"
    assert failed["density_bins"] == 24
    assert success["radial_r50_nm"] == 1.234567890123456
    assert success["noise_removed_groups"] == 0
    assert success["density_status"] == "warning"
    assert success["density_area_um2"] == .5
    assert success["density_volume_um3"] is None
    assert success["density_units"] == "groups/um^2"
    assert "Exploratory" in success["density_interpretation"]
    assert json.loads(success["additional_metadata_json"]) == {"future_diagnostic": [1, 2]}
    assert "median_radial_distance_nm" not in FIELDS
    path = tmp_path / "fields.csv"
    write_dictionary(path)
    with path.open(encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["field"] for row in rows] == list(FIELDS)
    assert all(row["description_zh"] and row["empty_meaning"] for row in rows)
    assert all(row["export_format_version"] == "2" for row in rows)

"""Stable, versioned CSV summary schema, independent of GUI presentation."""
import csv
import json

EXPORT_VERSION = "2"
COLUMNS = []


def _group(name, rows):
    for line in rows.strip().splitlines():
        field, source, description, unit, dimension = line.split("|")
        COLUMNS.append((field, source or field, name, description, unit, dimension))


_group("ROI 与降噪", """
ROI||ROI唯一编号||
roi_center_x_nm||ROI包围框中心X|nm|XY
roi_center_y_nm||ROI包围框中心Y|nm|XY
noise_mode|noise_rejection|共享降噪模式||
noise_input_groups|N_before_rejection|降噪前输入分组数|groups|
noise_invalid_xyz_groups|noise_excluded_nonfinite|共享降噪排除的无效XYZ数；None模式不适用|groups|XYZ
noise_removed_groups|N_noise_excluded|阈值排除的分组数，不含无效坐标|groups|XYZ
noise_retained_groups|N_after_rejection|降噪后分组数；None模式等于输入数|groups|
noise_cutoff_nm||实际筛选距离阈值|nm|XYZ
""")
_group("Radial Distance 核心结果", """
radial_status||模块状态：disabled/unavailable/complete/warning||
radial_dimension|dimension|距离分析维度||
radial_groups|N|参与距离分析的有效分组数|groups|
radial_length_nm|centerline_length_nm|中心线长度|nm|
radial_mean_nm|mean_radial_distance_nm|平均径向距离|nm|
""")
for percentile in (30, 50, 70, 90, 95, 100):
    _group("Radial Distance 核心结果",
           f"radial_r{percentile}_nm|R{percentile}_nm|径向距离P{percentile}；R50为中位数|nm|")
_group("Labeling Density 核心结果", """
density_status||模块状态：disabled/unavailable/complete/warning；完成仍为探索性估计||
density_dimension||密度分析维度||
density_retained_groups|N_grouped|密度分析有效分组数|groups|
density_geometry_groups|N_grouped_in_geometry|几何域内用于密度分子的分组数|groups|
density_area_um2|Area_um2|平面填充条带面积|um^2|2D
density_volume_um3|Volume_um3|填充圆管体积|um^3|3D
density_group_density|grouped_event_density|分组事件密度，非独立分子密度|见density_units|
density_units||分组事件密度单位：groups/um^2或groups/um^3||
density_nyquist_nm|grouped_event_nyquist_nm|探索性Nyquist估计，非实测分辨率|nm|
""")
_group("补充统计与诊断", """
radial_rms_nm|rms_radial_distance_nm|径向距离均方根|nm|
radial_sd_nm|sd_nm|径向距离总体标准差|nm|
radial_iqr_nm|iqr_nm|径向距离四分位距|nm|
radial_invalid_groups|excluded_nonfinite|当前距离维度的无效坐标数|groups|
radial_precision_replacements|precision_replacements|定位精度回退次数||
radial_populated_bins|populated_bins|实际有效拟合分箱数||
radial_coefficients|coefficients|样条系数数量||
radial_refit_shift_nm|rejection_refit_shift_nm|可用时记录降噪重拟合最大位移|nm|
density_length_nm||密度中心线长度|nm|
density_invalid_groups|density_excluded_nonfinite|当前密度维度的无效坐标数|groups|
density_endpoint_excluded_groups|N_endpoint_groups_excluded|因端点投影排除的分组数|groups|
density_segment_sensitivity||分段尺度变化的相对面积或体积变化|fraction|
density_area_coarse_um2||粗分段面积|um^2|2D
density_volume_coarse_um3||粗分段体积|um^3|3D
density_effective_cutoff_nm||考虑重拟合位移后的诊断阈值，不同于实际筛选阈值|nm|
density_independent_labels|N_independent_estimated|独立标签数；未校准时未知||
""")
_group("分析参数", """
noise_bins||共享降噪拟合分箱设置||
noise_smoothing||共享降噪平滑设置||
radial_bins|bins|距离拟合分箱设置||
radial_smoothing|smoothing|距离拟合平滑设置||
density_bins||密度独立拟合分箱设置||
density_smoothing||密度独立拟合平滑设置||
density_segment_length_nm||最大分段长度|nm|
density_min_events||每段最少分组数|groups|
density_sensitivity_limit||尺度敏感性上限|fraction|
density_radial_policy||径向质量控制策略||
density_point_influence_limit||单点影响上限|fraction|
""")
_group("警告与来源", """
radial_warning|distance_warning|距离分析失败或警告原因||
density_warning|density_warnings|密度质量说明或失败原因||
density_model||密度几何模型||
density_scope||计数和几何域定义||
density_interpretation||探索性结果解释及限制||
density_independent_reason|independent_reason|不能推断独立标签数的原因||
roi_bounds|source_roi|捕获的ROI几何描述||
roi_layer|layer|来源图层索引||
roi_layer_name|layer_name|来源图层名称||
roi_source_path|source_path|来源文件路径||
roi_selection|selection|选择描述||
roi_grouped|grouped|是否使用分组表||
roi_render_grouped|render_grouped|Render分组显示状态||
roi_filter_ranges|filter_ranges|过滤范围JSON||
roi_slab|slab|切片限制||
roi_grouping_parameters|grouping_parameters|分组参数JSON||
roi_grouping_reason|grouping_reason|分组计算说明||
radial_algorithm|algorithm|距离拟合算法||
plugin_version||插件版本||
export_format_version||CSV格式版本，当前为2||
additional_metadata_json||未纳入标准列的扩展元数据JSON，防止丢失信息||
""")
FIELDS = tuple(column[0] for column in COLUMNS)
REDUNDANT = {"median_radial_distance_nm", "radial_enabled", "density_enabled",
             "distance_status", "density_noise_cutoff_nm", "density_noise_excluded",
             "density_rejected_count", "density_pre_rejection_groups", "nyquist_units"}


def summary_row(values):
    row = {field: values.get(source) for field, source, *_ in COLUMNS}
    radial = values.get("radial_enabled", True)
    density = values.get("density_enabled", "density_status" in values)
    raw_radial = values.get("distance_status", "complete")
    raw_density = values.get("density_status", "unavailable")
    row["radial_status"] = "disabled" if not radial else raw_radial
    row["density_status"] = ("disabled" if not density else
        {"exploratory": "complete", "exploratory_warning": "warning"}.get(raw_density, raw_density))
    for prefix, enabled in (("radial_", radial), ("density_", density)):
        if not enabled:
            for field in FIELDS:
                if field.startswith(prefix) and field != prefix + "status":
                    row[field] = None
    if radial and raw_radial == "unavailable":
        for field, _, group, *_ in COLUMNS:
            if field.startswith("radial_") and group in ("Radial Distance 核心结果", "补充统计与诊断"):
                if field not in ("radial_status", "radial_dimension"):
                    row[field] = None
    if values.get("noise_rejection") == "None":
        for field in ("noise_invalid_xyz_groups", "noise_cutoff_nm", "noise_bins", "noise_smoothing"):
            row[field] = None
    if row["density_dimension"] == "2D":
        row["density_volume_um3"] = row["density_volume_coarse_um3"] = None
    elif row["density_dimension"] == "3D":
        row["density_area_um2"] = row["density_area_coarse_um2"] = None
    if density:
        dim = row["density_dimension"]
        row["density_units"] = row["density_units"] or {"2D": "groups/um^2", "3D": "groups/um^3"}.get(dim)
        row["density_interpretation"] = row["density_interpretation"] or "Exploratory sampling estimate; not measured image resolution or independent-label density."
    row["export_format_version"] = EXPORT_VERSION
    known = {column[1] for column in COLUMNS} | REDUNDANT
    extra = {key: value for key, value in values.items() if key not in known}
    row["additional_metadata_json"] = json.dumps(extra, ensure_ascii=False, default=str) if extra else None
    return row


def write_dictionary(path):
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["field", "group", "description_zh", "unit", "dimension", "empty_meaning", "legacy_field", "export_format_version"])
        for field, source, group, description, unit, dimension in COLUMNS:
            writer.writerow([field, group, description, unit, dimension or "按模块维度",
                "空白表示未启用、不可用、未记录或不适用；查看模块状态和警告。0为实际零值。",
                source if source != field else "", EXPORT_VERSION])

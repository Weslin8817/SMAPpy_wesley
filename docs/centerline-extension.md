# Centerline / Labeling Density 移植记录

## 范围和来源

- 只读来源：`G:/SMAPpy_wes`，`main`，commit
  `ed55cbf9433edc2fe028561ad13a97f400efcf55`，集成插件算法版本 `9`。
- 目标：`G:/SMAPpy_wesley`，分支
  `feature/Radial-diatance-&-labeling-density`；官方基线 commit
  `8839434f7286cd1f1d9f8720ba95171a1256b205`，SMAPpy `0.5.10`。
- 未切换、合并、重置、提交或推送分支；源目录未执行 Python、构建或测试，未写入文件。
- 目标 pyproject 明确要求 Python >=3.12，因此使用项目已有的 `.venv` Python
  3.12.16。未新增第三方依赖。用现有 MSVC 构建官方五个 C++ 扩展。

## 文件清单

| 新增文件 | 作用 |
| --- | --- |
| `src/smappy_centerline/__init__.py` | 独立命名空间，不向 `smappy` 注入模块 |
| `src/smappy_centerline/centerline.py` | 原版拟合、连续投影、噪声排除、统计、ROI 累计及导出组织 |
| `src/smappy_centerline/centerline_export.py` | 原版固定 CSV 字段、版本 2 导出、字段中文说明 |
| `src/smappy_centerline/tube_density.py` | 原版 filled-strip / filled-tube 分段、计数、质量控制和稳定性诊断 |
| `src/smappy_centerline/labeling_density.py` | 原版 Nyquist 与几何工具；保留旧模块及其测试 |
| `src/smappy_centerline/plugin.py` | 集成参数、原版计算顺序、绘图、独立模块失败处理，以及新版数据快照适配 |
| `src/smappy_centerline/gui.py` | 隔离的三张联动结果表、详情、删除、清空、CSV 操作及表单分节 |
| `user_plugins/Analysis/Measure/centerline_radial.py` | 官方插件扫描器可发现的薄入口 |
| `user_plugins/Analysis/Measure/centerline_radial.md` | 迁移的完整算法和使用帮助，保留版本 9 定义 |
| `user_plugins/README.md` | 安装入口、启动与操作说明 |
| `scripts/run_centerline.py` | 自动配置仓库内插件目录的启动器，配置也仅写入目标仓库 |
| `tests/centerline_extension/test_centerline_radial.py` | 迁移的拟合、精度、统计、ROI 累计与 CSV 测试 |
| `tests/centerline_extension/test_centerline_export.py` | 迁移的字段、空值、数值精度与导出版本测试 |
| `tests/centerline_extension/test_labeling_density.py` | 迁移的密度几何/Nyquist 测试 |
| `tests/centerline_extension/test_tube_density.py` | 迁移的分段真值、QC、独立维度、共享拒绝与模块失败测试 |
| `tests/centerline_extension/test_gui_centerline_radial.py` | 迁移的真实 Qt 面板、工作线程、表格交互和绘图测试 |
| `tests/centerline_extension/test_compatibility.py` | 新增发现/加载/帮助/运行、只读快照、slab、图片层拒绝和线程测试 |
| `tests/centerline_extension/source_manifest.json` | 来源 commits、算法文件 SHA256、计算流程 AST 摘要及默认值 |
| `docs/centerline-extension.md` | 本记录 |

原版四个算法/导出模块在移植时与源文件逐字节相同。回归测试检查规范化换行后的
SHA256，以兼容 Git 的 Windows CRLF 转换；同时锁定 `run`、`plots`、`active`、
ROI 坐标和统计函数的 AST 及全部参数默认值。改变这些定义时须明确更新来源基准，
不可通过降低数值容差绕过检查。

## 隔离程度及四处官方修改

没有复制旧版 GUI 框架、插件管理器或数据核心。没有修改任何官方内置插件、
拟合器、渲染算法、漂移校正或 HDF5 读写实现。额外插件根目录负责发现，
`src/` 下独立包由现有 setuptools 包发现机制处理，不需修改项目依赖或打包配置。

新版 `Result` 支持文字和图形，但不提供三表累计、多选删除和自定义 CSV 操作；
`Plugin.window` 需要由主程序注册窗口，不能单靠外部插件完成。
因此添加以下默认无操作的通用接入点，专用逻辑全部留在扩展包内，未使用全局 monkey patch。

| 修改的官方文件 | 必要改动 |
| --- | --- |
| `src/smappy/plugins/__init__.py` | 增加默认返回原 context 的 `prepare_context(ctx, settings)`、默认返回 None 的 `create_result_widget(panel)`，并自动发现源码仓库的 `user_plugins` |
| `src/smappy/gui/plugin_panel.py` | 面板建好后装载可选组件；GUI 线程启动 worker 前调用快照钩子并处理异常；成功展示结果后发出 `result_updated` 信号 |
| `src/smappy/session.py` | 为当前 Render 控件所选 layer 提供 `render_layer_index = 0` 初值 |
| `src/smappy/gui/render_tab.py` | `_bind_layer` 时更新上述索引，以使用用户选中的 layer，而非固定 layer 0 |

四个官方文件保留小型接入点。普通启动另在插件 roots() 中自动扫描当前源码仓库的 user_plugins 目录；路径由已导入的包定位，不依赖工作目录，不写个人配置。用户目录及显式配置目录仍有更高覆盖优先级。

## 数据适配及行为

旧版的 `layer.state.sets`、`grouped_stale`、每层 `group_settings` 在新版已被替换。
新版 grouping 使用 Session 共享的 `Tables` 和 `Session.group_settings`，filter 使用
`Layer.ranges`。本插件遵循新版当前 grouping 参数，不引入第二套分组算法。

1. 在 GUI 线程验证当前 Render ROI 和选中层；图片层明确拒绝。
2. 使用官方 `Session.scratch()`，在副本上调用 `table(index, grouping="grouped")`。
   它复用可用的分组结果；缺失或过期时在副本重建，不给原始表附加或更改分组列。
3. 官方 selection 负责 layer/file filter、Render axes 的 ROI 和适用 slab 筛选。
   将选中行、原分组表行号、group IDs、精度列、ROI 和元数据复制到独立 `Context`。
4. worker 只收到无 Session 引用的快照；之后拖动 ROI、修改 filter 或改变输入数组
   不会改变本次计算。插件不返回 locs、ROI、bounds 或渲染状态修改。

`Plugin.grouping = "ungrouped"` 只防止通用 Context 在快照前对原表隐式执行 grouping；
实际计算仍始终使用副本上的 **grouped** 表，保持原版分析口径。
新版 GUI 自身的正常渲染和后台 merging 仍由官方代码管理。

保留的数值/交互约定：

- 默认 Radial 开启、2D；Density 开启、3D；两个模块拟合参数独立。
- 手动与自动噪声排除均在 XYZ 中共享执行一次，再分别分析；保留原阈值、
  顺序、refit 位移处理和统计定义。
- 一个模块失败时另一模块仍可输出；保持原版错误类型和两模块均失败时的处理。
- R30、R50、R70、R90、R95、R100 都是径向距离，单位 nm；不转换为直径。
- 保留定位精度权重、稳健中心线、每段三次曲线的连续最近点求解。
- 2D filled-strip / 3D filled-tube 是连续填充、单条细长结构模型；保留稀疏段、
  径向尾部、影响点、截断、弯曲、非局部重叠和分段敏感性诊断。
- grouped event 不等于 independent labeling site。没有独立位点校准时，
  independent-label 密度与 Nyquist 保持不可用；不把探索性输出表述为真实标记密度。
- 每次一个 ROI，结果只在该插件实例内累计；选择、删除、清空不改变数据。
  编号不复用。保留三张联动表、详情、2D/3D 拟合图、径向分布和密度分段图。
- 使用官方 `Result` / `Plot` / ResultWindow 绘图；只有累计表格使用专用组件。

CSV 文件名、顺序、字段、单位、BOM 编码、空值语义及版本保持原版：
`centerline_radial_summary.csv`、`centerline_radial_all_distances.csv`、
`centerline_density_segments.csv`、`centerline_density_2d_distances.csv`、
`centerline_density_3d_distances.csv`、`centerline_summary_fields.csv`。
与原版一样，按存在的密度维度和分段结果写出相应文件，空/旧文件处理由原导出器负责。

## 启动和使用

普通启动会自动发现此插件。更新代码后须关闭旧 SMAPpy 进程，再在目标仓库根目录运行：

```powershell
.\.venv\python.exe -m smappy.cli.gui
```

以下专用启动器仍可用于独立配置：

```powershell
.\.venv\python.exe scripts/run_centerline.py
```

此方式自动注册仓库内 `user_plugins`，使用 `build/centerline-config` 中的独立配置。
本源码仓库沿用正常配置时无需手动注册。其他安装位置可在 **Preferences → Plugin folders** 添加
`G:/SMAPpy_wesley/user_plugins`，之后按原方式启动即可；无需复制插件到个人目录。
移动仓库后须更新 Preferences 中保存的绝对路径，启动器则自动使用自己的仓库路径。

在 Render 选择 localization layer、绘制一个 ROI，再打开
**Plugins → Analysis → Measure → Radial distance and density**。
配置三组参数后点击 **Analyze/Add current ROI**；换 ROI 后再次点击即可累计。
通过结果表选择查看、删除或清空，点击 **Export CSV** 导出所有累计 ROI。
结果不会自动写入 HDF5 或 GUI 配置，丢弃插件实例前应导出。

## 验证结果

- `python --version`：3.12.16；`python -m pip check`：无依赖冲突。
- `python setup.py build_ext --inplace`：官方五个 C++ 扩展构建成功。
- 移植和新增兼容性测试：**103 passed**。包括已知真值直/弯结构、2D/3D、
  刚体变换、各向异性精度、连续投影端点、噪声、半径/面积/体积/Nyquist、
  质量控制、独立模块失败、全部 CSV、三表操作、发现/加载/帮助及快照线程。
  自动 XYZ 噪声排除另验证仅执行一次、忽略手动阈值，并向不同维度/bins/smoothing
  的两个模块传递相同的保留点集。
- 官方相关回归：**173 passed, 2 skipped, 1 failed**。
  失败为 `test_render_tab_controls.py::test_the_quick_buttons_wrap_only_when_the_row_is_full`，
  Windows offscreen 下按钮未按预期换行；在 `build/official-baseline` 中提取未改动的
  目标 HEAD 源码、使用相同解释器和编译扩展运行该测试，同样失败。
  本次没有修改其布局实现或放宽其断言。两项 slow 测试默认跳过。
- 已通过实际扫描入口建立 `PluginPanel`，在 Qt worker 分析已知模拟 tube，
  Density 得到 `exploratory`；生成并检查面板、中心线、径向分布、分段图预览。
- 未宣称覆盖整个官方测试套件、GPU、长耗时性能测试或用户真实实验数据。

复现命令（在仓库根目录）：

```powershell
.\.venv\python.exe -m pytest tests/centerline_extension -q --basetemp=build/centerline-check
.\.venv\python.exe -m pytest tests/test_plugins.py tests/test_user_plugins.py tests/test_plugin_grouping.py tests/test_fitter.py tests/test_fit_plugin.py tests/test_drift.py tests/test_locs_io.py tests/test_render.py tests/test_gui_layers.py tests/test_gui_render_axes.py tests/test_render_tab_controls.py tests/test_gui_preview.py tests/test_session_restore.py -q --basetemp=build/centerline-regression-check
```

`--basetemp` 指向专用测试临时目录，pytest 会清理该目录，勿改成数据目录。
旧测试适配仅改变导入路径、已移除的 Layer 状态 API 和新版结果按钮名，
未降低数值断言；失效分组测试改用新版 `set_group_settings` 触发真实 cache 失效。

## 后续官方更新的维护点

1. 保留/适配上述四个官方文件中的小型接入点；官方若提供等效接口即可替换。
2. 数据侧核对 `Session.scratch/table/group_settings/selects_slab`、
   `Layer.ranges/state.settings.axes`、`Context/Selection/Localizations`。
   分组参数和 grouping 算法属于官方版本行为，迁移数值算法本身不重写它们。
3. GUI 侧核对 `create_result_widget`、`result_updated`、`ended`，以及隔离适配器用到的
   `SettingsForm.fields/row_layouts/set_active/set_shown`、`PluginPanel.busy/_thread`、
   `result/result_button/plot/_window/status`。
4. 插件扫描、旁置 Markdown 帮助、`Result.figures()` 和 `Plot.draw_into()`。
5. 更新后先运行扩展测试，再运行相关官方回归；保留源摘要检查，算法变化单独评审。

## 默认启动发现修复

此前专用启动器之外的正常启动没有注册额外目录，导致菜单缺少插件。现已在 roots() 中接入源码仓库目录，并增加从任意工作目录、空用户目录和无额外配置时仍可发现/加载插件及帮助的回归测试。

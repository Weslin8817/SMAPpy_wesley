# Radial distance and density / Labeling Density

仓库内可追踪的用户插件入口。算法、累计结果、CSV 和 GUI 适配实现位于
`src/smappy_centerline/`；这里通过新版官方的额外插件目录机制加载。

本源码仓库的普通启动现在会自动扫描 `user_plugins/`，无需额外注册。
关闭旧进程后，可直接使用项目 `.venv` 启动：

```powershell
.\.venv\python.exe -m smappy.cli.gui
```

也可使用独立配置启动器：

```powershell
.\.venv\python.exe scripts/run_centerline.py
```

也可在命令后传入 localization 文件。此启动器把配置保存在本仓库的
`build/centerline-config/`，自动注册 `user_plugins/`，不覆盖正常启动的个人配置。
该目录只存放运行配置，不是自定义代码的唯一副本。

本源码仓库无需修改个人配置。对于其他安装位置，仍可在 **Preferences → Plugin folders → Add folder...**
添加本仓库的 `user_plugins` 绝对路径，然后使用原来的启动方式。

入口：**Plugins → Analysis → Measure → Radial distance and density**。
在 Render 选中 localization layer 并画一个 ROI，设置两个模块及共享噪声筛选，
点击 **Analyze/Add current ROI**。也可通过标签页的 `+` 添加插件。

每次分析一个 ROI，结果在当前插件实例内累计。三表联动选择，支持多选删除、
清空、查看图和 CSV 导出；关闭并丢弃插件实例后结果不自动持久化，需先导出。

密度和 Nyquist 是 **grouped-event** 口径，不是已校准的 independent labeling site
密度；缺少独立位点校准时对应字段保持不可用。

详细功能说明见 [插件帮助](Analysis/Measure/centerline_radial.md)，
移植架构、文件清单、测试及维护说明见 [迁移记录](../docs/centerline-extension.md)。

"""Launch this checkout with its repository-tracked measurement plugin.

Run with .venv/python.exe scripts/run_centerline.py [FILE.hdf5].
Configuration stays in build/centerline-config, not the normal user profile.
"""
from pathlib import Path
import os


def configure():
    root = Path(__file__).resolve().parents[1]
    os.environ["SMAPPY_CONFIG_DIR"] = str(root / "build" / "centerline-config")
    from smappy import config
    config.load(reload=True)
    plugin_root = root / "user_plugins"
    roots = config.plugin_roots()
    if plugin_root not in roots:
        config.set_plugin_roots([*roots, plugin_root])
    return plugin_root


if __name__ == "__main__":
    configure()
    from smappy.cli.gui import main
    main()

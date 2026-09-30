`pip install eon` declares PyYAML, which `eon.config` imports at start-up.
Without it `python -m eon` failed with `ModuleNotFoundError: No module named
'yaml'` in any environment other than the pixi and conda ones, which already
listed it.

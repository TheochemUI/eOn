The basic-eon workflow runs the tests/ suite on Unix after the
client install, including the `sh` package those process tests import.
Job modules import `version` from the `eon` package so a source-tree
checkout without generated `version.py` still dispatches.

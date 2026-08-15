import importlib as _importlib
from eon.server import server

try:
    from eon.version import version as __version__
except ModuleNotFoundError:
    # version.py is generated (tools/gitversion.py, see eon/meson.build) and
    # gitignored, so a checkout that has not been built through meson does not
    # carry one. Fall back to installed metadata, then to a placeholder, so
    # importing eon from a source tree works.
    from importlib.metadata import PackageNotFoundError, version as _pkg_version

    try:
        __version__ = _pkg_version("eon")
    except PackageNotFoundError:
        __version__ = "0.0.0+unknown"

# Job modules import `from eon import version`. That name is the same
# string as __version__, whether it came from generated version.py or
# the fallback above.
version = __version__

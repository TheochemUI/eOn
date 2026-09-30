`eon.atoms` loads without minimage again. The `Cell.wrap_many` binding for
readcon-ops ran at import time, so the aKMC server stopped with
`ModuleNotFoundError: No module named 'minimage'` wherever minimage, which
is not on PyPI or conda-forge, was not installed by hand. The binding now
runs just before the readcon-ops calls that need it, and is skipped when
minimage is absent.

pyeonclient 0.4.0: the potentials it exposes come from ``librgpot`` now
rather than eOn's in-tree kernels, so the wheels bundle that library and
no longer carry per-pot plugin objects. The Python API is unchanged --
the same ``potential`` names select the same physics.

Its wheel workflow gained a ``target`` input so an upload can be
rehearsed against TestPyPI before the irreversible one; a
``pyeonclient-v*`` tag still goes straight to PyPI.

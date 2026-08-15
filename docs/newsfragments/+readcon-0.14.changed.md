C++ wrap and Python pins move to readcon 0.14.5 together, so the client
can read the spec-3 files the 0.14 writer emits. x-only constraints
round-trip on the client and the Python Structure writer (readcon-core #25).
0.14.5 ships a win_amd64 wheel, so Windows CI does not build the sdist.

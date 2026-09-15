Shared plugin TUs include Parameters.h so accessors compile when
Potential.h only forward-declares Parameters. Potential constructors
that take Parameters live in eoncbase so plugin .so files do not need
eonclib. Duplicate leftover field-style option blocks are gone.
readcon is resolved before plugin subdirs so SocketNWChem links rkr_*.

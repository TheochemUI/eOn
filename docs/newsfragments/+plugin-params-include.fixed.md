Shared plugin TUs include Parameters.h so accessors compile when
Potential.h only forward-declared Parameters. Potential constructors
that take Parameters are inline so plugin .so files do not need
eonclib. Duplicate leftover field-style option blocks are gone.
readcon is resolved before plugin subdirs so SocketNWChem links rkr_*.

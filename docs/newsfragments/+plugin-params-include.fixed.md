Shared plugin TUs include Parameters.h so accessors compile when
Potential.h only forward-declares Parameters. Duplicate leftover
field-style option blocks are gone from INI/JSON/bindings.
readcon is resolved before plugin subdirs so SocketNWChem sees
readcon-core.hpp and links rkr_*.

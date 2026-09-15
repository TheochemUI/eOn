Shared plugin TUs include Parameters.h so rgpot/xtb/metatomic
accessors compile when Potential.h only forward-declares Parameters.
LammpsLoader drops a duplicate ensure_loaded declaration.
readcon is resolved before plugin subdirs so SocketNWChem sees
readcon-core.hpp and links rkr_*.

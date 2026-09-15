Shared plugin TUs include Parameters.h so rgpot/xtb/metatomic
accessors compile when Potential.h only forward-declares Parameters.
LammpsLoader drops a duplicate ensure_loaded declaration.

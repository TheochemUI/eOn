Parameters option groups now expose only const accessors. INI/JSON loaders and
bindings write through ``ParametersLoadAccess``; MPI comm/rank use dedicated
setters so ``ParametersMpi.h`` no longer takes the address of a temporary.

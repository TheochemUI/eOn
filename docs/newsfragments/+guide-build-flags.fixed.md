The guide records that `-Dwith_xtb`, `-Dwith_water`, `-Dwith_metatomic`, and `-Dwith_vasp` default false, and that `-Dwith_lammps` and `-Dwith_dftd3` are not eOn options.
MPI runs launch `eonclient` and read `EON_CLIENT_STANDALONE` plus `EON_NUMBER_OF_CLIENTS`; the version record is 3.5.0, and the metatomic extra depends on `rgpot>=3.2.0`.

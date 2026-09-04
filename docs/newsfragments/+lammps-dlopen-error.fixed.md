``LammpsLoader::require_loaded`` now says whether ``liblammps.so`` is
missing, failed ``dlopen`` (glibc), or opened without
``lammps_open_no_mpi`` (the eOn plugin is not LAMMPS).

`-Dwith_mpi=enabled` requires an embedded Python (`python3-embed`): the server rank of the MPI client runs `EON_SERVER_PATH` through `Py_Main`. Configure stops when the embed dependency is not found.

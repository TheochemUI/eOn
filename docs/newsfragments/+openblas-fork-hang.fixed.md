Importing eon parks the OpenBLAS pool at one thread for the duration of a fork and restores the previous count in the parent, so a later factorization does not wait on server_lock.

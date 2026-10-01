The MPI communicator reads the job path a client returns with `ndarray.tobytes`. `tostring` is gone in NumPy 2, so the server stopped with an `AttributeError` at the first returned job.

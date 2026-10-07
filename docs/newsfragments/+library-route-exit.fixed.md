The library-route exit is registered again after the engine shuts down, so process exit does not run library destructors on a live MPI world.

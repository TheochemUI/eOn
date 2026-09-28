The local communicator starts each client in its own session and kills the whole process group on exit, so ExtPot wrappers and the MPI launchers they start do not outlive the server.

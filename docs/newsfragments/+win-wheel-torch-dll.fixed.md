The Windows metatomic wheel import loads the metatensor DLLs that match the installed torch. An older copy was being opened first, so the import failed on a missing procedure.

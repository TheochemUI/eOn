`PotRegistry` is a process-lifetime heap singleton, so a `Potential` destroyed during interpreter finalization no longer calls into a destroyed registry or locks a destroyed mutex.

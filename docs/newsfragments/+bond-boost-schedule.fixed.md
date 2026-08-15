Bond-boost hyperdynamics advances the equilibration counter once per MD step.
`boost()` only evaluates the current bias, so ParallelReplica (which also
installs the bias potential on the trajectory) and SafeHyper share the same
`rmd_time` schedule.

NEB `minimize_endpoints` defaults to false, matching SVN and the
pyeonclient NEB API. The old true default added endpoint relaxations
to `total_force_calls` on fixtures that never asked for them.

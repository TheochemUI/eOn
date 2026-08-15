Copying a Matter leaves biasPotential null. Copy-assign used to skip
that pointer, so getBiasForces read an indeterminate BondBoost and
the Python suite died on main.

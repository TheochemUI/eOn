Dynamics saddle search with bias_potential set to bond_boost now calls setBiasPotential and BondBoost::advance once per MD step, so the trajectory includes the bond-boost force.

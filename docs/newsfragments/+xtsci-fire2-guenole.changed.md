Version 2.0 of the fast inertial relaxation engine (FIRE) now runs for `[Xtsci] method = fire2`.
It runs algorithm 2 of Guenole et al. (2020) with the parameters of Table 1 in that paper, and the client selects it through `rgmin_solver_set_fire_variant`.
Before, `fire2` ran the default rgmin `RGMIN_FIRE2` schedule with the FIRE 1.0 parameters and no half-step correction or initial delay.
`Optimizer::run` keeps `fire2` on the rgmin session, since `rgmin_minimize_eindir` takes no variant.
Against rgmin older than application binary interface (ABI) minor 27, the client refuses `fire2`.
`method = fire` and the native FIRE are unchanged.

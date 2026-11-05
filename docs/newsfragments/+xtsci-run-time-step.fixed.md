`Optimizer::run` with `opt_method = xtsci` now passes `time_step` as the time step of `fire`, `fire2` and `quickmin`, as `step()` already did.
With eindir-core linked, the call to `rgmin_minimize_eindir` used `max_move` instead, so `time_step` had no effect there.

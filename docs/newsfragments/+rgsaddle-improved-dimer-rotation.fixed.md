The improved dimer, `[Dimer] improved = true`, is the default.
With rgsaddle it now sends `rotation_tol = 0` to the min-mode session in place of `[Dimer] torque_min`.
Its rotation stops on `converged_angle` or `rotations_max`, as in eOn's own improved dimer.
The old value let rgsaddle stop at a rotational force of `torque_min`.
From rgsaddle 0.3.0 it also turned on the classic torque window, which can stop after one rotation.
The classic dimer and Lanczos keep `torque_min`.

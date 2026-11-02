`[Saddle Search]` gains `irc_after`, `irc_step_size`, `irc_converged_force` and `irc_max_steps`.
With `irc_after = true`, rgsaddle follows the intrinsic reaction coordinate (IRC) from a converged saddle.
It goes downhill in both directions when the saddle has a negative eigenvalue.
eOn's minimizer then relaxes both ends.
The job writes `irc.con`, `irc.dat`, `irc_min1.con` and `irc_min2.con`, plus seven `results.dat` entries.
The step defaults to 0.1 amu^1/2 Angstroms, where amu is the atomic mass unit.
The force threshold defaults to 0.05 electronvolts/Angstrom (eV/Angstrom).
A branch ends at a minimum when rgsaddle reports one.
rgsaddle does so only after the first step and where the largest atomic force is at or below the threshold.
The option is off by default, and with it off a saddle search writes the same files as before.
It needs eOn compiled against rgsaddle C application binary interface (ABI) minor 14 or newer.
Other builds refuse it before the search.
On the systems benchmarked the IRC reached the same minima as two minimizations displaced along the mode.
It always cost more force calls.
The dense model inside the rgsaddle IRC makes it slow for many free atoms.

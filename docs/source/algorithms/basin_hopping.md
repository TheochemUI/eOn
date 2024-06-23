# Basin Hopping

Basin hopping is a Monte Carlo method in which the energy of each configuration
is taken to be the energy of a local minimum. {cite}`walesGlobalOptimizationBasinHopping1997`

At each basin hopping step the client will print out the current energy
(current), the trial energy (trial), the lowest energy found (global min), the
number of force calls needed to minimize the structure (fc), the acceptance
ratio (ar), and the current max displacement (md).

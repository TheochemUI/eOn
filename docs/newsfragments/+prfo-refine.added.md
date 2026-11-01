`[Saddle Search]` gains `prfo_refine`, `prfo_refine_force` and `prfo_refine_max_iterations`.
They are opt-in and off by default, and every job without them writes the same outputs as before.
With `prfo_refine = true`, a converged dimer or Lanczos saddle and its lowest mode go to the rgsaddle Sella session, which refines the saddle with partitioned rational function optimization (P-RFO) to `prfo_refine_force` in max-atom force.
A climbing image nudged elastic band (NEB) job refines a copy of its climbing image along the band tangent.
It writes the result to `sp_refined.con` and leaves the band unchanged.
Measured against rgsaddle 0.3.2 with its proposed P-RFO eigensolve series, `benchmarks/data/one_pt_saddle_search` (one free atom) reached 1.1e-7 eV/Å in 43 force calls, against 54 calls for a dimer search with a 1e-5 eV/Å tolerance.
On the Pt heptamer and the 13-atom Lennard-Jones (LJ13) cluster, the dimer converged to 1e-5 eV/Å needed fewer calls than the refinement.
On `benchmarks/data/neb_morse_pt` (175 free atoms) the refinement converged in 465 calls and 1207 s, against 427 calls and 3.2 s for the band converged to 1e-5 eV/Å.
On `examples/neb-al` (401 free atoms) it did not finish within 3600 s, because the TS-BFGS update (a Broyden-Fletcher-Goldfarb-Shanno variant for transition states) diagonalizes the whole Hessian model after every force call.
At the pinned rgsaddle 0.2.2 every dense eigensolve is a Jacobi solve, and neither band refinement is practical.
The keys need eOn built with `-Dwith_rgsaddle=true` against rgsaddle application binary interface (ABI) 1.15 or newer.

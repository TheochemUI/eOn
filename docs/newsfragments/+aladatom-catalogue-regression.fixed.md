The built-in `eam_al` adatom AKMC templates (`eon/aladatom-akmc` at 250 K and
`eon/standard-aladatom-akmc` at 150 K) no longer lose catalogue diversity from
the first min-mode L-BFGS translation step. Those rows use the classic dimer
(`improved = false`). A mild negative mode again keeps the capped
mode-curvature scale (`H0 = 1/|C_tau|`) instead of the finite-difference
effective-force probe. Improved-dimer and Lanczos translations still take
that probe, including when the mode is mild, which is the scaling the LJ13
and Pt heptamer force-call counts follow. A steep classic-dimer overrun still
uses the probe.

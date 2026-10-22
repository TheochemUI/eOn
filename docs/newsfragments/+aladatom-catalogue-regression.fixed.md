The built-in `eam_al` adatom AKMC templates (`eon/aladatom-akmc` at 250 K and
`eon/standard-aladatom-akmc` at 150 K) no longer lose catalogue diversity from
the first min-mode L-BFGS translation step. For mild negative modes the search
again keeps the capped mode-curvature scale (`H0 = 1/|C_tau|`) instead of
forcing the finite-difference effective-force probe reserved for steep-mode
overruns. This restores the pre-regression adatom rate and catalogue accuracy.

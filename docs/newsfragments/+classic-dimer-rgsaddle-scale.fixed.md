A classic dimer (`improved = false`) translated with L-BFGS keeps the capped
mode-curvature scale when the rotation runs in rgsaddle. That build uses one
solver for every min-mode method, so the scale follows the dimer setting
rather than the solver type. Improved dimer and Lanczos stay on the
effective-force probe, which is the scaling the LJ13 and Pt heptamer
force-call counts follow.

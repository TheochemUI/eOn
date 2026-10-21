Classical dimer rotation now preserves the signed `atan` limit when the
angle-fit denominator is near zero, so near-singular rotations take the
finite-angle directional update instead of collapsing to a zero-angle
fallback. This reduces stalled min-mode rotations in bowl-breakout saddle
searches.

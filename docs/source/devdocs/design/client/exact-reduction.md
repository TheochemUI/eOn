---
myst:
  html_meta:
    "description": "Why eOn computes the dot products and norms behind optimizer decisions with correctly rounded reductions, and what the exact path guarantees."
    "keywords": "eOn reproducibility, exact dot product, correctly rounded, Ozaki scheme, Eigen reduction order, bitwise reproducible"
---

# Exact reductions

`include/eon/ExactReduction.h` provides two kinds of reduction.

- `eonc::exact::dot`, `sum`, `sumStrided`, `squaredNorm`, `norm` and
  `normalized` return the exact value rounded once to the nearest double,
  ties to even. The optimizers (L-BFGS, CG, FIRE, Quickmin, steepest descent)
  and the `norm` and `rms` convergence tests use them for every scalar they
  branch on. `sumStrided` is `sum` for a column of a row-major atom matrix.
- `eonc::exact::orderedDot` and `orderedNorm` add rounded products in a tree
  fixed in the source: four lanes, then `(l0 + l1) + (l2 + l3)`, compiled
  without contraction. They are not exact, but their bits do not depend on
  the compiler, the Eigen release or the CPU. The L-BFGS two-loop recursion
  (which sets the step vector) and the per-atom lengths of the step caps and
  push-apart use them.

## Why

Floating-point addition does not associate, so a sum depends on the order of
its terms. Eigen chooses that order, and the order changes with the release
and the SIMD width: for `dot` and `squaredNorm`, Eigen 3.4 adds two packet
accumulators and Eigen 5.0 four (`Eigen/src/Core/InnerProduct.h`). With 39
coordinates (LJ13) under SSE2 the two releases build different summation
trees.

Optimizers amplify a one-ulp difference. The basin-hopping LJ13 case in
`test_jobs` first diverged at the finite-difference curvature of the first
L-BFGS step (`C` in `LBFGS::getStep` differed by one ulp). The run then ended
at 1640 force calls with Eigen 3.4 and 1700 with Eigen 5.0, under gcc 14 and
gcc 15 alike. A plain dot product of length `n` carries an error of up to
`n u sum |x_i y_i|` (`u = 2^-53`), so two orders can differ by twice that;
near a convergence threshold or an L-BFGS curvature test, that difference
flips a decision.

The mean force removed when `removeNetForce` is set is that same exact sum
of each column, divided by the atom count. Eigen's `colwise().sum()` follows
the packet width, so the mean differed between x86 and NEON and the LJ13
basin hop took 1671 force calls on Apple clang instead of 1640. An exact
mean left 1658. Rewriting the L-BFGS curvature probe and the two-loop step
as a multiply followed by an add left the count at 1658, so that step was
not the remaining gap.

The gap is the LJ pair evaluation. On x86 without FMA in the compile ISA,
rgpot compiles the pair loop and the inlined kernel in `accumulateHw`:
contraction off, and an explicit fma only for the force update. ARM has FMA
in the base ISA, so that clone is not built. Apple clang's default
contraction then fuses `dx*dx + dy*dy + dz*dz` and the 12-6 products
(`b*(a-1) - shift`, `2*a - 1`). The LJ objects are compiled with contraction
off, the same flag the minimum-image objects already use, and the force
update stays `__builtin_fma`. That is the rounding of the Linux clone.

## How

Every finite double is an integer significand `m < 2^53` times `2^e`.
`x_i y_i` is therefore the 106-bit integer `m_x m_y` times `2^(e_x + e_y)`,
with no rounding. The products go into a fixed-point integer accumulator
that spans the whole exponent range, 135 words of 32 bits in signed 64-bit
storage. Integer addition associates, so the accumulator contains the same
exact sum for any order. One final rounding gives `RN(sum x_i y_i)`.

This is the integer form of the error-free splitting of the Ozaki scheme
([Ozaki, Ogita, Oishi and Rump 2012](https://doi.org/10.1007/s11075-011-9478-1);
[Ozaki, Uchino and Imamura 2025](https://doi.org/10.48550/arXiv.2504.08009)).
There, each operand splits into slices whose products are exact and whose
sums are exact in any order. Compensated algorithms such as `Dot2`
([Ogita, Rump and Oishi 2005](https://doi.org/10.1137/030601818)) improve
accuracy but stay order-dependent. Reproducible summation in the sense of
[Demmel and Nguyen 2015](https://doi.org/10.1109/TC.2014.2345391) is
order-independent without being exact. A correctly rounded result is both.

`dot` first tries a faster route to the same number. It runs a compensated
dot product (`Dot2`) in four lanes, with TwoProduct by FMA where the CPU has
it and Veltkamp-Dekker splitting where it does not; both give the same exact
error term. It also bounds the error of the compensation sum a posteriori.
When the exact value then lies strictly between the two midpoints around
`r = RN(s + c)`, `r` is the correctly rounded result. Near a tie, or when a
range guard fails, the integer accumulator decides. Both routes return
`RN(sum x_i y_i)`, so the route never shows in the bits.

Non-finite inputs fall back to the plain IEEE sum, so NaN and infinity
propagate unchanged. An exact zero returns `+0`.

## Cost

Measured on one Xeon Platinum 8358 node (Elja, gcc 15.3, Eigen 5.0.1), with
run directories on local disk, one warm-up run, then five interleaved runs;
the table gives medians.

| workload                           | before (s) | after (s) | change |
|------------------------------------|-----------:|----------:|-------:|
| LJ13 basin hopping, L-BFGS, 200 steps |      0.221 |     0.241 |   +9 % |
| basin-hopping example, CG          |      2.58  |     2.18  | within noise |
| Al NEB, 7 images, L-BFGS           |      0.957 |     0.998 |   +4 % |
| Pt saddle search, dimer            |      0.267 |     0.269 |   +1 % |

LJ13 has the cheapest force call of the four, so it shows the largest
relative cost. The trajectories differ from the old ones (14469 against
14348 force calls for LJ13), so part of each change is a different path,
not overhead. In isolation, `dot` costs 3.4 ns per term at n = 39 and 1.0 ns
per term at n = 1008, against 0.6 and 0.2 ns for Eigen; `orderedDot` matches
Eigen at n = 39 and costs 1.9 times as much at n = 1008.

## Checks

- `client/validation/exact_reduction_cases.py` writes reference cases from
  exact rational arithmetic: cancellation, ties, sticky bits, subnormal
  results, the overflow threshold, and optimizer-shaped vectors.
  `test_exact_reduction` checks each case bitwise and under 16 permutations.
- `client/validation/exact_reduction_derivations.py` (sympy) derives the Ozaki
  slice condition `beta >= (p + log2 n) / 2`, the portable 128-bit product
  and the accumulator piece split. It also lists the Eigen 3.4 and 5.0 trees
  for `n = 39`.
- `client/validation/exact_reduction_bounds.sollya` bounds the plain-sum
  error, the `norm` error (`1.5 u`) and the accumulator range.
- `formal/` (Lean 4 and Mathlib, `lake build`) checks four contracts. First,
  floating-point sums on a common grid are exact in every parenthesisation
  within the exponent range, so any two orders agree. Second, the Ozaki
  slice-product condition. Third, the final rounding step is round to
  nearest, ties to even. Fourth, accumulator words cannot overflow between
  carry normalisations.

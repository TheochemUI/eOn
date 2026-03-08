---
myst:
  html_meta:
    "description": "Reference for SafeMath.h guarded arithmetic utilities in eOn."
    "keywords": "eOn SafeMath, floating-point exceptions, SIGFPE, guarded division, numerical safety"
---

# SafeMath utilities

`client/SafeMath.h` provides lightweight guarded arithmetic to prevent
floating-point exceptions (SIGFPE) in numerical code. All functions live in
`eonc::safemath`.

## Scalar functions

All are `inline constexpr` or `inline` (for `<cmath>` calls).

| Function | Signature | Behavior |
|----------|-----------|----------|
| `safe_div` | `(double num, double denom, double fallback=0.0)` | Returns `fallback` when `abs(denom) < eps` |
| `safe_recip` | `(double x, double fallback=0.0)` | `safe_div(1.0, x, fallback)` |
| `safe_acos` | `(double x)` | `std::acos(std::clamp(x, -1.0, 1.0))` |
| `safe_sqrt` | `(double x)` | `std::sqrt(std::max(0.0, x))` |
| `safe_atan_ratio` | `(double num, double denom, double fallback=0.0)` | Guarded `std::atan(num/denom)` |

The epsilon (`1e-300`) is well below any physically meaningful value but above
the FPE trap threshold.

## Eigen template

Available only when `Eigen/Core` is already included (guarded by
`#ifdef EIGEN_CORE_H`):

| Function | Signature | Behavior |
|----------|-----------|----------|
| `safe_normalized` | `(Eigen::MatrixBase<D> const& v, double min_norm=eps)` | Returns zero matrix/vector when `v.norm() < min_norm` |

## Design principles

- **Fallback values match existing control flow.** Each call site chooses a
  fallback that triggers the same branch the code would take for degenerate
  input (reset, skip, clamp, etc.). Valid inputs produce bit-identical results.

- **No overhead on the hot path.** The guard is a single `abs(x) < eps`
  comparison, which the branch predictor will learn is almost never taken.

- **Header-only, no link dependencies.** Just `#include "SafeMath.h"` after
  your Eigen includes (if using `safe_normalized`).

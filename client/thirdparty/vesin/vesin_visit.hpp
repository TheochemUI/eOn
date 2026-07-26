// Header-only fused MIC pair visitation for vesin's CPU brute-force
// algorithm (eOn extension, upstream candidate).
//
// The C API's vesin_neighbors_visit pays an indirect call per in-cutoff
// pair, which forces register spills inside the scan loop. This template
// lets the compiler inline the visitor into the single O(n^2) scan, and
// collects the candidate pair list (as compact int32 index pairs) in the
// same pass. Orthorhombic / free-boundary boxes only; exactly one pair per
// (i, j) — the nearest periodic image, folded with the widths in `w`.
#ifndef VESIN_VISIT_HPP
#define VESIN_VISIT_HPP

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace vesin {
namespace cpu {

/// Round to nearest via truncate-cast: no libm floor call on -march
/// targets without roundsd. Half-integer inputs round away from zero.
inline double visit_round(double t) {
    return static_cast<double>(static_cast<long long>(t + std::copysign(0.5, t)));
}

/// Scan all i < j pairs of `points` (interleaved x,y,z, `n` atoms), folding
/// each dimension by width `w[k]` when `inv[k] = 1/w[k]` is nonzero
/// (`inv[k] == 0` disables the fold for that dimension).
///
/// Pairs with `r2 < list_cutoff2` append (i, j) to `pairs` (flat int32,
/// cleared first); pairs with `r2 <= visit_cutoff2` additionally invoke
/// `visit(i, j, dx, dy, dz, r2)` with the folded `r_j - r_i` vector. Pass a
/// negative `visit_cutoff2` for a collect-only scan.
template <typename Visitor>
inline void brute_force_visit(
    const double* points,
    std::size_t n,
    const double w[3],
    const double inv[3],
    double list_cutoff2,
    double visit_cutoff2,
    std::vector<int32_t>& pairs,
    Visitor&& visit
) {
    pairs.clear();
    if (pairs.capacity() == 0) {
        // First build: skip the push_back growth-realloc ladder. ~100
        // neighbours per atom covers dense metallic cutoffs; denser systems
        // grow once.
        pairs.reserve(n * 200);
    }
    const double w0 = w[0], w1 = w[1], w2 = w[2];
    const double i0 = inv[0], i1 = inv[1], i2 = inv[2];
    for (std::size_t i = 0; i < n; i++) {
        const double xi = points[3 * i];
        const double yi = points[3 * i + 1];
        const double zi = points[3 * i + 2];
        for (std::size_t j = i + 1; j < n; j++) {
            double dx = points[3 * j] - xi;
            double dy = points[3 * j + 1] - yi;
            double dz = points[3 * j + 2] - zi;
            dx -= w0 * visit_round(dx * i0);
            dy -= w1 * visit_round(dy * i1);
            dz -= w2 * visit_round(dz * i2);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < list_cutoff2) {
                pairs.push_back(static_cast<int32_t>(i));
                pairs.push_back(static_cast<int32_t>(j));
                if (r2 <= visit_cutoff2) {
                    visit(static_cast<int32_t>(i), static_cast<int32_t>(j),
                          dx, dy, dz, r2);
                }
            }
        }
    }
}

/// Evaluation-only variant: same scan and fold, no pair collection. Serves
/// first-sighting evaluations where list capture would be paid for nothing.
template <typename Visitor>
inline void brute_force_visit_only(
    const double* points,
    std::size_t n,
    const double w[3],
    const double inv[3],
    double visit_cutoff2,
    Visitor&& visit
) {
    const double w0 = w[0], w1 = w[1], w2 = w[2];
    const double i0 = inv[0], i1 = inv[1], i2 = inv[2];
    for (std::size_t i = 0; i < n; i++) {
        const double xi = points[3 * i];
        const double yi = points[3 * i + 1];
        const double zi = points[3 * i + 2];
        for (std::size_t j = i + 1; j < n; j++) {
            double dx = points[3 * j] - xi;
            double dy = points[3 * j + 1] - yi;
            double dz = points[3 * j + 2] - zi;
            dx -= w0 * visit_round(dx * i0);
            dy -= w1 * visit_round(dy * i1);
            dz -= w2 * visit_round(dz * i2);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 <= visit_cutoff2) {
                visit(static_cast<int32_t>(i), static_cast<int32_t>(j),
                      dx, dy, dz, r2);
            }
        }
    }
}

} // namespace cpu
} // namespace vesin

#endif

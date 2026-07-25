/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Shared neighbor lists via vesin (single NL backend for classical pots
** and Metatomic). Do not reimplement Verlet / link-cell / O(N^2) loops.
*/
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

// Vendored / system vesin C API
#include "vesin.h"

namespace eonc {

/// RAII wrapper around ``vesin_neighbors`` for Matter-style boxes
/// (``double[9]`` row-major 3×3 cell, same layout as ``Matter::cell``).
///
/// Pair vector convention matches vesin: ``r_ij = r_j - r_i + S @ H``.
/// Classical pair pots that used ``r_i - r_j`` must flip the sign.
class VesinNeighbors {
public:
  struct Options {
    double cutoff{0.0};
    bool full{false}; ///< true: both i→j and j→i; false: half list
    bool sorted{false};
    bool return_shifts{false};
    bool return_distances{true};
    bool return_vectors{true};
    std::array<bool, 3> periodic{{true, true, true}};
    /// vesin 0.6 native Verlet skin: > 0 lets vesin cache the topology in
    /// this list until an atom moves more than ``skin/2``. Single-geometry
    /// consumers only — multi-image callers go through PairListCache.
    double skin{0.0};
    /// Threads for the vesin build; 1 keeps small builds deterministic and
    /// spawn-free, 0 defers to OMP_NUM_THREADS / core count.
    int32_t n_threads{1};
    /// vesin pair-search algorithm (auto / brute-force MIC / cell list).
    VesinAlgorithm algorithm{VesinAutoAlgorithm};
  };

  VesinNeighbors() = default;
  VesinNeighbors(const VesinNeighbors &) = delete;
  VesinNeighbors &operator=(const VesinNeighbors &) = delete;
  VesinNeighbors(VesinNeighbors &&other) noexcept;
  VesinNeighbors &operator=(VesinNeighbors &&other) noexcept;
  ~VesinNeighbors();

  /// Build / rebuild the list. ``R`` is length ``3*n`` (x,y,z interleaved);
  /// ``box`` is length 9 (row-major cell vectors).
  ///
  /// Keep the same ``VesinNeighbors`` instance across force evaluations so
  /// vesin can re-use pair/distance/vector allocations (see ``vesin_neighbors``
  /// docs). Freeing before every ``compute`` is incorrect and expensive.
  void compute(const double *R, std::size_t n, const double *box,
               const Options &opt);

  [[nodiscard]] std::size_t size() const { return list_.length; }

  [[nodiscard]] std::size_t i(std::size_t p) const { return list_.pairs[p][0]; }
  [[nodiscard]] std::size_t j(std::size_t p) const { return list_.pairs[p][1]; }

  /// Distance |r_ij| (requires ``return_distances``).
  [[nodiscard]] double distance(std::size_t p) const {
    return list_.distances[p];
  }

  /// Vector r_j − r_i (+ PBC). Requires ``return_vectors``.
  [[nodiscard]] const double *vector(std::size_t p) const {
    return list_.vectors[p];
  }

  [[nodiscard]] const int32_t *shift(std::size_t p) const {
    return list_.shifts[p];
  }

  /// Raw list for advanced consumers (Metatomic tensor conversion).
  [[nodiscard]] const VesinNeighborList &raw() const { return list_; }
  [[nodiscard]] VesinNeighborList &raw() { return list_; }

  /// Steal ownership of the internal list (for Metatomic custom deleters).
  /// After release, this object is empty until the next ``compute``.
  VesinNeighborList *release();

private:
  VesinNeighborList list_{};
  bool owns_{false};

  void free_list();
};

/// Verlet-skin cached half pair list on top of vesin.
///
/// vesin builds the pair set at ``cutoff + skin`` (with cell shifts); between
/// rebuilds the cached pairs are re-evaluated against the *current* positions,
/// so distances and vectors are always exact. The standard Verlet guarantee
/// applies: while every atom has moved less than ``skin/2`` since the build,
/// the cached set is a superset of all pairs within ``cutoff``; the per-pair
/// ``r2 <= cutoff^2`` filter in ``forEach`` keeps the physics identical to a
/// fresh build.
///
/// The pairs live in the slot's own vesin buffers — no unsorted->sorted
/// permutation and no copy on rebuild, which is what a one-shot evaluation
/// (point job, first NEB iteration) actually pays for.
class CachedPairList {
public:
  struct Options {
    double cutoff{0.0};
    double skin{1.0};
    std::array<bool, 3> periodic{{true, true, true}};

    bool operator==(const Options &o) const {
      return cutoff == o.cutoff && skin == o.skin && periodic == o.periodic;
    }
  };

  /// True while the cached pair set is valid for positions ``R``: same atom
  /// count, options and box as the build, and max displacement < skin/2.
  [[nodiscard]] bool valid(const double *R, std::size_t n, const double *box,
                           const Options &opt) const;

  /// Rebuild via vesin at ``cutoff + skin`` into this slot's own buffers.
  void rebuild(const double *R, std::size_t n, const double *box,
               const Options &opt);

  /// Visit every cached pair within the true cutoff of the current
  /// positions. ``fn(i, j, dx, dy, dz, r2)`` receives the *historical* eOn
  /// convention ``d = r_i - r_j`` (minimum image / shift applied), i.e. the
  /// negated vesin vector.
  ///
  /// Two evaluation modes, chosen at rebuild:
  /// - MIC (orthorhombic box, true cutoff <= half the smallest periodic
  ///   width): vesin brute-force enumerates candidates and the fold happens
  ///   here per call, so an image swap between rebuilds resolves exactly as
  ///   the historical per-call MIC loops did.
  /// - shift (everything else): cell-list pairs carry their build-time cell
  ///   shift; every image is its own pair entry.
  template <typename Fn> void forEach(const double *R, Fn &&fn) const {
    const double cutoff2 = opt_.cutoff * opt_.cutoff;
    const VesinNeighborList &vl = nl_.raw();
    if (mic_) {
      if (micInv_[0] == 0.0 && micInv_[1] == 0.0 && micInv_[2] == 0.0) {
        // Free boundaries: no folds at all.
        for (std::size_t p = 0; p < vl.length; ++p) {
          const auto i = static_cast<int32_t>(vl.pairs[p][0]);
          const auto j = static_cast<int32_t>(vl.pairs[p][1]);
          const double dx = R[3 * i] - R[3 * j];
          const double dy = R[3 * i + 1] - R[3 * j + 1];
          const double dz = R[3 * i + 2] - R[3 * j + 2];
          const double r2 = dx * dx + dy * dy + dz * dz;
          if (r2 <= cutoff2) {
            fn(i, j, dx, dy, dz, r2);
          }
        }
        return;
      }
      const double w0 = boxref_[0], w1 = boxref_[4], w2 = boxref_[8];
      const double i0 = micInv_[0], i1 = micInv_[1], i2 = micInv_[2];
      for (std::size_t p = 0; p < vl.length; ++p) {
        const auto i = static_cast<int32_t>(vl.pairs[p][0]);
        const auto j = static_cast<int32_t>(vl.pairs[p][1]);
        double dx = R[3 * i] - R[3 * j];
        double dy = R[3 * i + 1] - R[3 * j + 1];
        double dz = R[3 * i + 2] - R[3 * j + 2];
        // Non-periodic dims store inv = 0, so the fold rounds 0 to 0.
        dx -= w0 * micRound(dx * i0);
        dy -= w1 * micRound(dy * i1);
        dz -= w2 * micRound(dz * i2);
        const double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 <= cutoff2) {
          fn(i, j, dx, dy, dz, r2);
        }
      }
      return;
    }
    const double *H = boxref_.data();
    for (std::size_t p = 0; p < vl.length; ++p) {
      const auto i = static_cast<int32_t>(vl.pairs[p][0]);
      const auto j = static_cast<int32_t>(vl.pairs[p][1]);
      double dx = R[3 * i] - R[3 * j];
      double dy = R[3 * i + 1] - R[3 * j + 1];
      double dz = R[3 * i + 2] - R[3 * j + 2];
      const int32_t sa = vl.shifts[p][0];
      const int32_t sb = vl.shifts[p][1];
      const int32_t sc = vl.shifts[p][2];
      if ((sa | sb | sc) != 0) { // wrapped pair: d -= S @ H
        dx -= sa * H[0] + sb * H[3] + sc * H[6];
        dy -= sa * H[1] + sb * H[4] + sc * H[7];
        dz -= sa * H[2] + sb * H[5] + sc * H[8];
      } else if (i == j) {
        continue; // degenerate zero-shift self pair
      }
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 <= cutoff2) {
        fn(i, j, dx, dy, dz, r2);
      }
    }
  }

  [[nodiscard]] std::size_t pairCount() const { return nl_.size(); }

private:
  /// Round to nearest via SSE2 truncate-cast: std::floor(t + 0.5) lowers to
  /// a libm call under -march targets without roundsd and costs ~20 cycles
  /// per fold; this is ~5. Half-integer inputs (distance exactly half the
  /// box) round away from zero instead of up — a measure-zero difference.
  static double micRound(double t) {
    return static_cast<double>(
        static_cast<long long>(t + std::copysign(0.5, t)));
  }

  VesinNeighbors nl_;
  std::vector<double> Rref_;
  std::array<double, 9> boxref_{};
  std::array<double, 3> micInv_{}; ///< 1/width per periodic dim, else 0
  Options opt_{};
  std::size_t n_{0};
  bool mic_{false};
  /// MIC-mode candidate set that is the complete pair graph: valid for any
  /// positions with the same atoms/box/options (the eval fold re-derives
  /// exact geometry, and no pair can ever leave a complete candidate set).
  bool complete_{false};
  bool built_{false};
};

/// Process-wide pool of CachedPairList slots, matched by geometry proximity.
/// NEB evaluates several images through one shared Potential instance — and
/// spawns fresh threads every iteration (NudgedElasticBand::updateForces), as
/// does ImprovedDimer for one endpoint — so neither a per-instance cache
/// (data race on the shared pot) nor thread_local storage (dies with each
/// iteration's threads) survives. Each ``ensure`` call finds the slot whose
/// reference positions are within skin/2 of ``R`` (MRU order) or builds a new
/// slot. Slots are immutable after build and handed out as shared_ptr, so
/// readers never race eviction; the pool mutex guards only the match/insert,
/// never the force loops.
class PairListCache {
public:
  /// Slot count covers NEB default image counts; more images degrade to
  /// occasional rebuilds rather than incorrect results.
  static constexpr std::size_t kMaxSlots = 8;

  std::shared_ptr<const CachedPairList>
  ensure(const double *R, std::size_t n, const double *box,
         const CachedPairList::Options &opt);

  /// Pool shared by all classical pair potentials.
  static PairListCache &global();

private:
  std::mutex mu_;
  std::vector<std::shared_ptr<CachedPairList>> slots_; // MRU first
};

} // namespace eonc

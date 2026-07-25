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
#include <cstddef>
#include <cstdint>
#include <memory>
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
/// Pair storage is split: pairs with a zero cell shift carry only indices;
/// wrapped pairs additionally carry the precomputed shift offset ``S @ H``
/// (valid until rebuild — any box change forces a rebuild).
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

  struct Pair {
    int32_t i, j;
  };
  struct ShiftPair {
    int32_t i, j;
    double ox, oy, oz; ///< S @ H at build time
  };

  /// True while the cached pair set is valid for positions ``R``: same atom
  /// count, options and box as the build, and max displacement < skin/2.
  [[nodiscard]] bool valid(const double *R, std::size_t n, const double *box,
                           const Options &opt) const;

  /// Rebuild via vesin at ``cutoff + skin``. ``scratch`` supplies the
  /// long-lived vesin buffers (pass the same instance every time so vesin
  /// re-uses its allocations).
  void rebuild(VesinNeighbors &scratch, const double *R, std::size_t n,
               const double *box, const Options &opt);

  /// Visit every cached pair within ``sqrt(cutoff2)`` of the current
  /// positions. ``fn(i, j, dx, dy, dz, r2)`` receives the *historical* eOn
  /// convention ``d = r_i - r_j`` (minimum image / shift applied), i.e. the
  /// negated vesin vector.
  template <typename Fn> void forEach(const double *R, Fn &&fn) const {
    const double cutoff2 = opt_.cutoff * opt_.cutoff;
    for (const auto &p : plain_) {
      const double dx = R[3 * p.i] - R[3 * p.j];
      const double dy = R[3 * p.i + 1] - R[3 * p.j + 1];
      const double dz = R[3 * p.i + 2] - R[3 * p.j + 2];
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 <= cutoff2) {
        fn(p.i, p.j, dx, dy, dz, r2);
      }
    }
    for (const auto &p : shifted_) {
      const double dx = R[3 * p.i] - R[3 * p.j] - p.ox;
      const double dy = R[3 * p.i + 1] - R[3 * p.j + 1] - p.oy;
      const double dz = R[3 * p.i + 2] - R[3 * p.j + 2] - p.oz;
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 <= cutoff2) {
        fn(p.i, p.j, dx, dy, dz, r2);
      }
    }
  }

  [[nodiscard]] std::size_t pairCount() const {
    return plain_.size() + shifted_.size();
  }

private:
  std::vector<Pair> plain_;
  std::vector<ShiftPair> shifted_;
  std::vector<double> Rref_;
  std::array<double, 9> boxref_{};
  Options opt_{};
  std::size_t n_{0};
  bool built_{false};
};

/// Small thread-local pool of CachedPairList slots, matched by geometry
/// proximity. NEB evaluates several images through one shared Potential
/// instance (serially or in parallel); a single cached list would thrash on
/// alternating geometries. Each ``ensure`` call finds the slot whose reference
/// positions are within skin/2 of ``R`` (MRU order) or rebuilds the least
/// recently used slot. thread_local storage keeps this race-free under
/// parallel image evaluation.
class PairListCache {
public:
  /// Slot count covers NEB default image counts; more images degrade to
  /// occasional rebuilds rather than incorrect results.
  static constexpr std::size_t kMaxSlots = 8;

  const CachedPairList &ensure(const double *R, std::size_t n,
                               const double *box,
                               const CachedPairList::Options &opt);

  /// Per-thread pool shared by all classical pair potentials.
  static PairListCache &local();

private:
  std::vector<std::unique_ptr<CachedPairList>> slots_; // MRU first
  VesinNeighbors scratch_;
};

} // namespace eonc

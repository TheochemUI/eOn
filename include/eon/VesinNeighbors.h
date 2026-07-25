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
#include <stdexcept>
#include <string>

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

} // namespace eonc

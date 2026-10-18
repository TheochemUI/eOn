// SPDX-License-Identifier: BSD-3-Clause
#pragma once
#include <array>
#include <cstddef>
#include <vector>

namespace eonc {
struct CutoffPair {
  int first;
  int second;
  std::array<int, 3> shift;
  double distanceSquared;
};
// Fully periodic atom-image rows; d = R_j - R_i + shift * cell.
std::vector<CutoffPair> periodicCutoffNeighbors(const double *positions,
                                                size_t atoms,
                                                const double *rowMajorCell,
                                                double cutoff,
                                                bool full = true);
} // namespace eonc

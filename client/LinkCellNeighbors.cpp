// SPDX-License-Identifier: BSD-3-Clause
#include "eon/CutoffNeighbors.h"
#include <linkcell.hpp>

namespace eonc {
std::vector<CutoffPair> periodicCutoffNeighbors(
    const double *positions, size_t atoms, const double *cell,
    double cutoff, bool full) {
  const auto box = linkcell::Cell::from_vectors(
      {cell[0], cell[1], cell[2]}, {cell[3], cell[4], cell[5]},
      {cell[6], cell[7], cell[8]});
  const auto rows = linkcell::pairs_within(positions, atoms, box, cutoff,
                                          nullptr, 0.0, !full);
  std::vector<CutoffPair> result;
  result.reserve(rows.size());
  for (const auto &pair : rows)
    result.push_back({pair.i, pair.j, pair.shift, pair.dist2});
  return result;
}
} // namespace eonc

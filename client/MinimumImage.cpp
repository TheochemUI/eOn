// SPDX-License-Identifier: BSD-3-Clause
#include "eon/PbcSimd.h"
#include <minimage.hpp>
#include <stdexcept>

namespace eonc::pbc {
void closestDisplacements(double *data, size_t coordinates,
                          const double *rowMajorCell) {
  if (coordinates % 3 != 0)
    throw std::invalid_argument("periodic differences require Cartesian triples");
  if (coordinates == 0) return;
  const auto cell = minimage::Cell::from_con(rowMajorCell);
  for (size_t j = 0; j < coordinates; j += 3) {
    const auto image = cell.displacement_euclidean(
        {0.0, 0.0, 0.0}, {data[j], data[j + 1], data[j + 2]});
    for (size_t k = 0; k < 3; ++k) data[j + k] = image[k];
  }
}
} // namespace eonc::pbc

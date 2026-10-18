/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once

#include <cstddef>

namespace eonc::pbc {

// In-place fractional wrap. Minimum image: x - floor(x + 0.5) in [-0.5, 0.5).
// Legacy unit interval: fmod(x + 1, 1), i.e. (x + 1) - trunc(x + 1).
// Both are elementwise on contiguous storage (row-major AtomMatrix).
void wrapMinimumImage(double *data, size_t n);
void wrapLegacyUnit(double *data, size_t n);

// Euclidean closest images of packed Cartesian triples in a row-major cell.
// Orthogonal half-cell ties retain the negative half-cell representative.
void closestDisplacements(double *data, size_t coordinates,
                          const double *rowMajorCell);

} // namespace eonc::pbc

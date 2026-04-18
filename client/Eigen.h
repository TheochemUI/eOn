/*
 * This file is part of eOn.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Copyright (c) 2010--present, eOn Development Team
 * All rights reserved.
 *
 * Repo:
 * https://github.com/TheochemUI/eOn
 */
#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cassert>
#include <vector>

// Import individual Eigen names rather than `using namespace Eigen` so that
// eOn's row-major matrix aliases (below) don't collide with Eigen's
// column-major defaults.  This keeps eOn binary-compatible with other
// Eigen-based libraries that use the default column-major layout.
using Eigen::Matrix;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::VectorXi;

// Row-major layout so that each row is one atom's (x, y, z) triplet and
// .data() yields [x0, y0, z0, x1, y1, z1, ...] which is what the Fortran
// potentials and VectorXd::Map round-trips expect.
constexpr int eOnStorageOrder = Eigen::RowMajor;

using MatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, eOnStorageOrder>;
using Matrix3d = Eigen::Matrix<double, 3, 3, eOnStorageOrder>;
using Matrix4d = Eigen::Matrix<double, 4, 4, eOnStorageOrder>;
using AtomMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3, eOnStorageOrder>;
using RotationMatrix = Eigen::Matrix<double, 3, 3, eOnStorageOrder>;

// Column-major layout matching Fortran's memory order: (3, nat).
// Row-Major (N,3) and Col-Major (3,N) share the same contiguous memory layout,
// so conversions between AtomMatrix and AtomMatrixF are zero-copy
// reinterpretations.
using AtomMatrixF = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>;

/// SIMD-optimized dot product for contiguous Eigen matrices.
/// Maps both operands as flat VectorXd and uses Eigen's optimized .dot()
/// which leverages SSE/AVX intrinsics, avoiding the element-wise temporary
/// that (a.array() * b.array()).sum() creates.
inline double matDot(const AtomMatrix &a, const AtomMatrix &b) {
  return Eigen::Map<const VectorXd>(a.data(), a.size())
      .dot(Eigen::Map<const VectorXd>(b.data(), b.size()));
}

namespace eonc {

/// Reconstruct AtomMatrix from a flat column-major vector (e.g. from Fortran).
/// RowMajor(N,3) and ColMajor(3,N) share the same contiguous memory layout,
/// so this is effectively a reinterpretation of the flat data.
inline AtomMatrix
from_fortran_layout_vector(const std::vector<double> &flat_colmajor, int nat) {
  assert(flat_colmajor.size() >= static_cast<size_t>(3 * nat));
  return Eigen::Map<const AtomMatrix>(flat_colmajor.data(), nat, 3);
}

} // namespace eonc

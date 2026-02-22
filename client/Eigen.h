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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

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

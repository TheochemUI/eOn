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
using namespace Eigen;

// Row-major layout so that each row is one atom's (x, y, z) triplet and
// .data() yields [x0, y0, z0, x1, y1, z1, ...] which is what the Fortran
// potentials and VectorXd::Map round-trips expect.
constexpr int eOnStorageOrder = Eigen::RowMajor;

typedef Eigen::Matrix<double, Eigen::Dynamic, 3, eOnStorageOrder> AtomMatrix;
typedef Eigen::Matrix<double, 3, 3, eOnStorageOrder> RotationMatrix;

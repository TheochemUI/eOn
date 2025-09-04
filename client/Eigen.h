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

#define EIGEN_DEFAULT_TO_ROW_MAJOR
// #define EIGEN2_SUPPORT
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;

typedef Eigen::Matrix<double, Eigen::Dynamic, 3> AtomMatrix;
typedef Eigen::Matrix<double, 3, 3> RotationMatrix;

/*
** This file is part of eON.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eON Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eON
*/

#pragma once
#define EIGEN_DEFAULT_TO_ROW_MAJOR
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using ScalarType = double;
// Vectors are always column major
using VectorType =
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 1, Eigen::ColMajor>;
template <typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;
template <int E>
using FixedVecType = Eigen::Matrix<ScalarType, E, 1, Eigen::ColMajor>;
template <typename T, int E> using VectorTRC = Eigen::Matrix<T, E, 1>;
// Everything else is RowMajor
using AtomMatrix =
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 3, Eigen::RowMajor>;
using RotationMatrix = Eigen::Matrix<ScalarType, 3, 3, Eigen::RowMajor>;
using MatrixType =
    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
template <typename T>
using Matrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
// Replaces Matrix3d
using Matrix3S = Eigen::Matrix<ScalarType, 3, 3, Eigen::RowMajor>;
template <int E> // Fixed square matrix of ScalarType
using FSqMatS = Eigen::Matrix<ScalarType, E, E, Eigen::RowMajor>;
template <typename T, int R, int C>
using MatrixTRC = Eigen::Matrix<T, R, C, Eigen::RowMajor>;

// Function to map a C-style array to an Eigen matrix
template <typename T>
Matrix<T> cvec_to_mat(const T *array_ptr, int rows, int cols) {
  return Eigen::Map<Matrix<T>>(const_cast<T *>(array_ptr), rows, cols);
}

// Function to map a C-style array to an Eigen vector
template <typename T> Vector<T> cvec_to_vec(const T *array_ptr, int size) {
  return Eigen::Map<Vector<T>>(const_cast<T *>(array_ptr), size);
}

// TODO(rg): Maybe also have a free function for .setZero()

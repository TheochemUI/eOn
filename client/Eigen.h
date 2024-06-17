//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef EIGEN_H
#define EIGEN_H
#define EIGEN_DEFAULT_TO_ROW_MAJOR
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using ScalarType = double;
using AtomMatrix =
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 3, Eigen::RowMajor>;
using RotationMatrix = Eigen::Matrix<ScalarType, 3, 3, Eigen::RowMajor>;
using MatrixType =
    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using VectorType =
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 1, Eigen::RowMajor>;
template <typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::RowMajor>;
template <typename T>
using Matrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
// Replaces Matrix3d
using Matrix3S = Eigen::Matrix<ScalarType, 3, 3, Eigen::RowMajor>;
// Most general replacement
template <typename T, int R, int C>
using MatrixTRC =
    Eigen::Matrix<T, R, C, Eigen::RowMajor>;

#endif

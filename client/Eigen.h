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
// #define EIGEN2_SUPPORT
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;

typedef Eigen::Matrix<double, Eigen::Dynamic, 3> AtomMatrix;
typedef Eigen::Matrix<double, 3, 3> RotationMatrix;
#endif

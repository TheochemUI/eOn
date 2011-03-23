#ifndef EIGEN_H
#define EIGEN_H
#define EIGEN_DEFAULT_TO_ROW_MAJOR
#define EIGEN2_SUPPORT
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
using namespace Eigen;

typedef Eigen::Matrix<double,Eigen::Dynamic,3> AtomMatrix;
#endif

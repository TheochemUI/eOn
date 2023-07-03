#ifndef MATRIXHELPERS_H
#define MATRIXHELPERS_H

#include "Matter.h"
#include "Parameters.h"

namespace helper_functions {

/**
 * \brief Check two eigen objects for equality
 *
 * @param Eigen::MatrixBase the underlying base class
 */
template <typename T>
bool eigenEquality(const Eigen::MatrixBase<T> &lhs,
                   const Eigen::MatrixBase<T> &rhs,
                   const double threshold = 1e-4) {
  return lhs.isApprox(rhs, threshold);
}
} // namespace helper_functions
#endif /* MATRIXHELPERS_H */

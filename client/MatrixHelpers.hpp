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

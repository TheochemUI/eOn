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
#include "eon/SurrogatePotential.h"

#include <span>

std::tuple<double, AtomMatrix, double>
SurrogatePotential::get_ef_var(const AtomMatrix pos, const VectorXi atmnrs,
                               const Matrix3d box) {
  double energy{std::numeric_limits<double>::infinity()};
  long nAtoms = static_cast<long>(pos.rows());
  AtomMatrix forces{MatrixXd::Zero(nAtoms, 3)};
  double var{0};
  const auto n = static_cast<size_t>(nAtoms);
  this->force(std::span<const double>(pos.data(), n * 3),
              std::span<const int>(atmnrs.data(), n),
              std::span<double>(forces.data(), n * 3), &energy, &var,
              std::span<const double>(box.data(), 9));
  return std::make_tuple(energy, forces, var);
};

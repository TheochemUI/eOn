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
#include "SurrogatePotential.h"

std::tuple<double, AtomMatrix, double>
SurrogatePotential::get_ef_var(const AtomMatrix pos, const VectorXi atmnrs,
                               const Matrix3d box) {
  double energy{std::numeric_limits<double>::infinity()};
  long nAtoms{pos.rows()};
  AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
  double var{0};
  this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy, &var,
              box.data());
  return std::make_tuple(energy, forces, var);
};

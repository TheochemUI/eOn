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

#include "Potential.h"

class SurrogatePotential : public Potential {

public:
  SurrogatePotential(PotType a_ptype, std::shared_ptr<Parameters> a_params)
      : Potential(a_ptype, a_params) {}
  virtual ~SurrogatePotential() = default;
  std::tuple<double, AtomMatrix, double> // energy, forces, energy variance
  get_ef_var(const AtomMatrix pos, const VectorXi atmnrs, const Matrix3d box);
  virtual void train_optimize(Eigen::MatrixXd a_features,
                              Eigen::MatrixXd a_targets) = 0;
};

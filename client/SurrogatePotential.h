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
  SurrogatePotential(PotType a_ptype)
      : Potential(a_ptype) {}
  virtual ~SurrogatePotential() = default;
  std::tuple<double, AtomMatrix, double> // energy, forces, energy variance
  get_ef_var(const AtomMatrix pos, const Vector<int> atmnrs,
             const Matrix3S box);
  virtual void train_optimize(MatrixType a_features, MatrixType a_targets) = 0;
};

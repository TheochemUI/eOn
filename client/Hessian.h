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
#include "Eigen.h"
#include "EonLogger.h"

#include "Matter.h"
#include "Parameters.h"

namespace eonc {

class Hessian {
public:
  Hessian(const Parameters &params, Matter *matter);
  ~Hessian() = default;

  MatrixXd getHessian(Matter *matterIn, const VectorXi &atomsIn);
  VectorXd getFreqs(Matter *matterIn, const VectorXi &atomsIn);
  //    VectorXd getModes(Matter *matterIn, VectorXi atomsIn);
  VectorXd removeZeroFreqs(const VectorXd &freqs);

private:
  Matter *matter;
  const Parameters &parameters;

  MatrixXd hessian;
  VectorXd freqs;

  VectorXi atoms;
  bool calculate();
  eonc::log::Scoped log;
};

} // namespace eonc

using eonc::Hessian;

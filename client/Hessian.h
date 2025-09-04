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

#include "Matter.h"
#include "Parameters.h"

class Hessian {
public:
  Hessian(Parameters *params, Matter *matter);
  ~Hessian();

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> getHessian(Matter *matterIn,
                                                            VectorXi atomsIn);
  VectorXd getFreqs(Matter *matterIn, VectorXi atomsIn);
  //    VectorXd getModes(Matter *matterIn, VectorXi atomsIn);
  VectorXd removeZeroFreqs(VectorXd freqs);

private:
  Matter *matter;
  Parameters *parameters;

  Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessian;
  //    VectorXd modes;
  VectorXd freqs;

  VectorXi atoms;
  bool calculate(void);
  std::shared_ptr<spdlog::logger> log;
};

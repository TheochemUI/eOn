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

#include "Matter.h"

namespace eonc {

class ObjectiveFunction {
protected:
  const Parameters &params;

public:
  ObjectiveFunction(const Parameters &paramsPassed)
      : params{paramsPassed} {}
  virtual ~ObjectiveFunction() {}
  virtual double getEnergy() = 0;
  virtual VectorXd getGradient(bool fdstep = false) = 0;
  virtual void setPositions(const VectorXd &x) = 0;
  virtual VectorXd getPositions() = 0;
  virtual int degreesOfFreedom() = 0;
  virtual bool isConverged() = 0;
  virtual double getConvergence() = 0;
  virtual VectorXd difference(const VectorXd &a, const VectorXd &b) = 0;
};

} // namespace eonc

using eonc::ObjectiveFunction;

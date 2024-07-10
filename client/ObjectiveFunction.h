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

namespace eonc {
class ObjectiveFunction {
public:
  ObjectiveFunction() {}
  virtual ~ObjectiveFunction() {}
  virtual double getEnergy() const = 0;
  virtual VectorType getGradient(bool fdstep = false) const = 0;
  virtual void setPositions(VectorType x) = 0;
  virtual VectorType getPositions() const = 0;
  virtual int degreesOfFreedom() const = 0;
  virtual bool isConverged() const = 0;
  virtual double getConvergence() const = 0;
  virtual VectorType difference(VectorType a, VectorType b) = 0;
};

} // namespace eonc

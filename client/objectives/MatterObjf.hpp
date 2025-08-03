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

#include "client/BaseStructures.h"
#include "client/ObjectiveFunction.h"
#include "client/matter/Matter.h"

namespace eonc::objf {
class MatterObjectiveFunction : public ObjectiveFunction {
public:
  struct Params {
    ConvergenceMeasure optCM{ConvergenceMeasure::NORM};
    ScalarType optConvergedForce{1e-3};
  };

  MatterObjectiveFunction(const Params &_params, const Matter &_matter)
      : m_p(std::cref(_params)),
        matter(std::cref(_matter)) {}

  ~MatterObjectiveFunction() override = default;

  double getEnergy() const override;
  VectorType getGradient(bool = false) const override;
  void setPositions(const VectorType &) override;
  VectorType getPositions() const override;
  int degreesOfFreedom() const override;
  bool isConverged() const override;
  double getConvergence() const override;
  VectorType difference(const VectorType &, const VectorType &) const override;

private:
  const Params m_p;
  const Matter &matter;
};
} // namespace eonc::objf

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

#include "client/ObjectiveFunction.h"
#include "client/matter/Matter.h"
namespace eonc {
class MatterObjectiveFunction : public ObjectiveFunction {
public:
  struct Params {
    std::string optConvergenceMetric{"norm"};
    double optConvergedForce{1e-4};
  };

  MatterObjectiveFunction(const Params &_params, const Matter &_matter)
      : convMetric(_params.optConvergenceMetric),
        convForce(_params.optConvergedForce),
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
  const std::string convMetric;
  const double convForce;
  const Matter &matter;
};

class RelaxJob {
public:
  RelaxJob() { m_log = spdlog::get("combi"); }
  ~RelaxJob(void) = default;
  bool runImpl(Matter &);

private:
  std::shared_ptr<spdlog::logger> m_log;
};

} // namespace eonc

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

#include "Job.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
namespace eonc {

class MinObjF : public ObjectiveFunction {
private:
  Matter &_mat;
  ObjectiveFunction::Conv _metric{ObjectiveFunction::Conv::NORM};
  double _convForce = 1e-6;

public:
  MinObjF(Matter &mat_a)
      : ObjectiveFunction(),
        _mat{mat_a} {}
  MinObjF(Matter &mat_a, ObjectiveFunction::Conv metric_a, double cforce_a)
      : ObjectiveFunction(),
        _mat{mat_a},
        _metric{metric_a},
        _convForce{cforce_a} {}
  ~MinObjF() = default;
  double getEnergy() const override;
  VectorType getGradient(bool fdstep = false) const override;
  VectorType getPositions() const override;
  int degreesOfFreedom() const override;
  bool isConverged() const override;
  double getConvergence() const override;
  void setPositions(VectorType x) override;
  VectorType difference(VectorType a, VectorType b) override;
};

class MinimizationJob : public Job {
public:
  MinimizationJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)),
        fcalls{0} {
    log = spdlog::get("combi");
  }
  ~MinimizationJob(void) = default;
  std::vector<std::string> run(void);

private:
  size_t fcalls;
  RunStatus status;
  std::shared_ptr<spdlog::logger> log;
};

} // namespace eonc

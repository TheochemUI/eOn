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
#include "EonLogger.h"

#include "Eigen.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Parameters.h"

namespace eonc {

class SteepestDescent final : public Optimizer {
public:
  SteepestDescent(std::shared_ptr<ObjectiveFunction> a_objf,
                  const Parameters &a_params)
      : Optimizer(std::move(a_objf), OptType::SD,
                  OptimizerConfig::fromParams(a_params)),
        iteration{0} {}
  ~SteepestDescent() = default;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;

private:
  eonc::log::FileScoped m_log{"sd", "_sd.log"};
  Eigen::VectorXd getStep(const Eigen::VectorXd &a_f);
  size_t iteration{0};
  Eigen::VectorXd m_rPrev;
  Eigen::VectorXd m_fPrev;
};

} // namespace eonc

using eonc::SteepestDescent;

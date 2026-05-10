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

#include "HelperFunctions.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Parameters.h"

#include <deque>

namespace eonc {

/// Curvature-update gate. \f$|s_0 \cdot y_0|\f$ below this means the
/// gradient barely changed across the step; pushing it through
/// \f$\rho = 1 / (s_0 \cdot y_0)\f$ would amplify denormals into the
/// L-BFGS history. 1e-30 is well below any chemistry-relevant
/// value and well above double-precision underflow.
inline constexpr double LBFGS_EPS = 1e-30;

class LBFGS final : public Optimizer {

public:
  LBFGS(std::shared_ptr<ObjectiveFunction> a_objf, const Parameters &a_params)
      : Optimizer(std::move(a_objf), OptType::LBFGS,
                  OptimizerConfig::fromParams(a_params)),
        m_iteration{0},
        m_memory{std::min(
            m_objf->degreesOfFreedom(),
            static_cast<int>(a_params.optimizer_options.lbfgs.memory))} {}

  ~LBFGS() = default;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;
  int update(const Eigen::VectorXd &a_r1, const Eigen::VectorXd &a_r0,
             const Eigen::VectorXd &a_f1, const Eigen::VectorXd &a_f0);
  void reset(void);

private:
  Eigen::VectorXd getStep(double a_maxMove, const Eigen::VectorXd &a_f);

  int m_iteration;
  int m_memory;

  std::deque<Eigen::VectorXd> m_s;
  std::deque<Eigen::VectorXd> m_y;
  std::deque<double> m_rho;

  Eigen::VectorXd m_rPrev;
  Eigen::VectorXd m_fPrev;
  eonc::log::FileScoped m_log{"lbfgs", "_lbfgs.log"};
};

} // namespace eonc

using eonc::LBFGS;

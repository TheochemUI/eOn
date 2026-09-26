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

class LBFGS final : public Optimizer {

public:
  LBFGS(std::shared_ptr<ObjectiveFunction> a_objf, const Parameters &a_params)
      : Optimizer(a_objf, OptType::LBFGS,
                  OptimizerConfig::fromParams(a_params)),
        m_memory{std::min(
            a_objf->degreesOfFreedom(),
            static_cast<int>(a_params.optimizer_options().lbfgs.memory))},
        m_rPrev{Eigen::VectorXd::Zero(a_objf->degreesOfFreedom())},
        m_fPrev{Eigen::VectorXd::Zero(a_objf->degreesOfFreedom())} {}

  ~LBFGS() = default;

  [[nodiscard]] int step(double a_maxMove) override;
  [[nodiscard]] int run(size_t a_maxIterations, double a_maxMove) override;
  [[nodiscard]] int update(const Eigen::VectorXd &a_r1,
                           const Eigen::VectorXd &a_r0,
                           const Eigen::VectorXd &a_f1,
                           const Eigen::VectorXd &a_f0, double a_e1,
                           double a_e0);
  void reset();

private:
  [[nodiscard]] Eigen::VectorXd getStep(double a_maxMove,
                                        const Eigen::VectorXd &a_f);
  [[nodiscard]] Eigen::Vector3d micRij(const Eigen::VectorXd &pos, int i,
                                       int j) const;
  [[nodiscard]] bool usesPrecon() const;
  [[nodiscard]] Eigen::MatrixXd buildPrecon(const Eigen::VectorXd &pos) const;
  [[nodiscard]] Eigen::VectorXd applyH0(const Eigen::VectorXd &q, double H0,
                                        const Eigen::VectorXd &pos) const;
  [[nodiscard]] Eigen::VectorXd hessianStep(double a_maxMove,
                                            const Eigen::VectorXd &a_f);

  int m_iteration{0};
  int m_memory{0};

  std::deque<Eigen::VectorXd> m_s;
  std::deque<Eigen::VectorXd> m_y;
  std::deque<double> m_rho;
  std::deque<double> m_eHist;

  Eigen::VectorXd m_rPrev;
  Eigen::VectorXd m_fPrev;
  double m_ePrev{0.0};
  eonc::log::FileScoped m_log{"lbfgs", "_lbfgs.log"};
};

} // namespace eonc

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

#define LBFGS_EPS 1e-30

class LBFGS final : public Optimizer {

public:
  LBFGS(std::shared_ptr<ObjectiveFunction> a_objf, const Parameters &a_params)
      : Optimizer(a_objf, OptType::LBFGS, a_params),
        m_iteration{0},
        m_memory{std::min(
            a_objf->degreesOfFreedom(),
            static_cast<int>(a_params.optimizer_options.lbfgs.memory))} {}

  ~LBFGS() = default;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;
  int update(const Eigen::VectorXd &a_r1, const Eigen::VectorXd &a_r0,
             const Eigen::VectorXd &a_f1, const Eigen::VectorXd &a_f0,
             double a_e1, double a_e0);
  void reset(void);

private:
  Eigen::VectorXd getStep(double a_maxMove, const Eigen::VectorXd &a_f);
  Eigen::Vector3d micRij(const Eigen::VectorXd &pos, int i, int j) const;
  Eigen::MatrixXd buildPrecon(const Eigen::VectorXd &pos) const;
  Eigen::VectorXd applyH0(const Eigen::VectorXd &q, double H0,
                          const Eigen::VectorXd &pos) const;

  int m_iteration;
  int m_memory;

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

using eonc::LBFGS;

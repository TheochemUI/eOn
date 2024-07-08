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

#include "HelperFunctions.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Parameters.h"

#define LBFGS_EPS 1e-30

class LBFGS final : public Optimizer {

public:
  LBFGS(std::shared_ptr<ObjectiveFunction> a_objf,
        std::shared_ptr<Parameters> a_params)
      : Optimizer(a_objf, OptType::LBFGS, a_params),
        m_iteration{0},
        m_memory(std::min(a_objf->degreesOfFreedom(),
                          static_cast<int>(a_params->optim.LBFGSMemory))) {

    if (spdlog::get("lbfgs")) {
      m_log = spdlog::get("lbfgs");
    } else {
      m_log = spdlog::basic_logger_mt("lbfgs", "_lbfgs.log", true);
    }
    m_log->set_pattern("[%l] [LBFGS] %v");
  }

  ~LBFGS() = default;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;
  int update(VectorType a_r1, VectorType a_r0, VectorType a_f1,
             VectorType a_f0);
  void reset(void);

private:
  VectorType getStep(double a_maxMove, VectorType a_f);

  int m_iteration;
  int m_memory;

  std::vector<VectorType> m_s;
  std::vector<VectorType> m_y;
  std::vector<double> m_rho;

  VectorType m_rPrev;
  VectorType m_fPrev;
  std::shared_ptr<spdlog::logger> m_log;
};

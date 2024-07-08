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
#ifndef SteepestDescent_H
#define SteepestDescent_H

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Parameters.h"

class SteepestDescent final : public Optimizer {
public:
  SteepestDescent(std::shared_ptr<ObjectiveFunction> a_objf,
                  std::shared_ptr<Parameters> a_params)
      : Optimizer(a_objf, OptType::SD, a_params),
        iteration{0} {
    if (spdlog::get("sd")) {
      m_log = spdlog::get("sd");
    } else {
      m_log = spdlog::basic_logger_mt("sd", "_sd.log", true);
    }
    m_log->set_pattern("[%l] [SD] %v");
  }
  ~SteepestDescent() = default;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;

private:
  std::shared_ptr<spdlog::logger> m_log;
  VectorType getStep(VectorType a_f);
  size_t iteration;
  VectorType m_rPrev;
  VectorType m_fPrev;
};

#endif

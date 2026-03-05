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
#include "HelperFunctions.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Parameters.h"

class SteepestDescent final : public Optimizer {
public:
  SteepestDescent(std::shared_ptr<ObjectiveFunction> a_objf,
                  const Parameters &a_params)
      : Optimizer(a_objf, OptType::SD, a_params),
        iteration{0} {
    m_log = quill::Frontend::create_or_get_logger(
        "sd",
        quill::Frontend::create_or_get_sink<quill::FileSink>(
            "_sd.log",
            []() {
              quill::FileSinkConfig cfg;
              cfg.set_open_mode('w');
              return cfg;
            }(),
            quill::FileEventNotifier{}),
        quill::PatternFormatterOptions{"%(message)"},
        quill::ClockSourceType::System);
  }
  ~SteepestDescent() = default;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;

private:
  quill::Logger *m_log{nullptr};
  Eigen::VectorXd getStep(Eigen::VectorXd a_f);
  size_t iteration;
  Eigen::VectorXd m_rPrev;
  Eigen::VectorXd m_fPrev;
};

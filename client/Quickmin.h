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

#include "Matter.h"
#include "Optimizer.h"
#include "Parameters.h"

// when changing away from final, remember to mark the destructor as virtual
class Quickmin final : public Optimizer {

public:
  Quickmin(std::shared_ptr<ObjectiveFunction> a_objf,
           const Parameters &a_params)
      : Optimizer(a_objf, OptType::QM, a_params),
        m_dt{a_params.optimizer_options.time_step},
        m_dt_max{a_params.optimizer_options.max_time_step},
        m_max_move{a_params.optimizer_options.max_move},
        m_vel{Eigen::VectorXd::Zero(a_objf->degreesOfFreedom())},
        m_iteration{0},
        m_max_iter{a_params.optimizer_options.max_iterations} {
    m_log = quill::Frontend::create_or_get_logger(
        "qm",
        quill::Frontend::create_or_get_sink<quill::FileSink>(
            "_qm.log",
            []() {
              quill::FileSinkConfig cfg;
              cfg.set_open_mode('w');
              return cfg;
            }(),
            quill::FileEventNotifier{}),
        quill::PatternFormatterOptions{"%(message)"},
        quill::ClockSourceType::System);
  }
  ~Quickmin() = default;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;

private:
  double m_dt, m_dt_max, m_max_move;
  Eigen::VectorXd m_vel;
  size_t m_iteration, m_max_iter;
  quill::Logger *m_log{nullptr};
};

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
#include "Optimizer.h"
#include "Parameters.h"

class FIRE : public Optimizer {

public:
  FIRE(std::shared_ptr<ObjectiveFunction> a_objf, const Parameters &a_params)
      : Optimizer(a_objf, OptType::FIRE, a_params),
        m_dt{a_params.optimizer_options.time_step},
        m_dt_max{a_params.optimizer_options.max_time_step},
        m_max_move{a_params.optimizer_options.max_move},
        m_N_min{5},
        m_N{0},
        m_vel{Eigen::VectorXd::Zero(a_objf->degreesOfFreedom())},
        m_alpha_start{0.1},
        m_alpha{m_alpha_start},
        m_f_inc{1.1},
        m_f_dec{0.5},
        m_f_a{0.99},
        m_iteration{0} {
    m_log = quill::Frontend::create_or_get_logger(
        "fire",
        quill::Frontend::create_or_get_sink<quill::FileSink>(
            "_fire.log",
            []() {
              quill::FileSinkConfig cfg;
              cfg.set_open_mode('w');
              return cfg;
            }(),
            quill::FileEventNotifier{}),
        quill::PatternFormatterOptions{"%(message)"},
        quill::ClockSourceType::System);
  }
  virtual ~FIRE() = default;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;

private:
  double m_dt, m_dt_max, m_max_move;
  size_t m_N_min, m_N;
  Eigen::VectorXd m_vel;
  double m_alpha_start;
  double m_alpha;
  double m_f_inc;
  double m_f_dec;
  double m_f_a;
  size_t m_iteration;
  quill::Logger *m_log{nullptr};
};

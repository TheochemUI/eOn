#ifndef FIRE_H
#define FIRE_H

#include "Matter.h"
#include "Optimizer.h"
#include "Parameters.h"

class FIRE : public Optimizer {

public:
  FIRE(std::shared_ptr<ObjectiveFunction> a_objf,
       std::shared_ptr<Parameters> a_params)
      : m_objf{a_objf},
        m_params{a_params},
        m_dt{a_params->optTimeStep},
        m_dt_max{a_params->optMaxTimeStep},
        m_max_move{a_params->optMaxMove},
        m_N_min{5},
        m_N{0},
        m_v{Eigen::VectorXd::Zero(a_objf->degreesOfFreedom())},
        m_alpha_start{0.1},
        m_alpha{m_alpha_start},
        m_f_inc{1.1},
        m_f_dec{0.5},
        m_f_a{0.99},
        m_iteration{0} {
    if (spdlog::get("fire")) {
      m_log = spdlog::get("fire");
    } else {
      m_log = spdlog::basic_logger_st("fire", "_fire.log", true);
    }
    m_log->set_pattern("[%l] [FIRE] %v");
  }
  ~FIRE() = default;

  int step(double a_maxMove);
  int run(int a_maxIterations, double a_maxMove);

private:
  std::shared_ptr<ObjectiveFunction> m_objf;
  std::shared_ptr<Parameters> m_params;
  double m_dt, m_dt_max, m_max_move;
  size_t m_N_min, m_N;
  Eigen::VectorXd m_v;
  double m_alpha_start;
  double m_alpha;
  double m_f_inc;
  double m_f_dec;
  double m_f_a;
  size_t m_iteration;
  shared_ptr<spdlog::logger> m_log;
};

#endif

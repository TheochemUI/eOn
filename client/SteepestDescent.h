#ifndef SteepestDescent_H
#define SteepestDescent_H

#include "Eigen.h"
#include "HelperFunctions.h"
#include "Matter.h"
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
  shared_ptr<spdlog::logger> m_log;
  Eigen::VectorXd getStep(Eigen::VectorXd a_f);
  size_t iteration;
  Eigen::VectorXd m_rPrev;
  Eigen::VectorXd m_fPrev;
};

#endif

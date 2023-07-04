#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "HelperFunctions.h"
#include "Matter.h"
#include "Parameters.h"

// dimer method to find the lowest curvature mode
class MonteCarlo {

public:
  MonteCarlo(std::shared_ptr<Matter> const matterIn,
             std::shared_ptr<Parameters> paramsIn)
      : matter{matterIn}, params{paramsIn} {
    log = spdlog::get("combi");
  }
  ~MonteCarlo() = default;

  void run(int numSteps, double temperature, double stepSize);

private:
  std::shared_ptr<Matter> matter;
  std::shared_ptr<Parameters> params;
  shared_ptr<spdlog::logger> log;
};

#endif

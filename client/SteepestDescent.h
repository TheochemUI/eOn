#ifndef SteepestDescent_H
#define SteepestDescent_H

#include "Eigen.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Parameters.h"

class SteepestDescent : public Optimizer {

public:
  SteepestDescent(ObjectiveFunction *objf, Parameters *parameters);
  ~SteepestDescent() = default;

  int step(double maxMove);
  int run(int maxIterations, double maxMove);

private:
  shared_ptr<spdlog::logger> log;
  VectorXd getStep(VectorXd f);
  Parameters *parameters;
  ObjectiveFunction *objf;

  int iteration;

  VectorXd rPrev;
  VectorXd fPrev;
};

#endif

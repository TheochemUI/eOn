#ifndef LBFGS_H
#define LBFGS_H

#include "Eigen.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Parameters.h"
#include <vector>

#define LBFGS_EPS 1e-30

class LBFGS : public Optimizer {

public:
  LBFGS(ObjectiveFunction *objf, Parameters *parameters);
  ~LBFGS() = default;

  int step(double maxMove);
  int run(int maxIterations, double maxMove);
  int update(VectorXd r1, VectorXd r0, VectorXd f1, VectorXd f0);
  void reset(void);

private:
  VectorXd getStep(double maxMove, VectorXd f);
  Parameters *parameters;
  ObjectiveFunction *objf;

  int iteration;
  int memory;

  std::vector<VectorXd> s;
  std::vector<VectorXd> y;
  std::vector<double> rho;

  VectorXd rPrev;
  VectorXd fPrev;
  shared_ptr<spdlog::logger> log;
};

#endif

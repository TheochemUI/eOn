#ifndef FIRE_H
#define FIRE_H

#include "Matter.h"
#include "Optimizer.h"
#include "Parameters.h"

class FIRE : public Optimizer {

public:
  double dt;

  FIRE(ObjectiveFunction *objf, Parameters *parameters);
  ~FIRE() = default;

  int step(double maxMove);
  int run(int maxIterations, double maxMove);

private:
  ObjectiveFunction *objf;
  Parameters *parameters;
  VectorXd v;
  double alpha;
  double alpha_start;
  int N, N_min;
  double f_inc;
  double f_dec;
  double f_a;
  double dt_max;
  int iteration;
  shared_ptr<spdlog::logger> log;
};

#endif

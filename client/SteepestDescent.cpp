
// Based on the SteepestDescent minimizer written in ASE.

#include "SteepestDescent.h"

SteepestDescent::SteepestDescent(ObjectiveFunction *objfPassed,
                                 Parameters *parametersPassed) {
  objf = objfPassed;
  parameters = parametersPassed;

  iteration = 0;
  log = spdlog::basic_logger_st("sd", "_sd.log", true);
  log->set_pattern("[%l] [SD] %v");
}

int SteepestDescent::step(double maxMove) {
  VectorXd r = objf->getPositions();
  VectorXd f = -objf->getGradient();

  VectorXd dr;
  double alpha = parameters->optSDAlpha;
  if (parameters->optSDTwoPoint == true && iteration > 0) {
    VectorXd dx = r - rPrev;
    VectorXd dg = -f + fPrev;
    alpha = dx.dot(dx) / dx.dot(dg);
    if (alpha < 0) {
      alpha = parameters->optSDAlpha;
    }
    SPDLOG_LOGGER_DEBUG(log, "[SD] alpha: {:.4e}", alpha);
  }

  dr = alpha * f;
  dr = helper_functions::maxAtomMotionAppliedV(dr, maxMove);

  objf->setPositions(r + dr);

  rPrev = r;
  fPrev = f;

  iteration++;

  //    return objf->isConverged();
  if (objf->isConverged())
    return 1;
  return 0;
}

int SteepestDescent::run(int maxSteps, double maxMove) {
  while (!objf->isConverged() && iteration < maxSteps) {
    step(maxMove);
  }
  //    return objf->isConverged();
  if (objf->isConverged())
    return 1;
  return 0;
}

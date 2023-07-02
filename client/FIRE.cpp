#include "FIRE.h"
#include "HelperFunctions.h"
FIRE::FIRE(ObjectiveFunction *objfPassed, Parameters *parametersPassed) {
  objf = objfPassed;
  parameters = parametersPassed;
  dt = parametersPassed->optTimeStep;
  dt_max = parametersPassed->optMaxTimeStep;
  N_min = 5;
  N = 0;
  f_inc = 1.1;
  f_dec = 0.5;
  alpha_start = 0.1;
  alpha = alpha_start;
  f_a = 0.99;
  v.resize(objf->degreesOfFreedom());
  v.setZero();
  iteration = 0;
  log = spdlog::basic_logger_st("fire", "_fire.log", true);
  log->set_pattern("%v");
}

FIRE::~FIRE() { return; }

int FIRE::step(double maxMove) {
  double P = 0;
  // Check convergence.
  if (objf->isConverged()) {
    return 1;
  }

  // Velocity Verlet
  VectorXd f = -objf->getGradient();
  VectorXd x = objf->getPositions();

  v += f * dt;
  VectorXd dx = v * dt;

  dx = helper_functions::maxAtomMotionAppliedV(dx, parameters->optMaxMove);
  objf->setPositions(x + dx);

  f = -objf->getGradient();
  VectorXd f_unit = f / f.norm();

  // FIRE
  P = f.dot(v);
  v = (1 - alpha) * v + alpha * f_unit * v.norm();
  // SPDLOG_LOGGER_DEBUG(log, "[FIRE] P: {:.4f}, v: {:.4f}, dt: {:.4f}, alpha:
  // {:.4f}, N: {}", P, v, dt, alpha, N);
  if (P >= 0) {
    N++;
    if (N > N_min) {
      dt = min(dt * f_inc, dt_max);
      alpha = alpha * f_a;
    }
  } else {
    dt = dt * f_dec;
    v = v * 0.0;
    alpha = alpha_start;
    N = 0;
  }

  // add a sanity check on dt
  if (dt < 1e-6) {
    SPDLOG_LOGGER_CRITICAL(log, "[FIRE] [critical] dt is too small: {:.4f}",
                           dt);
    std::exit(1);
  }

  iteration++;
  return objf->isConverged();
}

int FIRE::run(int maxSteps, double maxMove) {
  while (!objf->isConverged() && iteration < maxSteps) {
    step(maxMove);
  }
  return objf->isConverged();
}

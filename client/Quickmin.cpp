#include "Quickmin.h"
#include "HelperFunctions.h"

Quickmin::Quickmin(ObjectiveFunction *objfPassed,
                   Parameters *parametersPassed) {
  objf = objfPassed;
  parameters = parametersPassed;
  dt = parametersPassed->optTimeStep;
  velocity.resize(objf->degreesOfFreedom());
  velocity.setZero();
  iteration = 0;
  if (spdlog::get("qm")) {
    log = spdlog::get("qm");
  } else {
    log = spdlog::basic_logger_st("qm", "_qm.log", true);
  }
  log->set_pattern("[%l] [QM] %v");
}

Quickmin::~Quickmin() { return; }

int Quickmin::step(double maxMove) {
  VectorXd force = -objf->getGradient();
  if (parameters->optQMSteepestDecent) {
    velocity.setZero();
  } else {
    if (velocity.dot(force) < 0) {
      velocity.setZero();
    } else {
      VectorXd f_unit = force / force.norm();
      velocity = velocity.dot(f_unit) * f_unit;
    }
  }

  velocity += force * dt;
  VectorXd dr = helper_functions::maxAtomMotionAppliedV(velocity * dt,
                                                        parameters->optMaxMove);
  // SPDLOG_LOGGER_INFO(log, "{} Velocity is {}", iteration,
  // fmt::streamed(velocity));
  objf->setPositions(objf->getPositions() + dr);
  iteration++;
  return objf->isConverged();
}

int Quickmin::run(int maxSteps, double maxMove) {
  while (!objf->isConverged() && iteration < maxSteps) {
    step(maxMove);
  }
  return objf->isConverged();
}

#include "Optimizer.h"
#include "ConjugateGradients.h"
#include "FIRE.h"
#include "LBFGS.h"
#include "Quickmin.h"
#include "SteepestDescent.h"

Optimizer *Optimizer::getOptimizer(ObjectiveFunction *objf,
                                   Parameters *parameters) {
  Optimizer *mizer = NULL;
  if (parameters->optMethod == "cg") {
    mizer = new ConjugateGradients(objf, parameters);
  } else if (parameters->optMethod == "qm") {
    mizer = new Quickmin(objf, parameters);
  } else if (parameters->optMethod == "lbfgs") {
    mizer = new LBFGS(objf, parameters);
  } else if (parameters->optMethod == "fire") {
    mizer = new FIRE(objf, parameters);
  } else if (parameters->optMethod == "sd") {
    mizer = new SteepestDescent(objf, parameters);
  } else {
    auto log = spdlog::get("_traceback");
    SPDLOG_LOGGER_CRITICAL(log, "Unknown optMethod: {}", parameters->optMethod);
    std::exit(1);
  }
  return mizer;
}

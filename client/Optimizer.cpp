#include "Optimizer.h"
#include "ConjugateGradients.h"
#include "FIRE.h"
#include "LBFGS.h"
#include "Quickmin.h"
#include "SteepestDescent.h"

Optimizer *Optimizer::getOptimizer(ObjectiveFunction *objf,
                                   Parameters *parameters, bool refine) {
  Optimizer *mizer = nullptr;
  std::string meth = "NONE"s;
  if (refine) {
    if (parameters->refineOptMethod == "NONE"s) {
      SPDLOG_CRITICAL(
          "refine was passed to getOptimizer when it shouldn't have been");
      std::exit(1);
      meth = parameters->optMethod;
    } else {
      meth = parameters->refineOptMethod;
    }
  } else {
    meth = parameters->optMethod;
  }
  if (meth == "cg") {
    mizer = new ConjugateGradients(objf, parameters);
  } else if (meth == "qm") {
    mizer = new Quickmin(objf, parameters);
  } else if (meth == "lbfgs") {
    mizer = new LBFGS(objf, parameters);
  } else if (meth == "fire") {
    mizer = new FIRE(objf, parameters);
  } else if (meth == "sd") {
    mizer = new SteepestDescent(objf, parameters);
  } else {
    auto log = spdlog::get("_traceback");
    SPDLOG_LOGGER_CRITICAL(log, "Unknown optMethod: {}", parameters->optMethod);
    std::exit(1);
  }
  return mizer;
}

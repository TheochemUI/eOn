#include "Optimizer.h"
#include "ConjugateGradients.h"
#include "FIRE.h"
#include "LBFGS.h"
#include "Quickmin.h"
#include "SteepestDescent.h"

std::unique_ptr<Optimizer>
Optimizer::getOptimizer(std::shared_ptr<ObjectiveFunction> a_objf,
                        std::shared_ptr<Parameters> a_params, bool a_refine) {
  std::unique_ptr<Optimizer> mizer = nullptr;
  std::string meth = "NONE"s;
  if (a_refine) {
    if (a_params->refineOptMethod == "NONE"s) {
      auto log = spdlog::get("_traceback");
      SPDLOG_LOGGER_CRITICAL(
          log, "refine was passed to getOptimizer when it shouldn't have been");
      std::exit(1);
      meth = a_params->optMethod;
    } else {
      meth = a_params->refineOptMethod;
    }
  } else {
    meth = a_params->optMethod;
  }
  if (meth == "cg") {
    mizer = std::make_unique<ConjugateGradients>(a_objf, a_params);
  } else if (meth == "qm") {
    mizer = std::make_unique<Quickmin>(a_objf, a_params);
  } else if (meth == "lbfgs") {
    mizer = std::make_unique<LBFGS>(a_objf, a_params);
  } else if (meth == "fire") {
    mizer = std::make_unique<FIRE>(a_objf, a_params);
  } else if (meth == "sd") {
    mizer = std::make_unique<SteepestDescent>(a_objf, a_params);
  } else {
    auto log = spdlog::get("_traceback");
    SPDLOG_LOGGER_CRITICAL(log, "Unknown optMethod or a_refineOptMethod: {}",
                           meth);
    std::exit(1);
  }
  return mizer;
}

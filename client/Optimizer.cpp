#include "Optimizer.h"
#include "BaseStructures.h"
#include "ConjugateGradients.h"
#include "FIRE.h"
#include "LBFGS.h"
#include "Quickmin.h"
#include "SteepestDescent.h"

namespace helpers::create {
std::unique_ptr<Optimizer> mkOptim(std::shared_ptr<ObjectiveFunction> a_objf,
                                   OptType a_otype,
                                   std::shared_ptr<Parameters> a_params) {
  switch (a_otype) {
  case OptType::FIRE: {
    return std::make_unique<FIRE>(a_objf, a_params);
  }
  case OptType::QM: {
    return std::make_unique<Quickmin>(a_objf, a_params);
  }
  case OptType::CG: {
    return std::make_unique<ConjugateGradients>(a_objf, a_params);
  }
  case OptType::LBFGS: {
    return std::make_unique<LBFGS>(a_objf, a_params);
  }
  case OptType::SD: {
    return std::make_unique<SteepestDescent>(a_objf, a_params);
  }
  case OptType::None: {
    throw std::runtime_error("[Optimizer] Cannot create None");
  }
  case OptType::Unknown: {
    throw std::runtime_error("[Optimizer] Cannot create Unknown");
  }
  default: {
    throw std::runtime_error("Unsupported optimization type");
  }
  }
}
} // namespace helpers::create

/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "eon/Optimizer.h"
#include "eon/BaseStructures.h"
#include "eon/ConjugateGradients.h"
#include "eon/FIRE.h"
#include "eon/LBFGS.h"
#include "eon/Quickmin.h"
#include "eon/SteepestDescent.h"
#ifdef WITH_RGMIN
#include "eon/XtsciOptimizer.h"
#endif

namespace eonc::helpers::create {

std::optional<std::string> rgminMethodFor(const ObjectiveFunction &a_objf,
                                          OptType a_otype,
                                          const Parameters &a_params) {
  if (!a_objf.gradientMatchesEnergy()) {
    return std::nullopt;
  }
  return rgminBackendMethod(a_otype, a_params);
}

std::optional<std::string> rgminBackendMethod(OptType a_otype,
                                              const Parameters &a_params) {
  if (a_params.optimizer_options().backend != "rgmin") {
    return std::nullopt;
  }
  switch (a_otype) {
  case OptType::LBFGS:
    return "lbfgs";
  case OptType::FIRE:
    return "fire";
  case OptType::CG:
    return "polak_ribiere";
  case OptType::SD:
    return "steepest";
  case OptType::QM:
    return "quickmin";
  default:
    return std::nullopt;
  }
}

std::unique_ptr<Optimizer> mkOptim(std::shared_ptr<ObjectiveFunction> a_objf,
                                   OptType a_otype,
                                   const Parameters &a_params) {
  if (auto method = rgminMethodFor(*a_objf, a_otype, a_params)) {
#ifdef WITH_RGMIN
    return std::make_unique<XtsciOptimizer>(a_objf, a_params,
                                            std::move(*method));
#else
    throw std::runtime_error(
        "[Optimizer] opt_backend=rgmin needs a build that links rgmin "
        "(-Dwith_rgmin=enabled)");
#endif
  }
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
#ifdef WITH_RGMIN
  case OptType::XTSCI: {
    return std::make_unique<XtsciOptimizer>(a_objf, a_params);
  }
#else
  case OptType::XTSCI: {
    throw std::runtime_error(
        "[Optimizer] opt_method=xtsci needs a build that links rgmin "
        "(-Dwith_rgmin=enabled)");
  }
#endif
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
} // namespace eonc::helpers::create

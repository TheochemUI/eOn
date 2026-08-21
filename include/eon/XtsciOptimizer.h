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
#pragma once

#include "eon/Optimizer.h"

namespace eonc {

/// xtsci-optimize backend. Newton / RFO live in Rust; this class is the
/// ObjectiveFunction + pair-Hessian waist.
class XtsciOptimizer final : public Optimizer {
public:
  XtsciOptimizer(std::shared_ptr<ObjectiveFunction> a_objf,
                 const Parameters &a_params)
      : Optimizer(a_objf, OptType::XTSCI,
                  OptimizerConfig::fromParams(a_params)) {}

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;
};

} // namespace eonc

using eonc::XtsciOptimizer;

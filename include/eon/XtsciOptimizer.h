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

struct xts_solver_t;

namespace eonc {

/// xtsci-optimize backend. Newton / RFO live in Rust; this class is the
/// ObjectiveFunction + pair-Hessian waist.
class XtsciOptimizer final : public Optimizer {
public:
  XtsciOptimizer(std::shared_ptr<ObjectiveFunction> a_objf,
                 const Parameters &a_params)
      : Optimizer(a_objf, OptType::XTSCI,
                  OptimizerConfig::fromParams(a_params)) {}
  ~XtsciOptimizer() override;

  XtsciOptimizer(const XtsciOptimizer &) = delete;
  XtsciOptimizer &operator=(const XtsciOptimizer &) = delete;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;

private:
  // Opaque xtsci session. L-BFGS pairs / NLCG directions live here.
  xts_solver_t *m_solver{nullptr};

  void ensureSolver(double a_maxMove);
};

} // namespace eonc

using eonc::XtsciOptimizer;

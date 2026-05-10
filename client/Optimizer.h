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

#include "Eigen.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"

#include "EonLogger.h"

namespace eonc {

/** @defgroup Optimizers
 *
 * \brief ClientEON methods for optimizing atomic structures
 *
 * This page provides links to all the available optimizers that can be run by
 * the ClientEON, as well as documentation on the optimizer class.
 *
 */

/**
 * @file
 * @ingroup Optimizers
 *
 * \brief The optimizer class is used to serve as an abstract class for all
 * optimizers, as well as to call an optimizer at runtime based off of the
 * passed in paramters.
 *
 * The set of optimizers are methods for optimizing atomic structures, solving
 * unconstrained energy minimization. Only a certain set of job types that
 * ClientEON runs can take advantage of the numeric optimizer (SEE OVERVIEW) and
 * are documented in their own files accordingly.
 *
 */

/**
 * Declaration of the optimizer class
 */

/// Configuration extracted from Parameters for the Optimizer hierarchy.
struct OptimizerConfig {
  Parameters::optimizer_options_t opts;
  double finiteDifference{0.01};
  bool bowlBreakout{false};

  static OptimizerConfig fromParams(const Parameters &p) {
    return {p.optimizer_options, p.main_options.finiteDifference,
            p.saddle_search_options.confine_positive.bowl_breakout};
  }
};

class Optimizer {
private:
  const OptType m_otype;

protected:
  const OptimizerConfig m_optConfig;
  std::shared_ptr<ObjectiveFunction> m_objf;

public:
  Optimizer(std::shared_ptr<ObjectiveFunction> a_objf,
            const OptimizerConfig &a_config)
      : m_otype{a_config.opts.method},
        m_optConfig{a_config},
        m_objf{std::move(a_objf)} {
    EONC_LOG_WARNING(
        "You should explicitly set an optimizer while constructing the "
        "optimizer!!\n Defaulting to opt_method from the parameters");
  }
  Optimizer(std::shared_ptr<ObjectiveFunction> a_objf, OptType a_optype,
            const OptimizerConfig &a_config)
      : m_otype{a_optype},
        m_optConfig{a_config},
        m_objf{std::move(a_objf)} {}

  // Backward-compat constructors
  [[deprecated("Pass OptimizerConfig directly")]]
  Optimizer(std::shared_ptr<ObjectiveFunction> a_objf,
            const Parameters &a_params)
      : Optimizer(std::move(a_objf), OptimizerConfig::fromParams(a_params)) {}
  [[deprecated("Pass OptimizerConfig directly")]]
  Optimizer(std::shared_ptr<ObjectiveFunction> a_objf, OptType a_optype,
            const Parameters &a_params)
      : Optimizer(std::move(a_objf), a_optype,
                  OptimizerConfig::fromParams(a_params)) {}

  virtual ~Optimizer() {};
  virtual int step(double a_maxMove) = 0;
  virtual int run(size_t a_maxIterations, double a_maxMove) = 0;
};

namespace helpers::create {
std::unique_ptr<Optimizer>
mkOptim(const std::shared_ptr<ObjectiveFunction> &a_objf, OptType a_otype,
        const Parameters &a_params);
}

} // namespace eonc

using eonc::Optimizer;

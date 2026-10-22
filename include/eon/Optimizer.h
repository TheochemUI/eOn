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

#include <optional>
#include <string>

namespace eonc {

/** @defgroup Optimizers Structure optimizers
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
  PotType potential{PotType::LJ};

  static OptimizerConfig fromParams(const Parameters &p) {
    return {p.optimizer_options(), p.main_options().finiteDifference,
            p.saddle_search_options().confine_positive.bowl_breakout,
            p.potential_options().potential};
  }
};

/** Abstract base of the structure optimizers. */
class Optimizer {
private:
  const OptType m_otype;

protected:
  const OptimizerConfig m_optConfig;
  std::shared_ptr<ObjectiveFunction> m_objf;

public:
  Optimizer(std::shared_ptr<ObjectiveFunction> a_objf,
            const OptimizerConfig &a_config)
      : m_otype{a_config.opts.method}, m_optConfig{a_config}, m_objf{a_objf} {
    EONC_LOG_WARNING(
        "You should explicitly set an optimizer while constructing the "
        "optimizer!!\n Defaulting to opt_method from the parameters");
  }
  Optimizer(std::shared_ptr<ObjectiveFunction> a_objf, OptType a_optype,
            const OptimizerConfig &a_config)
      : m_otype{a_optype}, m_optConfig{a_config}, m_objf{a_objf} {}

  // Backward-compat constructors
  [[deprecated("Pass OptimizerConfig directly")]]
  Optimizer(std::shared_ptr<ObjectiveFunction> a_objf,
            const Parameters &a_params)
      : Optimizer(a_objf, OptimizerConfig::fromParams(a_params)) {}
  [[deprecated("Pass OptimizerConfig directly")]]
  Optimizer(std::shared_ptr<ObjectiveFunction> a_objf, OptType a_optype,
            const Parameters &a_params)
      : Optimizer(a_objf, a_optype, OptimizerConfig::fromParams(a_params)) {}

  virtual ~Optimizer() = default;
  virtual int step(double a_maxMove) = 0;
  virtual int run(size_t a_maxIterations, double a_maxMove) = 0;
  // A min-mode search calls this after it shortens or restores a
  // translation that crossed max_energy. CG and FIRE drop the memory
  // that pointed along that step. L-BFGS stored the step's start, so
  // the next secant pair stays consistent and the default keeps it.
  virtual void onTranslationPulledBack() {}
};

namespace helpers::create {
/// The rgmin method token mkOptim runs for a_otype when [Optimizer]
/// opt_backend = rgmin. Empty for any other backend, for a method rgmin
/// does not stand in for, and for an objective whose gradient is not the
/// derivative of its energy (band and min-mode translation), which rgmin's
/// line searches would read as a potential.
[[nodiscard]] std::optional<std::string>
rgminMethodFor(const ObjectiveFunction &a_objf, OptType a_otype,
               const Parameters &a_params);

/// rgminMethodFor without the objective test: what a minimization of a
/// potential energy runs under opt_backend = rgmin.
[[nodiscard]] std::optional<std::string>
rgminBackendMethod(OptType a_otype, const Parameters &a_params);

[[nodiscard]] std::unique_ptr<Optimizer>
mkOptim(std::shared_ptr<ObjectiveFunction> a_objf, OptType a_otype,
        const Parameters &a_params);
}

} // namespace eonc

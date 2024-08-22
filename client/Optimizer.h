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

#include "BaseStructures.h"
#include "Eigen.h"
#include "ObjectiveFunction.h"
#include "client/thirdparty/toml.hpp"
#include <spdlog/spdlog.h>
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

using namespace std::string_literals;
class OptimBase {
public:
  struct Params {                                       // [Optimizer]
    OptType optM{OptType::CG};                          // opt_method
    ConvergenceMeasure optCM{ConvergenceMeasure::NORM}; // convergence_metric
    ScalarType optConvergedForce{1e-3};                 // converged_force
    size_t optMaxIter{1000L};                           // max_iterations
    ScalarType optMaxMove{0.2};                         // max_move
  };

public:
  OptimBase(const ObjectiveFunction &objectiveFunction)
      : m_objf(std::cref(objectiveFunction)) {}

  virtual ~OptimBase() = default;

  virtual bool step(ScalarType a_maxMove) = 0;
  virtual bool runOpt(size_t a_maxIter, ScalarType a_maxMove) = 0;

protected:
  const ObjectiveFunction &m_objf;
  size_t iters{0};
};

template <typename T> class Optimizer : public OptimBase {
public:
  using OptimBase::OptimBase;

protected:
  bool step(ScalarType maxMove) override final {
    return static_cast<T *>(this)->stepImpl(maxMove);
  }

  bool runOpt(size_t maxIter, ScalarType maxMove) override final {
    return static_cast<T *>(this)->runOptImpl(maxIter, maxMove);
  }
};

namespace helpers::create {
std::unique_ptr<OptimBase> mkOptim(const ObjectiveFunction &,
                                   const toml::table &);
}

} // namespace eonc

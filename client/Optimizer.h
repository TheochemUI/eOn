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
#include "client/thirdparty/toml.hpp"
#include <fmt/core.h>
#include <iostream>
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

class OptimBase {
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
    while (!m_objf.isConverged() && iters < maxIter) {
      auto pos = m_objf.getPositions();
      step(maxMove);
      // if (!static_cast<T *>(this)->stepImpl(maxMove)) {
      //   return false; // Stop if the step fails
      // }
      iters++;
      // double stepSize =
      //     helper_functions::maxAtomMotion(pbc(m_objf.getPositions() - pos));
      std::cout << fmt::format("{} {}  {} ", "[Matter]\n", iters,
                               m_objf.getConvergence(), m_objf.getEnergy());
    }
    return m_objf.isConverged();
  }
};

namespace helpers::create {
std::unique_ptr<OptimBase> mkOptim(const ObjectiveFunction &,
                                   const toml::table &);
}

} // namespace eonc

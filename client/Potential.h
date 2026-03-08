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
#include "EonLogger.h"

#include "Eigen.h"
#include "Parameters.h"
#include "PotRegistry.h"
#include <memory>

namespace eonc {

class Potential {
protected:
  PotType ptype;
  const Parameters &m_params;

private:
  uint64_t m_registry_id;
  PotRegistry::TimePoint m_created_at;

public:
  size_t forceCallCounter;

  // Main Constructor
  Potential(PotType a_ptype, const Parameters &a_params)
      : ptype{a_ptype},
        m_params{a_params},
        m_registry_id{PotRegistry::get().on_created(a_ptype)},
        m_created_at{PotRegistry::Clock::now()},
        forceCallCounter{0} {}

  // Delegating Constructor
  Potential(const Parameters &a_params)
      : Potential(a_params.potential_options.potential, a_params) {}

  virtual ~Potential() {
    PotRegistry::get().on_destroyed(m_registry_id, ptype, forceCallCounter,
                                    m_created_at);
  }

  // Does not take into account the fixed / free atoms
  // Variance here is null when not needed and that's OK
  void virtual force(long nAtoms, const double *positions, const int *atomicNrs,
                     double *forces, double *energy, double *variance,
                     const double *box) = 0;

  std::tuple<double, AtomMatrix>
  get_ef(const AtomMatrix &pos, const VectorXi &atmnrs, const Matrix3d &box);

  [[nodiscard]] PotType getType() { return this->ptype; }
};

namespace helpers {
std::shared_ptr<Potential> makePotential(const Parameters &params);
std::shared_ptr<Potential> makePotential(PotType ptype,
                                         const Parameters &params);
} // namespace helpers

} // namespace eonc

using eonc::Potential;

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
#include <algorithm>
#include <limits>
#include <memory>
#include <optional>

class Potential {
protected:
  PotType ptype;
  const Parameters &m_params;
  eonc::log::FileScoped m_log{"_potcalls", "_potcalls.log"};

public:
  size_t forceCallCounter;

  // Main Constructor
  Potential(PotType a_ptype, const Parameters &a_params)
      : ptype{a_ptype},
        m_params{a_params},
        forceCallCounter{0} {
    initializeLogger();
  }

  // Delegating Constructor
  Potential(const Parameters &a_params)
      : Potential(a_params.potential_options.potential, a_params) {}

  virtual ~Potential() {
    if (m_log) {
      LOG_TRACE_L1(m_log, "[{}] destroyed after {} calls",
                   magic_enum::enum_name<PotType>(getType()), forceCallCounter);
    } else {
      std::cerr << "Logger is not initialized\n";
    }
  }

  static int fcalls;
  static int fcallsTotal;
  static int wu_fcallsTotal;
  static double totalUserTime;

  // Does not take into account the fixed / free atoms
  // Variance here is null when not needed and that's OK
  void virtual force(long nAtoms, const double *positions, const int *atomicNrs,
                     double *forces, double *energy, double *variance,
                     const double *box) = 0;

  std::tuple<double, AtomMatrix>
  get_ef(const AtomMatrix &pos, const VectorXi &atmnrs, const Matrix3d &box);

  [[nodiscard]] PotType getType() { return this->ptype; }

  // Logger initialization
  void initializeLogger() {
    if (m_log) {
      LOG_TRACE_L1(m_log, "[{}] created",
                   magic_enum::enum_name<PotType>(getType()));
    }
  }
};

namespace helper_functions {
std::shared_ptr<Potential> makePotential(const Parameters &params);
std::shared_ptr<Potential> makePotential(PotType ptype,
                                         const Parameters &params);
} // namespace helper_functions

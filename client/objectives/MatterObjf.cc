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
#include "client/objectives/MatterObjf.hpp"
#include "client/BaseStructures.h"

namespace eonc::objf {
double MatterObjectiveFunction::getEnergy() const {
  return matter.getPotentialEnergy();
}

VectorType MatterObjectiveFunction::getGradient(bool fdstep) const {
  return -matter.getForcesFreeV();
}

void MatterObjectiveFunction::setPositions(const VectorType &x) {
  // Only mutable point
  const_cast<Matter &>(matter).setPositionsFreeV(x);
}

VectorType MatterObjectiveFunction::getPositions() const {
  return matter.getPositionsFreeV();
}

int MatterObjectiveFunction::degreesOfFreedom() const {
  return 3 * matter.numberOfFreeAtoms();
}

bool MatterObjectiveFunction::isConverged() const {
  return getConvergence() < m_p.optConvergedForce;
}

double MatterObjectiveFunction::getConvergence() const {
  switch (m_p.optCM) {
  case eonc::ConvergenceMeasure::NORM: {
    return matter.getForcesFreeV().norm();
  }
  case eonc::ConvergenceMeasure::MAX_ATOM: {
    return matter.getMaxForce();
  }
  case eonc::ConvergenceMeasure::MAX_COMPONENT: {
    return matter.getForcesFree().maxCoeff();
  }
  default: {
    SPDLOG_CRITICAL("{} Unknown opt_convergence_metric: {}", "[Matter]",
                    magic_enum::enum_name(m_p.optCM));
    throw std::runtime_error("Failure, cannot continue");
  }
  }
}

VectorType MatterObjectiveFunction::difference(const VectorType &a,
                                               const VectorType &b) const {
  return matter.pbcV(a - b);
}

} // namespace eonc::objf

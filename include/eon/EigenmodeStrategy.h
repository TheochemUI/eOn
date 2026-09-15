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

#include "ImprovedDimer.h"
#include "LowestEigenmode.h"
#include <memory>

namespace eonc {

class Matter;
class Parameters;
class Potential;

/// Type-erased eigenmode solver. Layout is the same with or without gprd.
using EigenmodeStrategy = LowestEigenmode;

std::shared_ptr<LowestEigenmode>
buildEigenmodeStrategy(std::shared_ptr<Matter> matter, const Parameters &params,
                       std::shared_ptr<Potential> pot);

inline void eigenmodeCompute(LowestEigenmode &s,
                             std::shared_ptr<Matter> matter,
                             AtomMatrix direction) {
  s.compute(matter, direction);
}

inline double eigenmodeGetEigenvalue(LowestEigenmode &s) {
  return s.getEigenvalue();
}

inline AtomMatrix eigenmodeGetEigenvector(LowestEigenmode &s) {
  return s.getEigenvector();
}

inline ImprovedDimer *asImprovedDimer(LowestEigenmode &s) {
  return dynamic_cast<ImprovedDimer *>(&s);
}

inline long eigenmodeTotalForceCalls(LowestEigenmode &s) {
  return s.totalForceCalls;
}

inline double eigenmodeStatsTorque(LowestEigenmode &s) { return s.statsTorque; }

inline double eigenmodeStatsAngle(LowestEigenmode &s) { return s.statsAngle; }

inline long eigenmodeStatsRotations(LowestEigenmode &s) {
  return s.statsRotations;
}

inline long eigenmodeTotalIterations(LowestEigenmode &s) {
  return s.totalIterations;
}

} // namespace eonc

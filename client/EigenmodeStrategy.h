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

#include "Davidson.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "LowestEigenmode.h"
#include <variant>

#ifdef WITH_GPRD
#include "AtomicGPDimer.h"
#endif

namespace eonc {

#ifdef WITH_GPRD
using EigenmodeStrategy =
    std::variant<Dimer, ImprovedDimer, Lanczos, Davidson, AtomicGPDimer>;
#else
using EigenmodeStrategy = std::variant<Dimer, ImprovedDimer, Lanczos, Davidson>;
#endif

/// Build the eigenmode solver from parameters.
/// min_mode_method: dimer (rotation CG) | lanczos | davidson | gprdimer
inline std::shared_ptr<EigenmodeStrategy>
buildEigenmodeStrategy(std::shared_ptr<Matter> matter, const Parameters &params,
                       std::shared_ptr<Potential> pot) {
  if (params.saddle_search_options.minmode_method ==
      LowestEigenmode::MINMODE_DIMER) {
    if (params.dimer_options.improved) {
      return std::make_shared<EigenmodeStrategy>(
          ImprovedDimer(matter, params, pot));
    }
    return std::make_shared<EigenmodeStrategy>(Dimer(matter, params, pot));
  } else if (params.saddle_search_options.minmode_method ==
             LowestEigenmode::MINMODE_LANCZOS) {
    return std::make_shared<EigenmodeStrategy>(Lanczos(matter, params, pot));
  } else if (params.saddle_search_options.minmode_method ==
             LowestEigenmode::MINMODE_DAVIDSON) {
    return std::make_shared<EigenmodeStrategy>(Davidson(matter, params, pot));
  }
#ifdef WITH_GPRD
  else if (params.saddle_search_options.minmode_method ==
           LowestEigenmode::MINMODE_GPRDIMER) {
    return std::make_shared<EigenmodeStrategy>(
        AtomicGPDimer(matter, params, pot));
  }
#endif
  // Default to improved dimer
  return std::make_shared<EigenmodeStrategy>(
      ImprovedDimer(matter, params, pot));
}

/// Dispatch compute() to the active variant.
inline void eigenmodeCompute(EigenmodeStrategy &s,
                             std::shared_ptr<Matter> matter,
                             AtomMatrix direction) {
  std::visit([&](auto &impl) { impl.compute(matter, direction); }, s);
}

/// Dispatch getEigenvalue() to the active variant.
inline double eigenmodeGetEigenvalue(EigenmodeStrategy &s) {
  return std::visit([](auto &impl) { return impl.getEigenvalue(); }, s);
}

/// Dispatch getEigenvector() to the active variant.
inline AtomMatrix eigenmodeGetEigenvector(EigenmodeStrategy &s) {
  return std::visit([](auto &impl) { return impl.getEigenvector(); }, s);
}

/// Access ImprovedDimer-specific features. Returns nullptr if not
/// ImprovedDimer.
inline ImprovedDimer *asImprovedDimer(EigenmodeStrategy &s) {
  return std::get_if<ImprovedDimer>(&s);
}

/// Read stats from any variant (all inherit LowestEigenmode stats fields).
inline long eigenmodeTotalForceCalls(EigenmodeStrategy &s) {
  return std::visit([](auto &impl) { return impl.totalForceCalls; }, s);
}

inline double eigenmodeStatsTorque(EigenmodeStrategy &s) {
  return std::visit([](auto &impl) { return impl.statsTorque; }, s);
}

inline double eigenmodeStatsAngle(EigenmodeStrategy &s) {
  return std::visit([](auto &impl) { return impl.statsAngle; }, s);
}

inline long eigenmodeStatsRotations(EigenmodeStrategy &s) {
  return std::visit([](auto &impl) { return impl.statsRotations; }, s);
}

inline long eigenmodeTotalIterations(EigenmodeStrategy &s) {
  return std::visit([](auto &impl) { return impl.totalIterations; }, s);
}

} // namespace eonc

using eonc::EigenmodeStrategy;

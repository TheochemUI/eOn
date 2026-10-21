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

/// Explicit Hessian-range contract for rgpot kernels.
///
/// Coloring a finite-difference Hessian is exact only when every block
/// between atoms farther apart than the reported distance is zero.
/// The primary template reports that the range is unknown. A kernel
/// opts in by specializing HessianRangeOf. The presence of
/// config().cutoff is not an opt-in: many-body and message-passing
/// kernels have a cutoff and still couple atoms beyond it.

#include <cmath>

#include "rgpot/LennardJones/LJClusterPot.hpp"
#include "rgpot/LennardJones/LJPot.hpp"
#include "rgpot/Morse/MorsePot.hpp"
#include "rgpot/ZBL/ZBLPot.hpp"
#include "rgpot/fortran/FortranPots.hpp"

#if defined(__has_include) && __has_include("rgpot/UmaPot/UmaPot.hpp")
#include "rgpot/UmaPot/UmaPot.hpp"
#endif
#if defined(__has_include) && __has_include("rgpot/Water/WaterPots.hpp")
#include "rgpot/Water/WaterPots.hpp"
#endif
#if defined(__has_include) && __has_include("rgpot/EAM/EAMPot.hpp")
#include "rgpot/EAM/EAMPot.hpp"
#endif
#if defined(__has_include) && __has_include("rgpot/EMT/EMTPot.hpp")
#include "rgpot/EMT/EMTPot.hpp"
#endif

namespace eonc {

/// How a kernel's energy couples atoms, which fixes the Hessian range.
enum class HessianLocality {
  /// No proven finite range. Coloring stays off.
  Unknown,
  /// Pair potential. Blocks vanish beyond the interaction cutoff.
  Pair,
  /// Atom-centred many-body (EAM, Gupta, bond order, three-body).
  /// Blocks vanish beyond twice the interaction cutoff.
  AtomCentred,
  /// Molecular site model. The range is the intermolecular cutoff
  /// plus twice the largest intramolecular extent, when that extent
  /// is a hard bound.
  MolecularSite,
  /// Message passing. The range is the layer count times the edge
  /// cutoff, when the layer count is part of the kernel.
  MessagePassing,
};

/// Reported Hessian range for one rgpot kernel type.
///
/// range() is the distance beyond which every Hessian block is exactly
/// zero. 0 means unknown, and the finite-difference Hessian then takes
/// one column per coordinate.
template <class RPot> struct HessianRangeOf {
  static constexpr HessianLocality locality = HessianLocality::Unknown;
  /// Interaction cutoff the range is built from. 0 when that cutoff is
  /// not a compile-time constant. Pair kernels with a runtime cutoff
  /// leave this at 0 and return the config value from range().
  static constexpr double interactionCutoff = 0.0;
  static constexpr double range(const RPot &) noexcept { return 0.0; }
};

/// Drop a non-positive or non-finite range so coloring stays off.
[[nodiscard]] inline double acceptedHessianRange(double range) noexcept {
  if (!(range > 0.0) || !std::isfinite(range)) {
    return 0.0;
  }
  return range;
}

template <> struct HessianRangeOf<rgpot::LJPot> {
  static constexpr HessianLocality locality = HessianLocality::Pair;
  /// config().cutoff is chosen per potential, so there is no compile-time
  /// cutoff. range() returns that runtime value.
  static constexpr double interactionCutoff = 0.0;
  static double range(const rgpot::LJPot &pot) noexcept {
    return pot.config().cutoff;
  }
};

template <> struct HessianRangeOf<rgpot::LJClusterPot> {
  static constexpr HessianLocality locality = HessianLocality::Pair;
  /// config().cutoff is chosen per potential, so there is no compile-time
  /// cutoff. range() returns that runtime value.
  static constexpr double interactionCutoff = 0.0;
  static double range(const rgpot::LJClusterPot &pot) noexcept {
    return pot.config().cutoff;
  }
};

template <> struct HessianRangeOf<rgpot::MorsePot> {
  static constexpr HessianLocality locality = HessianLocality::Pair;
  /// config().cutoff is chosen per potential, so there is no compile-time
  /// cutoff. range() returns that runtime value.
  static constexpr double interactionCutoff = 0.0;
  static double range(const rgpot::MorsePot &pot) noexcept {
    return pot.config().cutoff;
  }
};

/// ZBL is a pair potential. Energy and force both vanish at cut_global
/// (ZBLPot.cc), so the block range is that truncation, not cut_inner.
template <> struct HessianRangeOf<rgpot::ZBLPot> {
  static constexpr HessianLocality locality = HessianLocality::Pair;
  /// cut_global is a config field. range() returns that runtime value.
  static constexpr double interactionCutoff = 0.0;
  static double range(const rgpot::ZBLPot &pot) noexcept {
    return pot.config().cut_global;
  }
};

/// Fortran EAM aluminium. rcut = 7.1 is the default in rgpot_eam_al.f90
/// and the C++ face has no setter. Pairs and densities are zero inside
/// that radius, so a block past twice the radius is zero.
template <> struct HessianRangeOf<rgpot::fortranpots::EAMAlPot> {
  static constexpr HessianLocality locality = HessianLocality::AtomCentred;
  static constexpr double interactionCutoff = 7.1;
  static constexpr double range(const rgpot::fortranpots::EAMAlPot &) noexcept {
    return 2.0 * interactionCutoff;
  }
};

/// Stillinger-Weber. The neighbour table is built at alpha * sigma
/// (rgpot_sw.f90). The three-body term couples the two ends of a
/// triplet, so the block range is twice that cutoff.
template <> struct HessianRangeOf<rgpot::fortranpots::SWPot> {
  static constexpr HessianLocality locality = HessianLocality::AtomCentred;
  static constexpr double kSigma = 2.0951;
  static constexpr double kAlpha = 1.8;
  static constexpr double interactionCutoff = kAlpha * kSigma;
  static constexpr double range(const rgpot::fortranpots::SWPot &) noexcept {
    return 2.0 * interactionCutoff;
  }
};

/// Tersoff. The taper is zero at S = 3 Å (rgpot_tersoff.f90). The bond
/// order on a pair depends on neighbours inside S, so blocks reach 2 S.
template <> struct HessianRangeOf<rgpot::fortranpots::TersoffPot> {
  static constexpr HessianLocality locality = HessianLocality::AtomCentred;
  static constexpr double interactionCutoff = 3.0;
  static constexpr double
  range(const rgpot::fortranpots::TersoffPot &) noexcept {
    return 2.0 * interactionCutoff;
  }
};

/// EDIP. The neighbour table uses max(a, b), and a and b are both
/// 3.1213820 Å (rgpot_edip.f90). Coordination and the three-body term
/// stay inside that radius.
template <> struct HessianRangeOf<rgpot::fortranpots::EDIPPot> {
  static constexpr HessianLocality locality = HessianLocality::AtomCentred;
  static constexpr double interactionCutoff = 3.1213820;
  static constexpr double range(const rgpot::fortranpots::EDIPPot &) noexcept {
    return 2.0 * interactionCutoff;
  }
};

/// Lenosky silicon. The pair and density splines share a 4.5 Å
/// neighbour cutoff (rgpot_lenosky.f90).
template <> struct HessianRangeOf<rgpot::fortranpots::LenoskyPot> {
  static constexpr HessianLocality locality = HessianLocality::AtomCentred;
  static constexpr double interactionCutoff = 4.5;
  static constexpr double
  range(const rgpot::fortranpots::LenoskyPot &) noexcept {
    return 2.0 * interactionCutoff;
  }
};

/// Fe-He. Pair terms vanish at their own cutoffs, and the largest is the
/// He-He pair at 5.4 Å. d-band density vanishes at 4.2 Å and only Fe
/// embeds it. s-band density vanishes at 4.10 Å and both species embed
/// it. A block is zero past the larger of those pair cutoffs and twice
/// the largest density cutoff (rgpot_fehe.f90).
template <> struct HessianRangeOf<rgpot::fortranpots::FeHePot> {
  static constexpr HessianLocality locality = HessianLocality::AtomCentred;
  static constexpr double kMaxPair = 5.4;
  static constexpr double kMaxDensity = 4.2;
  /// Neighbour-table cutoff, fehe_cutoff(). The block range is wider.
  static constexpr double interactionCutoff = kMaxPair;
  static constexpr double range(const rgpot::fortranpots::FeHePot &) noexcept {
    return kMaxPair > 2.0 * kMaxDensity ? kMaxPair : 2.0 * kMaxDensity;
  }
};

/// Cu-H2 EAM. phi and rho are shifted to vanish at rcut = 6.1 Å
/// (rgpot_cuh2.f90).
template <> struct HessianRangeOf<rgpot::fortranpots::CuH2Pot> {
  static constexpr HessianLocality locality = HessianLocality::AtomCentred;
  static constexpr double interactionCutoff = 6.1;
  static constexpr double range(const rgpot::fortranpots::CuH2Pot &) noexcept {
    return 2.0 * interactionCutoff;
  }
};

/// potH-H2O. rgpot_water_h.f90 states there is no cutoff: every site
/// interacts with the lone H however far away. Coloring stays off.
template <> struct HessianRangeOf<rgpot::fortranpots::WaterHPot> {
  static constexpr HessianLocality locality = HessianLocality::Unknown;
  static constexpr double interactionCutoff = 0.0;
  static constexpr double
  range(const rgpot::fortranpots::WaterHPot &) noexcept {
    return 0.0;
  }
};

#if defined(__has_include) && __has_include("rgpot/UmaPot/UmaPot.hpp")
/// UMA builds edges at config().cutoff (UmaPot.cc) and then runs the
/// compiled graph. The layer count lives in the .pt2 package, not in
/// the kernel, so n_layers * cutoff is not derivable here.
template <> struct HessianRangeOf<rgpot::UmaPot> {
  static constexpr HessianLocality locality = HessianLocality::MessagePassing;
  static constexpr double interactionCutoff = 0.0;
  static constexpr double range(const rgpot::UmaPot &) noexcept { return 0.0; }
};
#endif

#if defined(__has_include) && __has_include("rgpot/Water/WaterPots.hpp")
/// Flexible TIP4P, SPC/E and TIP4P-Pt. WaterConfig.cutoff is the
/// intermolecular cutoff. The CCL intramolecular term
/// (ccl.hpp, intramolecular) has no distance at which it is identically
/// zero, and SPC/E also has a full Coulomb path (spce_ccl.hpp,
/// coulombFull). The extent of a molecule is not a hard bound, so
/// rc + 2 * extent is not an exact range.
template <> struct HessianRangeOf<rgpot::TIP4PPot> {
  static constexpr HessianLocality locality = HessianLocality::MolecularSite;
  static constexpr double interactionCutoff = 0.0;
  static constexpr double range(const rgpot::TIP4PPot &) noexcept {
    return 0.0;
  }
};

template <> struct HessianRangeOf<rgpot::SPCEPot> {
  static constexpr HessianLocality locality = HessianLocality::MolecularSite;
  static constexpr double interactionCutoff = 0.0;
  static constexpr double range(const rgpot::SPCEPot &) noexcept {
    return 0.0;
  }
};

template <> struct HessianRangeOf<rgpot::TIP4PPtPot> {
  static constexpr HessianLocality locality = HessianLocality::MolecularSite;
  static constexpr double interactionCutoff = 0.0;
  static constexpr double range(const rgpot::TIP4PPtPot &) noexcept {
    return 0.0;
  }
};
#endif

#if defined(__has_include) && __has_include("rgpot/EAM/EAMPot.hpp")
/// Cell-list EAM. EAMCellConfig has no cutoff. The element table
/// (EAM/Parameters.h) lists r_cut of 5.555 (Al), 5.412 (Pd) and
/// 5.5155 (Au). Blocks vanish beyond twice the largest of those.
template <> struct HessianRangeOf<rgpot::EAMPot> {
  static constexpr HessianLocality locality = HessianLocality::AtomCentred;
  static constexpr double interactionCutoff = 5.555;
  static constexpr double range(const rgpot::EAMPot &) noexcept {
    return 2.0 * interactionCutoff;
  }
};
#endif

#if defined(__has_include) && __has_include("rgpot/EMT/EMTPot.hpp")
/// Asap EMT. EMTConfig selects a parameter set and does not carry a
/// cutoff. calc_cutoff() builds one from the elements present, and the
/// weight is a Fermi function cut by the neighbour list. There is no
/// single distance that is exact for every system, so coloring stays off.
template <> struct HessianRangeOf<rgpot::EMTPot> {
  static constexpr HessianLocality locality = HessianLocality::AtomCentred;
  static constexpr double interactionCutoff = 0.0;
  static constexpr double range(const rgpot::EMTPot &) noexcept { return 0.0; }
};
#endif

} // namespace eonc

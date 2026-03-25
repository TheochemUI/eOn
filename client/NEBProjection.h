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
#include "NEBSpringForce.h"
#include "Parameters.h"
#include <variant>

namespace eonc::neb {

/// Data for a single image needed by projection strategies.
struct ImageForceData {
  const AtomMatrix &force;
  const AtomMatrix &tangent;
  const SpringResult &springResult;
  long nFreeAtoms;
  long nAtoms;
};

/// Mills, Jonsson, Schenter, Surf. Sci. 324:305, 1995.
/// Plain elastic band: full spring force + full true force (no nudging).
struct PlainEB {
  AtomMatrix project(const ImageForceData &d) const;
};

/// Jonsson, Mills, Jacobsen 1998 (World Scientific).
/// Standard NEB: parallel spring force + perpendicular true force.
struct NEB_Projection {
  AtomMatrix project(const ImageForceData &d) const;
};

/// Trygubenko & Wales, JCP 120:2082, 2004.
/// Doubly-nudged elastic band: NEB projection + perpendicular spring
/// correction.
struct DNEB_Projection {
  AtomMatrix project(const ImageForceData &d) const;
};

using ProjectionStrategy =
    std::variant<PlainEB, NEB_Projection, DNEB_Projection>;

/// Build the projection strategy from parameters.
/// DNEB is incompatible with OM and weighted springs; falls back to
/// NEB_Projection.
ProjectionStrategy buildProjectionStrategy(const Parameters &params);

/// Compute the DNEB force component for a given image.
/// Returned separately so it can be added to the CI force when DNEB is active.
AtomMatrix computeDNEBComponent(const AtomMatrix &forceSpring,
                                const AtomMatrix &tangent,
                                const AtomMatrix &forcePerp);

} // namespace eonc::neb

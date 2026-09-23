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
#include "eon/NEBProjection.h"
#include "eon/NEBForceProjection.h"

namespace eonc::neb {

AtomMatrix PlainEB::project(const ImageForceData &d) const {
  return d.springResult.forceSpring + d.force;
}

AtomMatrix NEB_Projection::project(const ImageForceData &d) const {
  AtomMatrix fPerp = eonc::neb::forcePerp(d.force, d.tangent);
  return d.springResult.forceSpringPar + fPerp;
}

AtomMatrix DNEB_Projection::project(const ImageForceData &d) const {
  AtomMatrix fPerp = eonc::neb::forcePerp(d.force, d.tangent);
  AtomMatrix forceDNEB = computeDNEBComponent(d.springResult.forceSpring,
                                              d.tangent, fPerp, use_switching);
  return d.springResult.forceSpringPar + fPerp + forceDNEB;
}

ProjectionStrategy buildProjectionStrategy(const Parameters &params) {
  bool omActive = params.neb_options().spring.om.enabled;
  bool weightedActive = params.neb_options().spring.weighting.enabled;

  if (params.neb_options().spring.use_elastic_band && !omActive &&
      !weightedActive) {
    return PlainEB{};
  }
  if (params.neb_options().spring.doubly_nudged && !omActive &&
      !weightedActive) {
    return DNEB_Projection{params.neb_options().spring.use_switching};
  }
  return NEB_Projection{};
}

AtomMatrix computeDNEBComponent(const AtomMatrix &forceSpring,
                                const AtomMatrix &tangent,
                                const AtomMatrix &fPerp, bool useSwitching) {
  return eonc::neb::computeDNEB(forceSpring, tangent, fPerp, useSwitching);
}

} // namespace eonc::neb

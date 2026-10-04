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
#include "../RpcCapabilitiesProbe.h"
#include "catch2/catch_amalgamated.hpp"

#ifndef EON_BUILD_VERSION
#define EON_BUILD_VERSION ""
#endif
#ifndef EON_BUILD_REVISION
#define EON_BUILD_REVISION ""
#endif

namespace {

void requireHandshake(const eonc::RpcCapabilitiesView &caps) {
  REQUIRE(caps.protocolFamily == "rgpot.potentials");
  REQUIRE(caps.protocolMajor == 1);
  REQUIRE(caps.protocolMinor == 0);
  REQUIRE(caps.schemaId == "0xbd1f89fa17369103");
  REQUIRE(caps.bridgeAbiMajor == 1);
  REQUIRE(caps.bridgeAbiMinor == 0);
  REQUIRE(caps.bridgeLayout == 1);
  REQUIRE(caps.dlpackMajor == 1);
  REQUIRE(caps.dlpackMinor == 0);
  REQUIRE(caps.bridgeFeatures == 0);
  REQUIRE(caps.backendName == "eon");
  REQUIRE(caps.available);
  REQUIRE(caps.buildVersion == EON_BUILD_VERSION);
  REQUIRE(caps.buildRevision == EON_BUILD_REVISION);
  REQUIRE(caps.servesEnergy);
  REQUIRE(caps.servesForces);
}

} // namespace

TEST_CASE("rpc servers implement getCapabilities", "[rpc]") {
  const auto both = eonc::readRpcServerCapabilities();
  requireHandshake(both.callbackServer);
  requireHandshake(both.pooledServer);
}

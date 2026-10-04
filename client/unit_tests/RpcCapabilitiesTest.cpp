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
#include "../RpcServerHandshake.h"
#include "catch2/catch_amalgamated.hpp"

#ifndef EON_BUILD_VERSION
#define EON_BUILD_VERSION ""
#endif
#ifndef EON_BUILD_REVISION
#define EON_BUILD_REVISION ""
#endif

namespace {

void requireHandshake(const eonc::RpcCapabilitiesView &caps) {
  REQUIRE(caps.protocolFamily == eonc::rpc::kProtocolFamily);
  REQUIRE(caps.protocolMajor == eonc::rpc::kProtocolMajor);
  REQUIRE(caps.protocolMinor == eonc::rpc::kProtocolMinor);
  REQUIRE(caps.schemaId == eonc::rpc::kSchemaId);
  REQUIRE(caps.bridgeAbiMajor == eonc::rpc::kBridgeAbiMajor);
  REQUIRE(caps.bridgeAbiMinor == eonc::rpc::kBridgeAbiMinor);
  REQUIRE(caps.bridgeLayout == eonc::rpc::kBridgeLayout);
  REQUIRE(caps.dlpackMajor == eonc::rpc::kDlpackMajor);
  REQUIRE(caps.dlpackMinor == eonc::rpc::kDlpackMinor);
  REQUIRE(caps.bridgeFeatures == eonc::rpc::kBridgeFeatures);
  REQUIRE(caps.backendName == eonc::rpc::kBackendName);
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

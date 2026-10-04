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

#include <cstdint>
#include <string>

namespace eonc {

/**
 * @brief Fields read from Potential.getCapabilities.
 *
 * The RPC test compares these with the handshake a client requires before
 * calculate. bridgeFeatures is 0 because the callback server does not
 * advertise eindir bridge feature bits.
 */
struct RpcCapabilitiesView {
  std::string protocolFamily;
  std::uint16_t protocolMajor = 0;
  std::uint16_t protocolMinor = 0;
  std::string schemaId;
  std::uint16_t bridgeAbiMajor = 0;
  std::uint16_t bridgeAbiMinor = 0;
  std::uint32_t bridgeLayout = 0;
  std::uint16_t dlpackMajor = 0;
  std::uint16_t dlpackMinor = 0;
  std::uint64_t bridgeFeatures = 0;
  std::string backendName;
  bool available = false;
  std::string buildVersion;
  std::string buildRevision;
  bool servesEnergy = false;
  bool servesForces = false;
};

struct RpcServerCapabilitiesPair {
  RpcCapabilitiesView callbackServer;
  RpcCapabilitiesView pooledServer;
};

/** Call getCapabilities on both server classes through an in-process RPC. */
RpcServerCapabilitiesPair readRpcServerCapabilities();

} // namespace eonc

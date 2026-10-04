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

namespace eonc::rpc {

/// Handshake values that Potential.getCapabilities reports. They are the
/// revisions of Potentials.capnp (the file id is its schema identity) and
/// of the eindir bridge that a client checks before calculate. The server
/// and the test read them from here.
inline constexpr const char *kProtocolFamily = "rgpot.potentials";
inline constexpr std::uint16_t kProtocolMajor = 1;
inline constexpr std::uint16_t kProtocolMinor = 0;
inline constexpr const char *kSchemaId = "0xbd1f89fa17369103";
inline constexpr std::uint16_t kBridgeAbiMajor = 1;
inline constexpr std::uint16_t kBridgeAbiMinor = 0;
inline constexpr std::uint32_t kBridgeLayout = 1;
inline constexpr std::uint16_t kDlpackMajor = 1;
inline constexpr std::uint16_t kDlpackMinor = 0;
/// This server does not advertise eindir bridge feature bits.
inline constexpr std::uint64_t kBridgeFeatures = 0;
inline constexpr const char *kBackendName = "eon";

} // namespace eonc::rpc

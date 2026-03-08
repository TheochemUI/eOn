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

/**
 * @file ServeMode.cpp
 * @brief Multi-model RPC serving for eOn potentials.
 *
 * Wraps any eOn Potential as a flat-array ForceCallback and serves it
 * over Cap'n Proto RPC. Uses ForceCallback (std::function taking flat
 * arrays) to completely avoid the AtomMatrix type collision between
 * eOn's Eigen-based AtomMatrix and rgpot's custom AtomMatrix.
 *
 * This translation unit includes eOn's Potential.h. The Cap'n Proto server
 * code is in ServeRpcServer.cpp (separate TU to avoid naming collision with
 * the capnp-generated `Potential` interface).
 */

#include "ServeMode.h"
#include "EonLogger.h"
#include "Potential.h"
#include "ServeRpcServer.h"

#include <algorithm>
#include <sstream>
#include <thread>
#include <vector>

namespace {

/// Create a ForceCallback that wraps an eOn Potential's force() method.
ForceCallback makeForceCallback(std::shared_ptr<::Potential> pot) {
  return [pot = std::move(pot)](long nAtoms, const double *positions,
                                const int *atomicNrs, double *forces,
                                double *energy, const double *box) {
    double variance = 0.0;
    pot->force(nAtoms, positions, atomicNrs, forces, energy, &variance, box);
  };
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Single-model serve
// ---------------------------------------------------------------------------

void serveMode(const Parameters &params, const std::string &host,
               uint16_t port) {
  auto pot_type = params.potential_options.potential;
  EONC_LOG_INFO("Creating potential: {}",
                std::string(magic_enum::enum_name(pot_type)));

  auto eon_pot = helper_functions::makePotential(params);
  if (!eon_pot) {
    EONC_LOG_ERROR("Failed to create potential of type {}",
                   std::string(magic_enum::enum_name(pot_type)));
    return;
  }

  auto callback = makeForceCallback(std::move(eon_pot));

  // Blocks until killed (runs Cap'n Proto event loop)
  startRpcServer(std::move(callback), host, port);
}

// ---------------------------------------------------------------------------
// Multi-model concurrent serve
// ---------------------------------------------------------------------------

void serveMultiple(const std::vector<ServeEndpoint> &endpoints,
                   const Parameters &base_params) {
  if (endpoints.empty()) {
    EONC_LOG_ERROR("No serve endpoints specified");
    return;
  }

  // Single endpoint: run in the main thread (no extra overhead)
  if (endpoints.size() == 1) {
    auto params = base_params;
    params.potential_options.potential = endpoints[0].potential;
    serveMode(params, endpoints[0].host, endpoints[0].port);
    return;
  }

  // Multiple endpoints: one thread per server
  EONC_LOG_INFO("Starting {} concurrent RPC servers", endpoints.size());

  std::vector<std::thread> threads;
  threads.reserve(endpoints.size());

  for (const auto &ep : endpoints) {
    threads.emplace_back([&base_params, ep]() {
      auto params = base_params;
      params.potential_options.potential = ep.potential;
      auto pot_name = std::string(magic_enum::enum_name(ep.potential));

      EONC_LOG_INFO("[{}:{}] Creating potential: {}", ep.host, ep.port,
                    pot_name);

      auto eon_pot = helper_functions::makePotential(params);
      if (!eon_pot) {
        EONC_LOG_ERROR("[{}:{}] Failed to create potential {}", ep.host,
                       ep.port, pot_name);
        return;
      }

      auto callback = makeForceCallback(std::move(eon_pot));
      startRpcServer(std::move(callback), ep.host, ep.port);
    });
  }

  // Wait for all threads (they block until killed)
  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }
}

// ---------------------------------------------------------------------------
// Replicated serve: N copies of same potential on sequential ports
// ---------------------------------------------------------------------------

void serveReplicated(const Parameters &params, const std::string &host,
                     uint16_t base_port, size_t replicas) {
  if (replicas == 0) {
    EONC_LOG_ERROR("Replicas must be >= 1");
    return;
  }
  if (replicas == 1) {
    serveMode(params, host, base_port);
    return;
  }

  EONC_LOG_INFO("Starting {} replicated servers on ports {}-{}", replicas,
                base_port, base_port + replicas - 1);

  std::vector<std::thread> threads;
  threads.reserve(replicas);

  for (size_t i = 0; i < replicas; ++i) {
    uint16_t port = static_cast<uint16_t>(base_port + i);
    threads.emplace_back(
        [&params, &host, port]() { serveMode(params, host, port); });
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }
}

// ---------------------------------------------------------------------------
// Gateway serve: single port backed by a pool of potential instances
// ---------------------------------------------------------------------------

void serveGateway(const Parameters &params, const std::string &host,
                  uint16_t port, size_t pool_size) {
  if (pool_size == 0) {
    EONC_LOG_ERROR("Pool size must be >= 1");
    return;
  }

  auto pot_type = params.potential_options.potential;
  EONC_LOG_INFO("Creating pool of {} {} instances for gateway on {}:{}",
                pool_size, std::string(magic_enum::enum_name(pot_type)), host,
                port);

  std::vector<ForceCallback> pool;
  pool.reserve(pool_size);

  for (size_t i = 0; i < pool_size; ++i) {
    auto eon_pot = helper_functions::makePotential(params);
    if (!eon_pot) {
      EONC_LOG_ERROR("Failed to create potential instance {}/{}", i + 1,
                     pool_size);
      return;
    }
    pool.push_back(makeForceCallback(std::move(eon_pot)));
  }

  EONC_LOG_INFO("Pool ready, starting gateway server");
  startPooledRpcServer(std::move(pool), host, port);
}

// ---------------------------------------------------------------------------
// Config-driven dispatch
// ---------------------------------------------------------------------------

void serveFromConfig(const Parameters &params) {
  const auto &opts = params.serve_options;

  // Multi-model endpoints take priority
  if (!opts.endpoints.empty()) {
    auto endpoints = parseServeSpec(opts.endpoints);
    if (endpoints.empty()) {
      EONC_LOG_ERROR("No valid endpoints in spec: {}", opts.endpoints);
      return;
    }
    serveMultiple(endpoints, params);
    return;
  }

  // Gateway mode
  if (opts.gateway_port > 0) {
    size_t pool = (opts.replicas > 0) ? opts.replicas : 1;
    serveGateway(params, opts.host, opts.gateway_port, pool);
    return;
  }

  // Replicated mode (default)
  serveReplicated(params, opts.host, opts.port, opts.replicas);
}

// ---------------------------------------------------------------------------
// Spec parser: "pot:port,pot:host:port,..."
// ---------------------------------------------------------------------------

std::vector<ServeEndpoint> parseServeSpec(const std::string &spec) {
  std::vector<ServeEndpoint> endpoints;
  std::istringstream stream(spec);
  std::string token;

  while (std::getline(stream, token, ',')) {
    // Trim whitespace
    token.erase(0, token.find_first_not_of(" \t"));
    token.erase(token.find_last_not_of(" \t") + 1);
    if (token.empty())
      continue;

    // Parse "potential:port" or "potential:host:port"
    size_t first_colon = token.find(':');
    if (first_colon == std::string::npos) {
      EONC_LOG_ERROR("Invalid serve spec '{}': expected 'potential:port'",
                     token);
      continue;
    }

    std::string pot_str = token.substr(0, first_colon);
    std::string rest = token.substr(first_colon + 1);

    // Trim parts after colon split
    pot_str.erase(0, pot_str.find_first_not_of(" \t"));
    pot_str.erase(pot_str.find_last_not_of(" \t") + 1);
    rest.erase(0, rest.find_first_not_of(" \t"));
    rest.erase(rest.find_last_not_of(" \t") + 1);

    // Lowercase the potential name
    std::transform(pot_str.begin(), pot_str.end(), pot_str.begin(), ::tolower);

    ServeEndpoint ep;
    ep.potential =
        magic_enum::enum_cast<PotType>(pot_str, magic_enum::case_insensitive)
            .value_or(PotType::UNKNOWN);

    if (ep.potential == PotType::UNKNOWN) {
      EONC_LOG_ERROR("Unknown potential type '{}'", pot_str);
      continue;
    }

    size_t second_colon = rest.find(':');
    if (second_colon != std::string::npos) {
      // "host:port" format
      ep.host = rest.substr(0, second_colon);
      ep.host.erase(0, ep.host.find_first_not_of(" \t"));
      ep.host.erase(ep.host.find_last_not_of(" \t") + 1);
      std::string port_str = rest.substr(second_colon + 1);
      port_str.erase(0, port_str.find_first_not_of(" \t"));
      port_str.erase(port_str.find_last_not_of(" \t") + 1);
      ep.port = static_cast<uint16_t>(std::stoi(port_str));
    } else {
      // "port" only
      ep.host = "localhost";
      ep.port = static_cast<uint16_t>(std::stoi(rest));
    }

    QUILL_LOG_INFO(eonc::log::get(), "Parsed endpoint: {} on {}:{}",
                   std::string(magic_enum::enum_name(ep.potential)), ep.host,
                   ep.port);
    endpoints.push_back(ep);
  }

  return endpoints;
}

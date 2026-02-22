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

#include "Parameters.h"
#include <cstdint>
#include <string>
#include <vector>

/**
 * @brief Configuration for a single serve endpoint.
 *
 * Maps a potential type to a host:port pair. Multiple endpoints can be
 * served concurrently, each in its own thread with its own potential instance.
 */
struct ServeEndpoint {
  PotType potential;
  std::string host;
  uint16_t port;
};

/**
 * @brief Start a single rgpot-compatible Cap'n Proto RPC server.
 *
 * Wraps eOn's Potential::force() as an rgpot PotentialBase and serves it
 * over Cap'n Proto RPC. Blocks until the server is killed.
 *
 * @param params  The eOn parameters (potential_options.potential must be set).
 * @param host    The hostname to listen on (e.g., "localhost" or "*").
 * @param port    The TCP port to listen on.
 */
void serveMode(const Parameters &params, const std::string &host,
               uint16_t port);

/**
 * @brief Serve multiple potentials concurrently on different ports.
 *
 * Each endpoint gets its own thread, its own potential instance, and its own
 * Cap'n Proto event loop. All threads block until SIGINT/SIGTERM.
 *
 * @param endpoints  List of {potential, host, port} configurations.
 * @param params     Base eOn parameters (potential type is overridden per
 * endpoint).
 */
void serveMultiple(const std::vector<ServeEndpoint> &endpoints,
                   const Parameters &params);

/**
 * @brief Serve N replicas of the same potential across sequential ports.
 *
 * Starts `replicas` threads, each serving the same potential type on
 * ports base_port, base_port+1, ..., base_port+replicas-1.
 *
 * @param params    The eOn parameters (potential_options.potential must be
 * set).
 * @param host      The hostname to listen on.
 * @param base_port The first port; replicas use base_port+0 .. base_port+N-1.
 * @param replicas  Number of concurrent server instances.
 */
void serveReplicated(const Parameters &params, const std::string &host,
                     uint16_t base_port, size_t replicas);

/**
 * @brief Start a gateway server backed by a pool of potential instances.
 *
 * Creates `pool_size` instances of the configured potential and serves them
 * behind a single gateway port. Incoming requests are dispatched round-robin
 * across the pool. This gives clients a single endpoint while spreading load.
 *
 * @param params     The eOn parameters (potential_options.potential must be
 * set).
 * @param host       The hostname to listen on.
 * @param port       The gateway port.
 * @param pool_size  Number of potential instances in the pool.
 */
void serveGateway(const Parameters &params, const std::string &host,
                  uint16_t port, size_t pool_size);

/**
 * @brief Start serve mode from config-file parameters.
 *
 * Reads serve_options from params and dispatches to the appropriate mode:
 * - If endpoints is set, uses multi-model serve (serveMultiple).
 * - If gateway_port > 0, uses gateway mode (serveGateway).
 * - Otherwise, uses replicated mode (serveReplicated).
 *
 * @param params  The eOn parameters with serve_options populated.
 */
void serveFromConfig(const Parameters &params);

/**
 * @brief Parse a serve configuration string into endpoints.
 *
 * Accepts comma-separated entries of the form "potential:port" or
 * "potential:host:port". Examples:
 *   "lj:12345"
 *   "metatomic:12345,lj:12346"
 *   "metatomic:0.0.0.0:12345"
 *
 * @param spec  The comma-separated endpoint specification.
 * @return A vector of ServeEndpoint structs.
 */
std::vector<ServeEndpoint> parseServeSpec(const std::string &spec);

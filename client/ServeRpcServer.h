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

#include "rgpot/Potential.hpp"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

/**
 * @brief Start a blocking Cap'n Proto RPC server for a given PotentialBase.
 *
 * This is an internal function used by serveMode(). It sets up the Cap'n Proto
 * event loop and blocks until the process is killed.
 *
 * Defined in a separate translation unit to avoid naming collision between
 * eOn's `class Potential` and the capnp-generated `Potential` interface.
 *
 * @param pot   Ownership of the rgpot PotentialBase to serve.
 * @param host  The hostname to listen on.
 * @param port  The TCP port to listen on.
 */
void startRpcServer(std::unique_ptr<rgpot::PotentialBase> pot,
                    const std::string &host, uint16_t port);

/**
 * @brief Start a blocking Cap'n Proto RPC server backed by a pool of
 *        PotentialBase instances dispatched round-robin.
 *
 * Creates a single gateway endpoint. Incoming RPC requests are dispatched
 * to the next available potential instance in a round-robin fashion. This
 * gives clients a single address while spreading computational load.
 *
 * @param pool  Vector of PotentialBase instances (ownership transferred).
 * @param host  The hostname to listen on.
 * @param port  The TCP port to listen on.
 */
void startPooledRpcServer(
    std::vector<std::unique_ptr<rgpot::PotentialBase>> pool,
    const std::string &host, uint16_t port);

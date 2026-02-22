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
 * @file ServeRpcServer.cpp
 * @brief Cap'n Proto RPC server for rgpot PotentialBase instances.
 *
 * This translation unit does NOT include eOn's Potential.h to avoid naming
 * collision with the capnp-generated `Potential` interface class.
 */

#include "ServeRpcServer.h"

#include <atomic>
#include <capnp/ez-rpc.h>
#include <kj/debug.h>
#include <mutex>
#include <spdlog/spdlog.h>

#include "rgpot/types/AtomMatrix.hpp"
#include "rgpot/types/adapters/capnp/capnp_adapter.hpp"

// Cap'n Proto generated header (from Potentials.capnp).
// This defines `class Potential` -- which collides with eOn's Potential class,
// hence the separate translation unit.
#include "Potentials.capnp.h"

namespace {

/**
 * @class GenericPotImpl
 * @brief Cap'n Proto RPC server implementation wrapping a PotentialBase.
 *
 * Same pattern as rgpot's potserv -- receives ForceInput over RPC, dispatches
 * to the polymorphic PotentialBase, returns PotentialResult.
 */
class GenericPotImpl final : public Potential::Server {
public:
  explicit GenericPotImpl(std::unique_ptr<rgpot::PotentialBase> pot)
      : m_potential(std::move(pot)) {}

  kj::Promise<void> calculate(CalculateContext context) override {
    auto fip = context.getParams().getFip();
    const size_t numAtoms = fip.getPos().size() / 3;

    KJ_REQUIRE(fip.getAtmnrs().size() == numAtoms, "AtomNumbers size mismatch");

    auto nativePositions =
        rgpot::types::adapt::capnp::convertPositionsFromCapnp(fip.getPos(),
                                                              numAtoms);
    auto nativeAtomTypes =
        rgpot::types::adapt::capnp::convertAtomNumbersFromCapnp(
            fip.getAtmnrs());
    auto nativeBoxMatrix =
        rgpot::types::adapt::capnp::convertBoxMatrixFromCapnp(fip.getBox());

    auto [energy, forces] =
        (*m_potential)(nativePositions, nativeAtomTypes, nativeBoxMatrix);

    auto result = context.getResults();
    auto pres = result.initResult();
    pres.setEnergy(energy);

    auto forcesList = pres.initForces(numAtoms * 3);
    rgpot::types::adapt::capnp::populateForcesToCapnp(forcesList, forces);

    return kj::READY_NOW;
  }

private:
  std::unique_ptr<rgpot::PotentialBase> m_potential;
};

} // anonymous namespace

void startRpcServer(std::unique_ptr<rgpot::PotentialBase> pot,
                    const std::string &host, uint16_t port) {
  spdlog::info("Starting Cap'n Proto RPC server on {}:{}", host, port);

  capnp::EzRpcServer server(kj::heap<GenericPotImpl>(std::move(pot)), host,
                             port);

  auto &waitScope = server.getWaitScope();
  spdlog::info("Server ready on port {}. Ctrl+C to stop.", port);
  kj::NEVER_DONE.wait(waitScope);
}

// ---------------------------------------------------------------------------
// Pooled (round-robin gateway) server
// ---------------------------------------------------------------------------

namespace {

/**
 * @class PooledPotImpl
 * @brief Cap'n Proto RPC server that dispatches across a pool of potentials.
 *
 * Incoming calculate() requests are assigned round-robin to the pool members.
 * Each pool member is guarded by its own mutex so concurrent RPC calls are
 * safe even though individual PotentialBase instances are not thread-safe.
 */
class PooledPotImpl final : public Potential::Server {
public:
  explicit PooledPotImpl(
      std::vector<std::unique_ptr<rgpot::PotentialBase>> pool)
      : m_pool(std::move(pool)), m_mutexes(m_pool.size()),
        m_next(0) {}

  kj::Promise<void> calculate(CalculateContext context) override {
    // Round-robin selection
    size_t idx =
        m_next.fetch_add(1, std::memory_order_relaxed) % m_pool.size();

    auto fip = context.getParams().getFip();
    const size_t numAtoms = fip.getPos().size() / 3;

    KJ_REQUIRE(fip.getAtmnrs().size() == numAtoms,
               "AtomNumbers size mismatch");

    auto nativePositions =
        rgpot::types::adapt::capnp::convertPositionsFromCapnp(fip.getPos(),
                                                              numAtoms);
    auto nativeAtomTypes =
        rgpot::types::adapt::capnp::convertAtomNumbersFromCapnp(
            fip.getAtmnrs());
    auto nativeBoxMatrix =
        rgpot::types::adapt::capnp::convertBoxMatrixFromCapnp(fip.getBox());

    // Lock the selected pool member for the duration of the force call
    std::lock_guard<std::mutex> lock(m_mutexes[idx]);
    auto [energy, forces] =
        (*m_pool[idx])(nativePositions, nativeAtomTypes, nativeBoxMatrix);

    auto result = context.getResults();
    auto pres = result.initResult();
    pres.setEnergy(energy);

    auto forcesList = pres.initForces(numAtoms * 3);
    rgpot::types::adapt::capnp::populateForcesToCapnp(forcesList, forces);

    return kj::READY_NOW;
  }

private:
  std::vector<std::unique_ptr<rgpot::PotentialBase>> m_pool;
  std::vector<std::mutex> m_mutexes;
  std::atomic<size_t> m_next;
};

} // anonymous namespace

void startPooledRpcServer(
    std::vector<std::unique_ptr<rgpot::PotentialBase>> pool,
    const std::string &host, uint16_t port) {
  spdlog::info("Starting pooled RPC gateway on {}:{} with {} instances", host,
               port, pool.size());

  capnp::EzRpcServer server(
      kj::heap<PooledPotImpl>(std::move(pool)), host, port);

  auto &waitScope = server.getWaitScope();
  spdlog::info("Gateway ready on port {}. Ctrl+C to stop.", port);
  kj::NEVER_DONE.wait(waitScope);
}

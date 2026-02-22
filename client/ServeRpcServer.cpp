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
 * @brief Cap'n Proto RPC server using flat-array force callbacks.
 *
 * This translation unit does NOT include eOn's Potential.h to avoid naming
 * collision with the capnp-generated `Potential` interface class.
 *
 * Uses a ForceCallback (std::function taking flat arrays) instead of
 * rgpot::PotentialBase to avoid the AtomMatrix type collision between
 * eOn's Eigen-based AtomMatrix and rgpot's custom AtomMatrix.
 */

#include "ServeRpcServer.h"

#include <atomic>
#include <capnp/ez-rpc.h>
#include <kj/debug.h>
#include <mutex>
#include <spdlog/spdlog.h>

// Cap'n Proto generated header (from Potentials.capnp).
// This defines `class Potential` -- which collides with eOn's Potential class,
// hence the separate translation unit.
#include "Potentials.capnp.h"

namespace {

/**
 * @class CallbackPotImpl
 * @brief Cap'n Proto RPC server wrapping a flat-array ForceCallback.
 *
 * Receives ForceInput over RPC, extracts flat arrays, calls the callback,
 * and returns PotentialResult.
 */
class CallbackPotImpl final : public Potential::Server {
public:
  explicit CallbackPotImpl(ForceCallback cb)
      : m_callback(std::move(cb)) {}

  kj::Promise<void> calculate(CalculateContext context) override {
    auto fip = context.getParams().getFip();
    const size_t numAtoms = fip.getPos().size() / 3;

    KJ_REQUIRE(fip.getAtmnrs().size() == numAtoms, "AtomNumbers size mismatch");

    // Extract flat arrays from capnp
    std::vector<double> positions(numAtoms * 3);
    auto capnpPos = fip.getPos();
    for (size_t i = 0; i < numAtoms * 3; ++i) {
      positions[i] = capnpPos[i];
    }

    std::vector<int> atomicNrs(numAtoms);
    auto capnpAtmnrs = fip.getAtmnrs();
    for (size_t i = 0; i < numAtoms; ++i) {
      atomicNrs[i] = capnpAtmnrs[i];
    }

    double box[9] = {};
    auto capnpBox = fip.getBox();
    for (size_t i = 0; i < 9 && i < capnpBox.size(); ++i) {
      box[i] = capnpBox[i];
    }

    // Call the force callback
    std::vector<double> forces(numAtoms * 3, 0.0);
    double energy = 0.0;
    m_callback(static_cast<long>(numAtoms), positions.data(), atomicNrs.data(),
               forces.data(), &energy, box);

    // Serialize result back to capnp
    auto result = context.getResults();
    auto pres = result.initResult();
    pres.setEnergy(energy);

    auto forcesList = pres.initForces(numAtoms * 3);
    for (size_t i = 0; i < numAtoms * 3; ++i) {
      forcesList.set(i, forces[i]);
    }

    return kj::READY_NOW;
  }

private:
  ForceCallback m_callback;
};

} // anonymous namespace

void startRpcServer(ForceCallback callback, const std::string &host,
                    uint16_t port) {
  spdlog::info("Starting Cap'n Proto RPC server on {}:{}", host, port);

  capnp::EzRpcServer server(kj::heap<CallbackPotImpl>(std::move(callback)),
                            host, port);

  auto &waitScope = server.getWaitScope();
  spdlog::info("Server ready on port {}. Ctrl+C to stop.", port);
  kj::NEVER_DONE.wait(waitScope);
}

// ---------------------------------------------------------------------------
// Pooled (round-robin gateway) server
// ---------------------------------------------------------------------------

namespace {

/**
 * @class PooledCallbackPotImpl
 * @brief Cap'n Proto RPC server dispatching across a pool of callbacks.
 */
class PooledCallbackPotImpl final : public Potential::Server {
public:
  explicit PooledCallbackPotImpl(std::vector<ForceCallback> pool)
      : m_pool(std::move(pool)),
        m_mutexes(m_pool.size()),
        m_next(0) {}

  kj::Promise<void> calculate(CalculateContext context) override {
    size_t idx = m_next.fetch_add(1, std::memory_order_relaxed) % m_pool.size();

    auto fip = context.getParams().getFip();
    const size_t numAtoms = fip.getPos().size() / 3;

    KJ_REQUIRE(fip.getAtmnrs().size() == numAtoms, "AtomNumbers size mismatch");

    std::vector<double> positions(numAtoms * 3);
    auto capnpPos = fip.getPos();
    for (size_t i = 0; i < numAtoms * 3; ++i) {
      positions[i] = capnpPos[i];
    }

    std::vector<int> atomicNrs(numAtoms);
    auto capnpAtmnrs = fip.getAtmnrs();
    for (size_t i = 0; i < numAtoms; ++i) {
      atomicNrs[i] = capnpAtmnrs[i];
    }

    double box[9] = {};
    auto capnpBox = fip.getBox();
    for (size_t i = 0; i < 9 && i < capnpBox.size(); ++i) {
      box[i] = capnpBox[i];
    }

    std::vector<double> forces(numAtoms * 3, 0.0);
    double energy = 0.0;

    std::lock_guard<std::mutex> lock(m_mutexes[idx]);
    m_pool[idx](static_cast<long>(numAtoms), positions.data(), atomicNrs.data(),
                forces.data(), &energy, box);

    auto result = context.getResults();
    auto pres = result.initResult();
    pres.setEnergy(energy);

    auto forcesList = pres.initForces(numAtoms * 3);
    for (size_t i = 0; i < numAtoms * 3; ++i) {
      forcesList.set(i, forces[i]);
    }

    return kj::READY_NOW;
  }

private:
  std::vector<ForceCallback> m_pool;
  std::vector<std::mutex> m_mutexes;
  std::atomic<size_t> m_next;
};

} // anonymous namespace

void startPooledRpcServer(std::vector<ForceCallback> pool,
                          const std::string &host, uint16_t port) {
  spdlog::info("Starting pooled RPC gateway on {}:{} with {} instances", host,
               port, pool.size());

  capnp::EzRpcServer server(kj::heap<PooledCallbackPotImpl>(std::move(pool)),
                            host, port);

  auto &waitScope = server.getWaitScope();
  spdlog::info("Gateway ready on port {}. Ctrl+C to stop.", port);
  kj::NEVER_DONE.wait(waitScope);
}

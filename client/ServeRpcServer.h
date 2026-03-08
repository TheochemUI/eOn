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
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace eonc {



/**
 * @brief Callback type for potential energy/force evaluation.
 *
 * Flat-array interface that avoids any AtomMatrix type dependency,
 * preventing collisions between eOn's Eigen-based AtomMatrix and
 * rgpot's custom AtomMatrix.
 *
 * @param nAtoms    Number of atoms.
 * @param positions Flat array [x1,y1,z1, x2,y2,z2, ...] (nAtoms*3).
 * @param atomicNrs Atomic numbers [Z1, Z2, ...] (nAtoms).
 * @param forces    Output forces [Fx1,Fy1,Fz1, ...] (nAtoms*3).
 * @param energy    Output energy (single double).
 * @param box       Simulation cell flat 3x3 row-major (9 doubles).
 */
using ForceCallback = std::function<void(long nAtoms, const double *positions,
                                         const int *atomicNrs, double *forces,
                                         double *energy, const double *box)>;

/**
 * @brief Start a blocking Cap'n Proto RPC server using a force callback.
 *
 * Defined in a separate translation unit to avoid naming collision between
 * eOn's `class Potential` and the capnp-generated `Potential` interface.
 *
 * @param callback  Force evaluation function.
 * @param host      The hostname to listen on.
 * @param port      The TCP port to listen on.
 */
void startRpcServer(ForceCallback callback, const std::string &host,
                    uint16_t port);

/**
 * @brief Start a blocking Cap'n Proto RPC server backed by a pool of
 *        force callbacks dispatched round-robin.
 *
 * @param pool  Vector of force callbacks (one per pool instance).
 * @param host  The hostname to listen on.
 * @param port  The TCP port to listen on.
 */
void startPooledRpcServer(std::vector<ForceCallback> pool,
                          const std::string &host, uint16_t port);

} // namespace eonc

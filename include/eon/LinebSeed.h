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

#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"

#include <memory>
#include <string>
#include <vector>

namespace eonc::lineb {

/// Angstrom tolerance for cell entries and frozen-substrate positions.
inline constexpr double kSeedMatchTol = 1e-8;

/// True when ref and image share atom count, periodic flags, cell,
/// atomic numbers, fixed masks, and (for fully fixed atoms) positions
/// under the minimum image of ref's cell.
bool matchingSeedImage(const Matter &ref, const Matter &image);

/// Refuse a mismatched seed image (identity, cell, masks, frozen substrate).
void requireMatchingSeedImage(const Matter &ref, const Matter &image);

/// Refuse mismatched or non-positive atomic masses between seed images.
void requireMatchingSeedMasses(const Matter &ref, const Matter &image);

/// Free atoms that are free on every axis. Partial per-axis freezes throw.
void collectFreeAtoms(const Matter &m, std::vector<long> &freeAtoms,
                      std::vector<double> &masses);

/// Pack free-atom coordinates into a flat [image][atom][xyz] buffer for the
/// lineb C API. Adjacent images are joined by the minimum-image step of the
/// reference cell when periodic; otherwise the stored Cartesian values are
/// used. Unwrapped coordinates are never written through Matter::setPositions
/// (which rewraps under PBC). Molecular seeds under PBC are refused as
/// ambiguous without a bond topology; Rust performs mass-weighted alignment.
void packSeedCartesian(const std::vector<std::shared_ptr<Matter>> &band,
                       const std::vector<long> &freeAtoms, bool molecular,
                       std::vector<double> &out);

/// Load a multi-frame seed and enforce identity, masses, masks, frozen
/// substrate, and finite positions. Does not evaluate potentials: endpoints
/// may sit at or below E; interior feasible width is checked by rgsaddle.
std::vector<std::shared_ptr<Matter>>
loadSeedBand(const Parameters &params, const std::shared_ptr<Potential> &pot,
             const std::string &path);

} // namespace eonc::lineb

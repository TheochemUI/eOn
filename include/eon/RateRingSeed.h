/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#pragma once

#include "eon/ConFileIO.h"
#include "eon/Matter.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc::tunneling {

/// Read closed-ring geometry on the reference's ordered physical coordinates.
/// Cached energies and forces do not participate in the seed.
inline std::vector<AtomMatrix> readRateRingSeed(const std::string &filename,
                                                const Matter &reference,
                                                long beads) {
  const auto frames = readcon::read_all_frames(filename);
  if (beads < 4 || frames.size() != static_cast<size_t>(beads)) {
    throw std::invalid_argument(
        "instanton: initial_ring must contain exactly beads frames, without "
        "a repeated closing frame");
  }
  std::vector<AtomMatrix> positions;
  positions.reserve(frames.size());
  for (const auto &frame : frames) {
    Matter image(reference);
    if (!io::io_ok(io::con2matter(image, frame))) {
      throw std::runtime_error("instanton: cannot read ring seed " + filename);
    }
    auto mismatch = [&]() {
      throw std::invalid_argument(
          "instanton: ring seed differs in ordered species, IDs, masses, "
          "cell, constraints or frozen coordinates");
    };
    if (image.numberOfAtoms() != reference.numberOfAtoms() ||
        image.getPeriodic() != reference.getPeriodic() ||
        !image.getCell().allFinite() || !reference.getCell().allFinite() ||
        (image.getCell() - reference.getCell()).cwiseAbs().maxCoeff() > 1e-8) {
      mismatch();
    }
    const AtomMatrix delta =
        reference.pbc(image.getPositions() - reference.getPositions());
    if (!delta.allFinite()) {
      mismatch();
    }
    for (long atom = 0; atom < reference.numberOfAtoms(); ++atom) {
      if (image.getAtomicNr(atom) != reference.getAtomicNr(atom) ||
          image.getAtomIndex(atom) != reference.getAtomIndex(atom) ||
          image.getFixedMask(atom) != reference.getFixedMask(atom) ||
          image.getMass(atom) != reference.getMass(atom) ||
          !std::isfinite(image.getMass(atom)) || !(image.getMass(atom) > 0.) ||
          (reference.getFixed(atom) &&
           delta.row(atom).cwiseAbs().maxCoeff() > 1e-8)) {
        mismatch();
      }
    }
    positions.push_back(image.getPositions());
  }
  return positions;
}

} // namespace eonc::tunneling

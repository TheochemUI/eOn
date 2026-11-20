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
#include "eon/LinebSeed.h"

#include "eon/ConFileIO.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace eonc::lineb {

bool matchingSeedImage(const Matter &ref, const Matter &image) {
  if (ref.numberOfAtoms() != image.numberOfAtoms() ||
      ref.getPeriodic() != image.getPeriodic() || !ref.getCell().allFinite() ||
      !image.getCell().allFinite() ||
      (ref.getCell() - image.getCell()).cwiseAbs().maxCoeff() > kSeedMatchTol ||
      !(ref.getAtomicNrs().array() == image.getAtomicNrs().array()).all()) {
    return false;
  }
  const AtomMatrix delta = ref.pbc(image.getPositions() - ref.getPositions());
  if (!delta.allFinite()) {
    return false;
  }
  for (long i = 0; i < ref.numberOfAtoms(); ++i) {
    if (ref.getFixedMask(i) != image.getFixedMask(i)) {
      return false;
    }
    if (ref.getFixed(i) && delta.row(i).cwiseAbs().maxCoeff() > kSeedMatchTol) {
      return false;
    }
  }
  return true;
}

void requireMatchingSeedImage(const Matter &ref, const Matter &image) {
  if (!matchingSeedImage(ref, image)) {
    throw std::runtime_error(
        "lineb: seed images differ in cell, species, periodic flags, frozen "
        "masks or frozen-substrate positions");
  }
}

void requireMatchingSeedMasses(const Matter &ref, const Matter &image) {
  if (ref.numberOfAtoms() != image.numberOfAtoms()) {
    throw std::runtime_error("lineb: seed images differ in atom count");
  }
  for (long i = 0; i < ref.numberOfAtoms(); ++i) {
    const double mRef = ref.getMass(i);
    const double mImg = image.getMass(i);
    if (!(mRef > 0.0) || !std::isfinite(mRef) || !(mImg > 0.0) ||
        !std::isfinite(mImg) || mRef != mImg) {
      throw std::runtime_error(
          "lineb: seed images differ in atomic masses (every atom needs a "
          "finite positive mass shared by every image)");
    }
  }
}

void collectFreeAtoms(const Matter &m, std::vector<long> &freeAtoms,
                      std::vector<double> &masses) {
  freeAtoms.clear();
  masses.clear();
  for (long i = 0; i < m.numberOfAtoms(); ++i) {
    const int fx = m.getFixed(i, 0);
    const int fy = m.getFixed(i, 1);
    const int fz = m.getFixed(i, 2);
    if (fx && fy && fz) {
      continue;
    }
    if (fx || fy || fz) {
      throw std::runtime_error(
          "lineb: per-axis frozen coordinates are unsupported; fix or free "
          "all three axes of each atom");
    }
    const double mass = m.getMass(i);
    if (!(mass > 0.0) || !std::isfinite(mass)) {
      throw std::runtime_error(
          "lineb: every free atom needs a finite positive mass");
    }
    freeAtoms.push_back(i);
    masses.push_back(mass);
  }
  if (freeAtoms.empty()) {
    throw std::runtime_error("lineb: every atom is fixed");
  }
}

void packSeedCartesian(const std::vector<std::shared_ptr<Matter>> &band,
                       const std::vector<long> &freeAtoms, bool molecular,
                       std::vector<double> &out) {
  if (band.empty()) {
    throw std::runtime_error("lineb: seed band is empty");
  }
  if (freeAtoms.empty()) {
    throw std::runtime_error("lineb: every atom is fixed");
  }
  const Matter &ref = *band.front();
  if (molecular) {
    if (ref.getPeriodic()) {
      throw std::runtime_error(
          "lineb: molecular seed under periodic boundaries is ambiguous "
          "without a bond topology; use a non-periodic cell");
    }
    if (static_cast<long>(freeAtoms.size()) != ref.numberOfAtoms()) {
      throw std::runtime_error(
          "lineb: molecular alignment needs every atom free (no frozen "
          "substrate atoms)");
    }
    if (freeAtoms.size() < 2) {
      throw std::runtime_error(
          "lineb: molecular alignment needs at least two free atoms");
    }
  }

  const size_t nImages = band.size();
  const size_t nFree = freeAtoms.size();
  out.assign(nImages * nFree * 3, 0.0);

  auto writeImage = [&](size_t image, const AtomMatrix &pos) {
    for (size_t a = 0; a < nFree; ++a) {
      const long atom = freeAtoms[a];
      const size_t base = (image * nFree + a) * 3;
      out[base + 0] = pos(atom, 0);
      out[base + 1] = pos(atom, 1);
      out[base + 2] = pos(atom, 2);
    }
  };

  // Keep a local continuous chain. Never call setPositions: Matter rewraps
  // under PBC and would erase the adjacent unwrap before packing.
  AtomMatrix prev = band.front()->getPositions();
  if (!prev.allFinite()) {
    throw std::runtime_error("lineb: seed image has non-finite positions");
  }
  writeImage(0, prev);

  for (size_t i = 1; i < nImages; ++i) {
    const AtomMatrix raw = band[i]->getPositions();
    if (!raw.allFinite()) {
      throw std::runtime_error("lineb: seed image has non-finite positions");
    }
    requireMatchingSeedImage(ref, *band[i]);
    const AtomMatrix step = ref.pbc(raw - prev);
    if (!step.allFinite()) {
      throw std::runtime_error(
          "lineb: seed path has a non-finite minimum-image displacement");
    }
    AtomMatrix next = prev + step;
    for (long a = 0; a < ref.numberOfAtoms(); ++a) {
      if (ref.getFixed(a)) {
        next.row(a) = ref.getPositions().row(a);
      }
    }
    writeImage(i, next);
    prev = std::move(next);
  }
}

std::vector<std::shared_ptr<Matter>>
loadSeedBand(const Parameters &params, const std::shared_ptr<Potential> &pot,
             const std::string &path) {
  if (path.empty()) {
    throw std::runtime_error(
        "lineb: [Lineb] initial_path must name a multi-frame seed band");
  }
  const auto frames = readcon::read_all_frames(path);
  if (frames.size() < 3) {
    throw std::runtime_error(
        "lineb: the seed band needs at least three images (endpoints and one "
        "interior image)");
  }
  std::vector<std::shared_ptr<Matter>> band;
  band.reserve(frames.size());
  for (const auto &frame : frames) {
    auto image = std::make_shared<Matter>(pot, params);
    if (!io::io_ok(io::con2matter(*image, frame))) {
      throw std::runtime_error("lineb: cannot read a frame of " + path);
    }
    if (!image->getPositions().allFinite()) {
      throw std::runtime_error("lineb: seed image has non-finite positions");
    }
    band.push_back(std::move(image));
  }

  const Matter &ref = *band.front();
  for (size_t i = 1; i < band.size(); ++i) {
    requireMatchingSeedImage(ref, *band[i]);
    requireMatchingSeedMasses(ref, *band[i]);
  }
  return band;
}

} // namespace eonc::lineb

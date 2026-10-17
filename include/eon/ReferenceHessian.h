/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#pragma once

#include "eon/ConFileIO.h"
#include "eon/Matter.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace eonc::tunneling {

/// Read a mass-weighted reference matrix bound to filename + ".con".
/// The caller verifies external-potential and finite-difference provenance.
inline MatrixXd readReferenceHessian(const std::string &filename,
                                     const Matter &reference, long dimension) {
  if (dimension <= 0 || reference.numberOfAtoms() <= 0) {
    throw std::invalid_argument("instanton: empty reference Hessian coordinates");
  }
  Matter recorded(reference);
  if (!io::io_ok(recorded.con2matter(filename + ".con"))) {
    throw std::runtime_error("instanton: cannot read Hessian geometry " + filename + ".con");
  }
  auto mismatch = [&]() {
    throw std::invalid_argument("instanton: Hessian reference geometry or mass metric differs: " + filename);
  };
  if (recorded.numberOfAtoms() != reference.numberOfAtoms() ||
      recorded.getPeriodic() != reference.getPeriodic() ||
      !recorded.getCell().allFinite() || !reference.getCell().allFinite() ||
      (recorded.getCell() - reference.getCell()).cwiseAbs().maxCoeff() > 1e-8) {
    mismatch();
  }
  const AtomMatrix delta = reference.pbc(recorded.getPositions() - reference.getPositions());
  if (!delta.allFinite() || delta.cwiseAbs().maxCoeff() > 1e-8) {
    mismatch();
  }
  long coordinates = 0;
  for (long i = 0; i < reference.numberOfAtoms(); ++i) {
    if (recorded.getAtomicNr(i) != reference.getAtomicNr(i) ||
        recorded.getAtomIndex(i) != reference.getAtomIndex(i) ||
        recorded.getFixedMask(i) != reference.getFixedMask(i) ||
        recorded.getMass(i) != reference.getMass(i) ||
        !(recorded.getMass(i) > 0.0) || !std::isfinite(recorded.getMass(i))) {
      mismatch();
    }
    if (!reference.getFixed(i)) {
      coordinates += 3;
    }
  }
  if (dimension <= 0 || dimension != coordinates) {
    throw std::invalid_argument("instanton: Hessian coordinate dimension differs: " + filename);
  }
  std::ifstream stream(filename);
  MatrixXd matrix(dimension, dimension);
  for (long i = 0; i < dimension; ++i) {
    for (long j = 0; j < dimension; ++j) {
      if (!(stream >> matrix(i, j)) || !std::isfinite(matrix(i, j))) {
        throw std::runtime_error("instanton: incomplete or non-finite reference Hessian: " + filename);
      }
    }
  }
  std::string extra;
  if (stream >> extra) {
    throw std::runtime_error("instanton: extra reference Hessian data: " + filename);
  }
  const double scale = std::max(1.0, matrix.cwiseAbs().maxCoeff());
  if ((matrix - matrix.transpose()).cwiseAbs().maxCoeff() >
      128.0 * std::numeric_limits<double>::epsilon() * scale) {
    throw std::invalid_argument("instanton: reference Hessian is asymmetric: " + filename);
  }
  return (0.5 * (matrix + matrix.transpose())).eval();
}

} // namespace eonc::tunneling

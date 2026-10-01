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
#include "eon/SolidStateNEB.h"

#include <cmath>
#include <stdexcept>

namespace eonc::neb {
namespace {

void zeroUpper(Matrix3d &cell) {
  cell(0, 1) = 0.0;
  cell(0, 2) = 0.0;
  cell(1, 2) = 0.0;
}

double wrapFrac(double value) {
  double wrapped = value - std::floor(value);
  if (wrapped >= 0.5) {
    wrapped -= 1.0;
  }
  return wrapped;
}

AtomMatrix fractionalCoordinates(const Matter &image) {
  return image.getPositions() * image.getCell().inverse();
}

} // namespace

double solidStateJacobian(double meanVolume, int nAtoms, double weight) {
  if (!(meanVolume > 0.0) || nAtoms < 1 || !(weight > 0.0)) {
    throw std::invalid_argument(
        "solid-state Jacobian needs a positive volume, atom count, and weight");
  }
  const double n = static_cast<double>(nAtoms);
  return std::cbrt(meanVolume / n) * std::sqrt(n) * weight;
}

bool orientCellLowerTriangular(Matrix3d &cell, AtomMatrix &positions) {
  const Vector3d a = cell.row(0);
  const Vector3d b = cell.row(1);
  const double na = a.norm();
  if (!(na > 1e-12)) {
    return false;
  }
  const Vector3d e1 = a / na;
  const Vector3d bPerp = b - b.dot(e1) * e1;
  const double nb = bPerp.norm();
  if (!(nb > 1e-12)) {
    return false;
  }
  const Vector3d e2 = bPerp / nb;
  const Vector3d e3 = e1.cross(e2);
  Matrix3d rotation;
  rotation.col(0) = e1;
  rotation.col(1) = e2;
  rotation.col(2) = e3;
  Matrix3d oriented = cell * rotation;
  zeroUpper(oriented);
  if (positions.rows() > 0) {
    positions = positions * rotation;
  }
  if (!(std::abs(oriented.determinant()) > 1e-18)) {
    return false;
  }
  cell = oriented;
  return true;
}

void orientSolidStateMatter(Matter &image) {
  Matrix3d cell = image.getCell();
  AtomMatrix positions = image.getPositions();
  if (!orientCellLowerTriangular(cell, positions)) {
    throw std::invalid_argument(
        "solid_state NEB could not put the cell in lower-triangular form");
  }
  image.setCell(cell);
  if (positions.rows() > 0) {
    image.setPositions(positions);
  }
}

void interpolateSolidStateLinear(std::vector<Matter> &images) {
  const auto n = static_cast<long>(images.size());
  if (n < 3) {
    throw std::invalid_argument(
        "solid_state interpolation needs two endpoints and one image");
  }
  const Matrix3d h0 = images.front().getCell();
  const Matrix3d h1 = images.back().getCell();
  const AtomMatrix s0 = fractionalCoordinates(images.front());
  AtomMatrix ds = fractionalCoordinates(images.back()) - s0;
  for (int atom = 0; atom < ds.rows(); ++atom) {
    for (int axis = 0; axis < 3; ++axis) {
      ds(atom, axis) = wrapFrac(ds(atom, axis));
    }
  }
  const double denom = static_cast<double>(n - 1);
  for (long i = 1; i < n - 1; ++i) {
    const double t = static_cast<double>(i) / denom;
    Matrix3d h = (1.0 - t) * h0 + t * h1;
    zeroUpper(h);
    const AtomMatrix fractional = s0 + t * ds;
    images[static_cast<size_t>(i)].setCell(h);
    images[static_cast<size_t>(i)].setPositions(fractional * h);
  }
}

JointBlock jointDisplacement(const Matter &from, const Matter &to,
                             double jacobian) {
  if (!(jacobian > 0.0)) {
    throw std::invalid_argument(
        "solid-state displacement needs a positive Jacobian");
  }
  const Matrix3d hFrom = from.getCell();
  const Matrix3d hTo = to.getCell();
  AtomMatrix frac = fractionalCoordinates(to) - fractionalCoordinates(from);
  for (int atom = 0; atom < frac.rows(); ++atom) {
    for (int axis = 0; axis < 3; ++axis) {
      frac(atom, axis) = wrapFrac(frac(atom, axis));
    }
  }
  const Matrix3d average = 0.5 * (hFrom + hTo);
  JointBlock out;
  out.atomic = frac * average;
  const Matrix3d dh = jacobian * (hTo - hFrom);
  out.cell = 0.5 * (hFrom.inverse() * dh + hTo.inverse() * dh);
  zeroUpper(out.cell);
  return out;
}

double jointNorm(const JointBlock &block) {
  return std::hypot(block.atomic.norm(), block.cell.norm());
}

Matrix3d cellNebForce(const Matrix3d &cauchy, double volume, double jacobian,
                      const Matrix3d &externalStress) {
  if (!(volume > 0.0) || !(jacobian > 0.0)) {
    throw std::invalid_argument(
        "solid-state cell force needs a positive volume and Jacobian");
  }
  Matrix3d force = -(volume / jacobian) * (cauchy + externalStress);
  zeroUpper(force);
  return force;
}

Matrix3d finiteDifferenceCauchyStress(const Matter &image, double strainStep) {
  if (!(strainStep > 0.0)) {
    throw std::invalid_argument(
        "stress finite difference needs a positive step");
  }
  const Matrix3d cell = image.getCell();
  const AtomMatrix positions = image.getPositions();
  const double volume = std::abs(cell.determinant());
  if (!(volume > 0.0)) {
    throw std::invalid_argument(
        "stress finite difference needs a nonzero cell");
  }
  Matrix3d sigma = Matrix3d::Zero();
  const int rows[6] = {0, 1, 1, 2, 2, 2};
  const int cols[6] = {0, 0, 1, 0, 1, 2};
  for (int comp = 0; comp < 6; ++comp) {
    Matrix3d strain = Matrix3d::Zero();
    strain(rows[comp], cols[comp]) = strainStep;
    const Matrix3d plus = Matrix3d::Identity() + strain;
    const Matrix3d minus = Matrix3d::Identity() - strain;
    Matter raised(image);
    raised.setCell(cell * plus);
    raised.setPositions(positions * plus);
    Matter lowered(image);
    lowered.setCell(cell * minus);
    lowered.setPositions(positions * minus);
    const double derivative =
        (raised.getPotentialEnergy() - lowered.getPotentialEnergy()) /
        (2.0 * strainStep);
    sigma(rows[comp], cols[comp]) = derivative / volume;
  }
  return sigma;
}

double solidStateEnthalpy(const Matter &image, const Matter &reference,
                          double pressure) {
  const double energy = image.getPotentialEnergy();
  if (pressure == 0.0) {
    return energy;
  }
  const Matrix3d h0 = reference.getCell();
  const Matrix3d strain = h0.inverse() * (image.getCell() - h0);
  const double volume = std::abs(h0.determinant());
  return energy + pressure * strain.trace() * volume;
}

CartesianStep solidStateCartesianStep(const Matter &image,
                                      const AtomMatrix &atomicForce,
                                      const Matrix3d &cellForce,
                                      double jacobian) {
  if (!(jacobian > 0.0)) {
    throw std::invalid_argument("solid-state step needs a positive Jacobian");
  }
  if (atomicForce.rows() != image.numberOfAtoms()) {
    throw std::invalid_argument("solid-state step force row count mismatch");
  }
  const Matrix3d cell = image.getCell();
  Matrix3d strain = cellForce / jacobian;
  zeroUpper(strain);
  CartesianStep step;
  step.cell = cell * strain;
  zeroUpper(step.cell);
  step.positions = atomicForce + image.getPositions() * strain;
  for (long atom = 0; atom < image.numberOfAtoms(); ++atom) {
    if (image.getFixed(atom)) {
      step.positions.row(atom).setZero();
    }
  }
  return step;
}

} // namespace eonc::neb

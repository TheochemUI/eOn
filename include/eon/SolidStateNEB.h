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

#include "Eigen.h"
#include "Matter.h"

#include <vector>

namespace eonc::neb {

/// J = (V/N)^{1/3} * N^{1/2} * weight, with V the mean endpoint volume.
/// doi:10.1063/1.3684549
double solidStateJacobian(double meanVolume, int nAtoms, double weight);

/// Rotate the lab frame so the cell is lower triangular: the first lattice
/// vector lies on x and the second lies in the xy plane. Fractional
/// coordinates and the sign of the volume are unchanged. Returns false when
/// the cell is singular or two lattice vectors are parallel.
bool orientCellLowerTriangular(Matrix3d &cell, AtomMatrix &positions);

void orientSolidStateMatter(Matter &image);

/// Replace interior images by a fractional linear interpolation of the
/// oriented endpoints. `images` includes the endpoints.
void interpolateSolidStateLinear(std::vector<Matter> &images);

struct JointBlock {
  AtomMatrix atomic;
  Matrix3d cell{Matrix3d::Zero()};
};

/// Displacement of `to` relative to `from` in the joint metric.
/// Atomic rows are fractional minimum-image displacements mapped through the
/// average cell. The cell block is the averaged right strain of the
/// Jacobian-scaled cell difference.
JointBlock jointDisplacement(const Matter &from, const Matter &to,
                             double jacobian);

double jointNorm(const JointBlock &block);

/// NEB force on the Jacobian-scaled strain.
/// Cauchy stress uses sigma = (1/V) dE/dε for h <- h (I+ε) at fixed
/// fractional coordinates. External stress is added before the -V factor
/// (positive hydrostatic pressure pushes the cell inward).
Matrix3d cellNebForce(const Matrix3d &cauchy, double volume, double jacobian,
                      const Matrix3d &externalStress);

/// Central difference of the potential energy on the six lower strain
/// components. Upper-triangle components stay zero.
Matrix3d finiteDifferenceCauchyStress(const Matter &image, double strainStep);

/// Potential energy plus P : (h0^{-1} (h-h0)) * V0. Pressure is hydrostatic,
/// in eV/Angstrom^3. A zero pressure returns the potential energy.
double solidStateEnthalpy(const Matter &image, const Matter &reference,
                          double pressure);

/// One steepest step in Cartesian coordinates and the cell matrix.
/// delta_h = h (F_cell / J), delta_r = F_atomic + r (F_cell / J).
/// Fixed atoms keep delta_r = 0. The upper triangle of delta_h is zero.
struct CartesianStep {
  AtomMatrix positions;
  Matrix3d cell{Matrix3d::Zero()};
};

CartesianStep solidStateCartesianStep(const Matter &image,
                                      const AtomMatrix &atomicForce,
                                      const Matrix3d &cellForce,
                                      double jacobian);

} // namespace eonc::neb

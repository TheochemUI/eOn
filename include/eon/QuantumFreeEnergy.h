#pragma once

#include "Eigen.h"
#include "Matter.h"
#include "Parameters.h"

#include <memory>
#include <vector>

namespace eonc {

/// kT ln(2 sinh(beta hbar omega / 2)) in eV. A non-positive temperature
/// is the zero-temperature limit, hbar omega / 2.
double harmonicModeFreeEnergy(double hbarOmegaEv, double temperature);

/// Eigenvalues of the mass-weighted Hessian after the tangent direction is
/// removed. The eigenvalue nearest zero, the projected tangent, is dropped.
/// `tangent` is mass-weighted and has one entry per Hessian row.
Eigen::VectorXd perpendicularEigenvalues(const Eigen::MatrixXd &hessian,
                                         const Eigen::VectorXd &tangent);

/// Sum of harmonicModeFreeEnergy over the positive perpendicular eigenvalues.
double perpendicularHarmonicFreeEnergy(const Eigen::MatrixXd &hessian,
                                       const Eigen::VectorXd &tangent,
                                       double temperature);

/// F(s) = V(s) + the perpendicular harmonic free energy, one entry per image.
/// The tangent is the Cartesian NEB tangent, mass-weighted inside.
std::vector<double> quantumFreeEnergies(
    const std::vector<std::shared_ptr<Matter>> &path,
    const std::vector<std::shared_ptr<AtomMatrix>> &tangent, double temperature,
    const Parameters &params);

} // namespace eonc

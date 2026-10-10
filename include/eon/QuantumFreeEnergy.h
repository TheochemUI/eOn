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

/// Eigenvalues of the mass-weighted Hessian on the orthogonal complement of
/// the columns of `removed`, any spanning set of the directions to leave
/// out: a zero column, or one in the span of the others, removes nothing
/// more. No column leaves the whole spectrum.
Eigen::VectorXd complementEigenvalues(const Eigen::MatrixXd &hessian,
                                      const Eigen::MatrixXd &removed);

/// Sum of harmonicModeFreeEnergy over the positive eigenvalues on the
/// complement of `removed`.
double complementHarmonicFreeEnergy(const Eigen::MatrixXd &hessian,
                                    const Eigen::MatrixXd &removed,
                                    double temperature);

/// Eigenvalues of the mass-weighted Hessian after the tangent direction is
/// removed: one fewer than the Hessian has rows. `tangent` is mass-weighted
/// and has one entry per Hessian row.
Eigen::VectorXd perpendicularEigenvalues(const Eigen::MatrixXd &hessian,
                                         const Eigen::VectorXd &tangent);

/// Sum of harmonicModeFreeEnergy over the positive perpendicular eigenvalues.
double perpendicularHarmonicFreeEnergy(const Eigen::MatrixXd &hessian,
                                       const Eigen::VectorXd &tangent,
                                       double temperature);

/// F(s) = V(s) + the harmonic free energy of the vibrations, one entry per
/// image. With no atom fixed, the translations, and each rotation the first
/// image's Hessian leaves null (a free cluster, periodic cell or not), are
/// rigid motions and carry no vibration. The end images are minima, where
/// every vibration counts, the one along the band included; between them
/// the band tangent is the reaction coordinate and leaves the sum. The
/// tangent is the Cartesian NEB tangent, mass-weighted inside.
std::vector<double> quantumFreeEnergies(
    const std::vector<std::shared_ptr<Matter>> &path,
    const std::vector<std::shared_ptr<AtomMatrix>> &tangent, double temperature,
    const Parameters &params);

} // namespace eonc

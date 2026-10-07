#include "eon/QuantumFreeEnergy.h"

#include "eon/Hessian.h"
#include "eon/Tunneling.h"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <stdexcept>

namespace eonc {
namespace {

double lnTwoSinh(double x) {
  if (!(x > 0.0) || !std::isfinite(x)) {
    throw std::invalid_argument("the sinh argument must be positive");
  }
  if (x > 16.0) {
    return x + std::log1p(-std::exp(-2.0 * x));
  }
  return std::log(2.0 * std::sinh(x));
}

Eigen::VectorXi mobileAtoms(const Matter &matter) {
  const long n = matter.numberOfAtoms();
  Eigen::VectorXi ids(n);
  long count = 0;
  for (long atom = 0; atom < n; ++atom) {
    if (!matter.getFixed(atom)) {
      ids(count++) = atom;
    }
  }
  ids.conservativeResize(count);
  return ids;
}

Eigen::VectorXd massWeightedTangent(const Matter &matter,
                                    const Eigen::VectorXi &mobile,
                                    const AtomMatrix &cartesian) {
  if (cartesian.rows() != matter.numberOfAtoms() || cartesian.cols() != 3) {
    throw std::invalid_argument("the tangent is not N by 3");
  }
  Eigen::VectorXd tangent(3 * mobile.size());
  for (long k = 0; k < mobile.size(); ++k) {
    const long atom = mobile(k);
    const double mass = matter.getMass(atom);
    if (!(mass > 0.0)) {
      throw std::invalid_argument("a perpendicular mode needs a positive mass");
    }
    const double scale = std::sqrt(mass);
    for (int dim = 0; dim < 3; ++dim) {
      tangent(3 * k + dim) = cartesian(atom, dim) * scale;
    }
  }
  return tangent;
}

AtomMatrix imageTangent(const std::vector<std::shared_ptr<Matter>> &path,
                        const std::vector<std::shared_ptr<AtomMatrix>> &tangent,
                        long index) {
  const long last = static_cast<long>(path.size()) - 1;
  if (index <= 0) {
    return path[0]->pbc(path[1]->getPositions() - path[0]->getPositions());
  }
  if (index >= last) {
    return path[static_cast<size_t>(last - 1)]->pbc(
        path[static_cast<size_t>(last)]->getPositions() -
        path[static_cast<size_t>(last - 1)]->getPositions());
  }
  if (static_cast<size_t>(index) < tangent.size() && tangent[index]) {
    return *tangent[index];
  }
  return path[static_cast<size_t>(index - 1)]->pbc(
      path[static_cast<size_t>(index + 1)]->getPositions() -
      path[static_cast<size_t>(index - 1)]->getPositions());
}

} // namespace

double harmonicModeFreeEnergy(double hbarOmegaEv, double temperature) {
  if (!(hbarOmegaEv > 0.0) || !std::isfinite(hbarOmegaEv)) {
    throw std::invalid_argument(
        "a perpendicular mode needs a positive hbar omega");
  }
  if (!(temperature > 0.0)) {
    return 0.5 * hbarOmegaEv;
  }
  const double x =
      hbarOmegaEv / (2.0 * tunneling::kBoltzmann * temperature);
  return tunneling::kBoltzmann * temperature * lnTwoSinh(x);
}

Eigen::VectorXd perpendicularEigenvalues(const Eigen::MatrixXd &hessian,
                                         const Eigen::VectorXd &tangent) {
  const long n = hessian.rows();
  if (hessian.cols() != n || tangent.size() != n || n < 1) {
    throw std::invalid_argument("the Hessian and the tangent sizes differ");
  }
  const double norm = tangent.norm();
  if (!(norm > 0.0) || !std::isfinite(norm)) {
    throw std::invalid_argument("the tangent is zero");
  }
  const Eigen::VectorXd unit = tangent / norm;
  const Eigen::MatrixXd projector =
      Eigen::MatrixXd::Identity(n, n) - unit * unit.transpose();
  const Eigen::MatrixXd projected = projector * hessian * projector;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(projected);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("the perpendicular Hessian did not diagonalize");
  }
  const Eigen::VectorXd values = solver.eigenvalues();
  long drop = 0;
  double nearest = std::abs(values(0));
  for (long i = 1; i < values.size(); ++i) {
    const double gap = std::abs(values(i));
    if (gap < nearest) {
      nearest = gap;
      drop = i;
    }
  }
  Eigen::VectorXd kept(values.size() - 1);
  long write = 0;
  for (long i = 0; i < values.size(); ++i) {
    if (i == drop) {
      continue;
    }
    kept(write++) = values(i);
  }
  return kept;
}

double perpendicularHarmonicFreeEnergy(const Eigen::MatrixXd &hessian,
                                       const Eigen::VectorXd &tangent,
                                       double temperature) {
  const Eigen::VectorXd lambda = perpendicularEigenvalues(hessian, tangent);
  double sum = 0.0;
  for (long i = 0; i < lambda.size(); ++i) {
    if (!(lambda(i) > 0.0)) {
      continue;
    }
    const double hw = tunneling::kHbar * std::sqrt(lambda(i));
    sum += harmonicModeFreeEnergy(hw, temperature);
  }
  return sum;
}

std::vector<double> quantumFreeEnergies(
    const std::vector<std::shared_ptr<Matter>> &path,
    const std::vector<std::shared_ptr<AtomMatrix>> &tangent, double temperature,
    const Parameters &params) {
  std::vector<double> freeEnergy;
  freeEnergy.reserve(path.size());
  for (long image = 0; image < static_cast<long>(path.size()); ++image) {
    Matter &matter = *path[static_cast<size_t>(image)];
    const double potential = matter.getPotentialEnergy();
    const Eigen::VectorXi mobile = mobileAtoms(matter);
    if (mobile.size() == 0) {
      freeEnergy.push_back(potential);
      continue;
    }
    Hessian hess(params, &matter);
    hess.writeHessianFile(false);
    const Eigen::MatrixXd weighted = hess.getHessian(&matter, mobile);
    if (weighted.rows() != 3 * mobile.size()) {
      throw std::runtime_error("the image Hessian is empty");
    }
    const AtomMatrix cart = imageTangent(path, tangent, image);
    const Eigen::VectorXd direction =
        massWeightedTangent(matter, mobile, cart);
    freeEnergy.push_back(
        potential +
        perpendicularHarmonicFreeEnergy(weighted, direction, temperature));
  }
  return freeEnergy;
}

} // namespace eonc

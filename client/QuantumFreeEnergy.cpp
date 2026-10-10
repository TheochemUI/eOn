#include "eon/QuantumFreeEnergy.h"

#include "eon/Hessian.h"
#include "eon/Tunneling.h"

#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <array>
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

// Three translations, then three rotations about the centre of mass, over
// the mobile atoms in mass-weighted coordinates. Positions are taken at
// their minimum image from atom 0, so a cluster across a cell face turns as
// one.
Eigen::MatrixXd rigidGenerators(const Matter &matter,
                                const Eigen::VectorXi &mobile) {
  const AtomMatrix raw = matter.getPositions();
  AtomMatrix diff = raw;
  diff.rowwise() -= raw.row(0);
  AtomMatrix whole = matter.pbc(diff);
  whole.rowwise() += raw.row(0);
  Eigen::RowVector3d com = Eigen::RowVector3d::Zero();
  double total = 0.0;
  for (long k = 0; k < mobile.size(); ++k) {
    const double mass = matter.getMass(mobile(k));
    com += mass * whole.row(mobile(k));
    total += mass;
  }
  com /= total;
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(3 * mobile.size(), 6);
  for (long k = 0; k < mobile.size(); ++k) {
    const double scale = std::sqrt(matter.getMass(mobile(k)));
    const Eigen::Vector3d x = (whole.row(mobile(k)) - com).transpose();
    for (int c = 0; c < 3; ++c) {
      out(3 * k + c, c) = scale;
      Eigen::Vector3d axis = Eigen::Vector3d::Zero();
      axis(c) = 1.0;
      out.block(3 * k, 3 + c, 3, 1) = scale * axis.cross(x);
    }
  }
  return out;
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
  const double x = hbarOmegaEv / (2.0 * tunneling::kBoltzmann * temperature);
  return tunneling::kBoltzmann * temperature * lnTwoSinh(x);
}

Eigen::VectorXd complementEigenvalues(const Eigen::MatrixXd &hessian,
                                      const Eigen::MatrixXd &removed) {
  const long n = hessian.rows();
  if (hessian.cols() != n || n < 1 ||
      (removed.cols() > 0 && removed.rows() != n)) {
    throw std::invalid_argument(
        "the Hessian and the removed directions differ in size");
  }
  const Eigen::MatrixXd sym = 0.5 * (hessian + hessian.transpose());
  // Unit columns, so that the rank decision compares directions and not
  // their lengths. A column at round-off of the longest (the rotation about
  // a linear molecule's axis) removes nothing, rather than a direction of
  // noise scaled up.
  Eigen::MatrixXd unit(n, removed.cols());
  long kept = 0;
  const double largest =
      removed.cols() > 0 ? removed.colwise().norm().maxCoeff() : 0.0;
  if (!std::isfinite(largest)) {
    throw std::invalid_argument("a removed direction is not finite");
  }
  for (long c = 0; c < removed.cols(); ++c) {
    const double norm = removed.col(c).norm();
    if (norm > 1e-8 * largest) {
      unit.col(kept++) = removed.col(c) / norm;
    }
  }
  Eigen::MatrixXd complement = Eigen::MatrixXd::Identity(n, n);
  if (kept > 0) {
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(unit.leftCols(kept));
    qr.setThreshold(1e-8);
    const long rank = qr.rank();
    const Eigen::MatrixXd q = qr.householderQ();
    complement = q.rightCols(n - rank);
  }
  if (complement.cols() == 0) {
    return Eigen::VectorXd();
  }
  const Eigen::MatrixXd reduced = complement.transpose() * sym * complement;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(
      0.5 * (reduced + reduced.transpose()), Eigen::EigenvaluesOnly);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("the projected Hessian did not diagonalize");
  }
  return solver.eigenvalues();
}

double complementHarmonicFreeEnergy(const Eigen::MatrixXd &hessian,
                                    const Eigen::MatrixXd &removed,
                                    double temperature) {
  const Eigen::VectorXd lambda = complementEigenvalues(hessian, removed);
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
  return complementEigenvalues(hessian, tangent / norm);
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

std::vector<double>
quantumFreeEnergies(const std::vector<std::shared_ptr<Matter>> &path,
                    const std::vector<std::shared_ptr<AtomMatrix>> &tangent,
                    double temperature, const Parameters &params) {
  std::vector<double> freeEnergy;
  freeEnergy.reserve(path.size());
  const long last = static_cast<long>(path.size()) - 1;
  std::array<bool, 3> rotations{{false, false, false}};
  bool rotationsSet = false;
  for (long image = 0; image <= last; ++image) {
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
    std::vector<Eigen::VectorXd> removed;
    if (mobile.size() == matter.numberOfAtoms()) {
      const Eigen::MatrixXd generators = rigidGenerators(matter, mobile);
      if (!rotationsSet) {
        rotations = tunneling::rotationZeroModes(weighted, generators).zero;
        rotationsSet = true;
      }
      for (int c = 0; c < 3; ++c) {
        removed.emplace_back(generators.col(c));
        if (rotations[static_cast<size_t>(c)]) {
          removed.emplace_back(generators.col(3 + c));
        }
      }
    }
    if (image > 0 && image < last) {
      removed.push_back(massWeightedTangent(
          matter, mobile, imageTangent(path, tangent, image)));
    }
    Eigen::MatrixXd columns(weighted.rows(), static_cast<long>(removed.size()));
    for (size_t k = 0; k < removed.size(); ++k) {
      columns.col(static_cast<long>(k)) = removed[k];
    }
    freeEnergy.push_back(potential + complementHarmonicFreeEnergy(
                                         weighted, columns, temperature));
  }
  return freeEnergy;
}

} // namespace eonc

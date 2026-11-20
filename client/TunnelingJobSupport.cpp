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
#include "eon/TunnelingJobSupport.h"

#include "eon/Hessian.h"
#include "eon/Potential.h"

#include <Eigen/QR>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <thread>

namespace eonc {
namespace {

void requireWholeAtomConstraints(const Matter &matter) {
  for (long atom = 0; atom < matter.numberOfAtoms(); ++atom) {
    bool anyFixed = false;
    bool anyFree = false;
    for (int axis = 0; axis < 3; ++axis) {
      anyFixed = anyFixed || matter.getFixed(atom, axis);
      anyFree = anyFree || !matter.getFixed(atom, axis);
    }
    if (anyFixed && anyFree) {
      throw std::invalid_argument(
          "tunnelling: partial Cartesian constraints are unsupported");
    }
  }
}

} // namespace

MassWeighted::MassWeighted(const Matter &reference)
    : reference_(reference) {
  for (long i = 0; i < reference.numberOfAtoms(); ++i) {
    if (reference.getFixed(i)) {
      continue;
    }
    const double mass = reference.getMass(i);
    if (!(mass > 0.0)) {
      throw std::invalid_argument(
          "tunnelling: a free atom does not have a positive mass");
    }
    free_.push_back(i);
    sqrtMass_.push_back(std::sqrt(mass));
  }
  if (free_.empty()) {
    throw std::invalid_argument("tunnelling: every atom is fixed");
  }
}

long MassWeighted::dimension() const {
  return 3 * static_cast<long>(free_.size());
}

VectorXi MassWeighted::freeAtoms() const {
  VectorXi out(static_cast<long>(free_.size()));
  for (size_t k = 0; k < free_.size(); ++k) {
    out(static_cast<long>(k)) = static_cast<int>(free_[k]);
  }
  return out;
}

std::vector<double>
MassWeighted::spreadAbout(const std::vector<VectorXd> &beads,
                          VectorXd &centroid, long atoms) const {
  if (beads.empty()) {
    throw std::invalid_argument(
        "tunnelling: cannot compute a spread without replicas");
  }
  centroid = VectorXd::Zero(dimension());
  for (const auto &q : beads) {
    centroid += q;
  }
  centroid /= static_cast<double>(beads.size());
  std::vector<double> out(static_cast<size_t>(3 * atoms), 0.0);
  for (size_t k = 0; k < free_.size(); ++k) {
    for (int c = 0; c < 3; ++c) {
      const long i = static_cast<long>(3 * k) + c;
      double sumSquares = 0.0;
      for (const auto &q : beads) {
        const double displacement = q(i) - centroid(i);
        sumSquares += displacement * displacement;
      }
      out[static_cast<size_t>(3 * free_[k] + c)] =
          std::sqrt(sumSquares / static_cast<double>(beads.size())) /
          sqrtMass_[k];
    }
  }
  return out;
}

VectorXd MassWeighted::toQ(const Matter &matter) const {
  const AtomMatrix displacement =
      reference_.pbc(matter.getPositions() - reference_.getPositions());
  VectorXd q(dimension());
  for (size_t k = 0; k < free_.size(); ++k) {
    for (int c = 0; c < 3; ++c) {
      q(static_cast<long>(3 * k) + c) =
          sqrtMass_[k] * displacement(free_[k], c);
    }
  }
  return q;
}

void MassWeighted::place(const VectorXd &q, Matter &matter) const {
  if (q.size() != dimension()) {
    throw std::invalid_argument(
        "tunnelling: mass-weighted coordinate has the wrong dimension");
  }
  AtomMatrix positions = reference_.getPositions();
  for (size_t k = 0; k < free_.size(); ++k) {
    for (int c = 0; c < 3; ++c) {
      positions(free_[k], c) += q(static_cast<long>(3 * k) + c) / sqrtMass_[k];
    }
  }
  matter.setPositions(positions);
}

MatrixXd MassWeighted::rigidGenerators(const Matter &matter) const {
  const long n = dimension();
  MatrixXd generators = MatrixXd::Zero(n, 6);
  const AtomMatrix positions = matter.getPositions();
  Eigen::RowVector3d centre = Eigen::RowVector3d::Zero();
  double totalMass = 0.0;
  for (size_t k = 0; k < free_.size(); ++k) {
    const double mass = sqrtMass_[k] * sqrtMass_[k];
    centre += mass * positions.row(free_[k]);
    totalMass += mass;
  }
  centre /= totalMass;
  for (size_t k = 0; k < free_.size(); ++k) {
    const long i = static_cast<long>(3 * k);
    const Eigen::Vector3d x = (positions.row(free_[k]) - centre).transpose();
    for (int c = 0; c < 3; ++c) {
      generators(i + c, c) = sqrtMass_[k];
      Eigen::Vector3d axis = Eigen::Vector3d::Zero();
      axis(c) = 1.0;
      generators.block(i, 3 + c, 3, 1) = sqrtMass_[k] * axis.cross(x);
    }
  }
  return generators;
}

void MassWeighted::markRotationZeroModes(
    const MatrixXd &hessian, const MatrixXd &generators,
    std::array<bool, 3> &keep, std::array<double, 3> &residual) const {
  const tunneling::RotationZeroModes modes =
      tunneling::rotationZeroModes(hessian, generators);
  keep = modes.zero;
  residual = modes.residual;
}

MatrixXd MassWeighted::rigidBasis(const Matter &matter,
                                  const std::array<bool, 3> &rotations) const {
  if (static_cast<long>(free_.size()) != matter.numberOfAtoms()) {
    return MatrixXd(dimension(), 0);
  }
  const MatrixXd generators = rigidGenerators(matter);
  const long n = generators.rows();
  std::vector<int> columns{0, 1, 2};
  for (int c = 0; c < 3; ++c) {
    if (rotations[static_cast<size_t>(c)]) {
      columns.push_back(3 + c);
    }
  }
  MatrixXd basis(n, static_cast<long>(columns.size()));
  for (size_t k = 0; k < columns.size(); ++k) {
    basis.col(static_cast<long>(k)) = generators.col(columns[k]);
  }
  const Eigen::ColPivHouseholderQR<MatrixXd> qr(basis);
  const long rank = qr.rank();
  if (rank <= 0) {
    return MatrixXd(n, 0);
  }
  return qr.householderQ() * MatrixXd::Identity(n, rank);
}

const std::vector<double> &MassWeighted::sqrtMasses() const {
  return sqrtMass_;
}

double MassWeighted::minimumSqrtMass() const {
  return *std::min_element(sqrtMass_.begin(), sqrtMass_.end());
}

VectorXd MassWeighted::referenceFree() const {
  const AtomMatrix positions = reference_.getPositions();
  VectorXd out(dimension());
  for (size_t k = 0; k < free_.size(); ++k) {
    for (int c = 0; c < 3; ++c) {
      out(static_cast<long>(3 * k) + c) = positions(free_[k], c);
    }
  }
  return out;
}

VectorXd MassWeighted::gradient(const AtomMatrix &forces) const {
  VectorXd out(dimension());
  for (size_t k = 0; k < free_.size(); ++k) {
    for (int c = 0; c < 3; ++c) {
      out(static_cast<long>(3 * k) + c) = -forces(free_[k], c) / sqrtMass_[k];
    }
  }
  return out;
}

void alignRigid(const Matter &reference, Matter &matter) {
  const long atoms = reference.numberOfAtoms();
  for (long i = 0; i < atoms; ++i) {
    if (reference.getFixed(i)) {
      return;
    }
  }
  AtomMatrix displacement =
      reference.pbc(matter.getPositions() - reference.getPositions());
  VectorXd masses(atoms);
  for (long i = 0; i < atoms; ++i) {
    masses(i) = reference.getMass(i);
  }
  const double totalMass = masses.sum();
  const Eigen::RowVector3d shift =
      (masses.transpose() * displacement) / totalMass;
  displacement.rowwise() -= shift;
  if (!reference.getPeriodic()) {
    const Eigen::RowVector3d centre =
        (masses.transpose() * reference.getPositions()) / totalMass;
    AtomMatrix referenceCentred = reference.getPositions();
    referenceCentred.rowwise() -= centre;
    const AtomMatrix movingCentred = referenceCentred + displacement;
    const Eigen::Matrix3d covariance =
        movingCentred.transpose() * masses.asDiagonal() * referenceCentred;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(covariance, Eigen::ComputeFullU |
                                                          Eigen::ComputeFullV);
    Eigen::Matrix3d determinantFix = Eigen::Matrix3d::Identity();
    determinantFix(2, 2) =
        (svd.matrixV() * svd.matrixU().transpose()).determinant() < 0.0 ? -1.0
                                                                        : 1.0;
    const Eigen::Matrix3d rotation =
        svd.matrixV() * determinantFix * svd.matrixU().transpose();
    displacement = movingCentred * rotation.transpose() - referenceCentred;
  }
  matter.setPositions(reference.getPositions() + displacement);
}

TunnelingJobContext::TunnelingJobContext(const Parameters &parameters,
                                         std::shared_ptr<Potential> potential,
                                         const Matter &reference)
    : TunnelingJobContext(parameters, std::move(potential), reference,
                          Options{}) {}

TunnelingJobContext::TunnelingJobContext(const Parameters &parameters,
                                         std::shared_ptr<Potential> potential,
                                         const Matter &reference,
                                         Options options)
    : parameters_(parameters),
      potential_(std::move(potential)),
      reference_(reference),
      coordinates_(reference_),
      forwardNonPeriodicCell_(options.forwardNonPeriodicCell),
      hessianProductThreads_(options.hessianProductThreads) {
  if (hessianProductThreads_ < 1) {
    throw std::invalid_argument(
        "tunnelling: Hessian-product thread count must be positive");
  }
  if (options.requireWholeAtomConstraints) {
    requireWholeAtomConstraints(reference);
  }
}

const MassWeighted &TunnelingJobContext::coordinates() const {
  return coordinates_;
}

void TunnelingJobContext::evaluate(const std::vector<VectorXd> &q,
                                   std::vector<double> &energies,
                                   std::vector<VectorXd> &gradients) {
  while (pool_.size() < q.size()) {
    pool_.push_back(std::make_unique<Matter>(reference_));
  }
  std::vector<Matter *> images;
  images.reserve(q.size());
  for (size_t j = 0; j < q.size(); ++j) {
    coordinates_.place(q[j], *pool_[j]);
    images.push_back(pool_[j].get());
  }
  if (forwardNonPeriodicCell_) {
    evaluateTogether(*potential_, images);
  } else if (potential_->supportsBatchEvaluation() && images.size() > 1) {
    const long atoms = reference_.numberOfAtoms();
    std::vector<VectorXi> atomicNumbers(images.size());
    std::vector<Matrix3d> boxes(images.size());
    std::vector<const double *> positionPointers;
    std::vector<const int *> numberPointers;
    std::vector<double *> forcePointers;
    std::vector<const double *> boxPointers;
    positionPointers.reserve(images.size());
    numberPointers.reserve(images.size());
    forcePointers.reserve(images.size());
    boxPointers.reserve(images.size());
    for (size_t j = 0; j < images.size(); ++j) {
      atomicNumbers[j] = images[j]->getAtomicNrs();
      boxes[j] = images[j]->getPeriodic() ? images[j]->getCell()
                                          : Matrix3d::Zero().eval();
      positionPointers.push_back(images[j]->getPositions().data());
      numberPointers.push_back(atomicNumbers[j].data());
      forcePointers.push_back(images[j]->forcesData());
      boxPointers.push_back(boxes[j].data());
    }
    std::vector<double> batchEnergies(images.size());
    std::vector<double> variances(images.size());
    potential_->forceBatch(static_cast<long>(images.size()), atoms,
                           positionPointers.data(), numberPointers.data(),
                           forcePointers.data(), batchEnergies.data(),
                           variances.data(), boxPointers.data());
    for (size_t j = 0; j < images.size(); ++j) {
      images[j]->setComputedPotential(batchEnergies[j], variances[j]);
    }
  } else {
    for (Matter *image : images) {
      image->getForcesRaw();
    }
  }
  energies.resize(q.size());
  gradients.resize(q.size());
  for (size_t j = 0; j < q.size(); ++j) {
    energies[j] = pool_[j]->getPotentialEnergy();
    gradients[j] = coordinates_.gradient(pool_[j]->getForces());
    if (!std::isfinite(energies[j]) || !gradients[j].allFinite()) {
      throw std::runtime_error(
          "tunnelling potential returned a non-finite energy or gradient");
    }
  }
}

MatrixXd TunnelingJobContext::hessian(const VectorXd &q,
                                      bool projectRotations) {
  Matter matter(reference_);
  coordinates_.place(q, matter);
  MatrixXd out;
  if (potential_->supportsMassWeightedHessianVectorProduct()) {
    const long freedom = coordinates_.dimension();
    out.resize(freedom, freedom);
    const VectorXi freeAtoms = coordinates_.freeAtoms();
    for (long column = 0; column < freedom; ++column) {
      const VectorXd product = potential_->massWeightedHessianVectorProduct(
          matter, VectorXd::Unit(freedom, column), freeAtoms);
      if (product.size() != freedom || !product.allFinite()) {
        throw std::runtime_error(
            "tunnelling: potential returned an invalid Hessian-vector "
            "product");
      }
      out.col(column) = product;
      products_->fetch_add(1, std::memory_order_relaxed);
    }
    out = 0.5 * (out + out.transpose()).eval();
  } else {
    reportFiniteDifferenceFallback();
    Hessian finiteDifference(parameters_, &matter);
    finiteDifference.writeHessianFile(false);
    out = finiteDifference.getHessian(&matter, coordinates_.freeAtoms());
  }
  if (out.rows() != coordinates_.dimension() ||
      out.cols() != coordinates_.dimension() || !out.allFinite()) {
    throw std::runtime_error("tunnelling: a replica Hessian failed");
  }
  if (!rotationsKnown_) {
    coordinates_.markRotationZeroModes(out,
                                       coordinates_.rigidGenerators(matter),
                                       rotationZero_, rotationResidual_);
    rotationsKnown_ = true;
  }
  const std::array<bool, 3> translationsOnly{{false, false, false}};
  const MatrixXd rigid = coordinates_.rigidBasis(
      matter, projectRotations ? rotationZero_ : translationsOnly);
  if (rigid.cols() > 0) {
    const MatrixXd projector =
        MatrixXd::Identity(out.rows(), out.cols()) - rigid * rigid.transpose();
    out = projector * out * projector;
  }
  return out;
}

MatrixXd TunnelingJobContext::hessianAt(const VectorXd &q) {
  return hessian(q, true);
}

MatrixXd TunnelingJobContext::beadHessianAt(const VectorXd &q) {
  if (!rotationsKnown_) {
    hessian(VectorXd::Zero(coordinates_.dimension()), true);
  }
  return hessian(q, false);
}

void TunnelingJobContext::reportFiniteDifferenceFallback() {
  if (!reportedFiniteDifferenceFallback_) {
    EONC_LOG_WARNING(
        "[Tunnelling] potential has no analytic mass-weighted "
        "Hessian-vector product; using finite-difference Hessians with "
        "finite_difference={:.8g}",
        parameters_.main_options().finiteDifference);
    reportedFiniteDifferenceFallback_ = true;
  }
}

tunneling::HessianProduct
TunnelingJobContext::hessianProductAt(const VectorXd &q) {
  if (!potential_->supportsMassWeightedHessianVectorProduct()) {
    auto matrix = std::make_shared<MatrixXd>(hessian(q, false));
    return [matrix = std::move(matrix), count = products_](
               const VectorXd &direction) {
      if (direction.size() != matrix->cols()) {
        throw std::invalid_argument(
            "tunnelling: Hessian direction has the wrong size");
      }
      count->fetch_add(1, std::memory_order_relaxed);
      return (*matrix) * direction;
    };
  }
  auto matter = std::make_shared<Matter>(reference_);
  coordinates_.place(q, *matter);
  return [matter, potential = potential_, freeAtoms = coordinates_.freeAtoms(),
          freedom = coordinates_.dimension(),
          count = products_](const VectorXd &direction) -> VectorXd {
    VectorXd out = potential->massWeightedHessianVectorProduct(
        *matter, direction, freeAtoms);
    count->fetch_add(1, std::memory_order_relaxed);
    if (out.size() != freedom || !out.allFinite()) {
      throw std::runtime_error(
          "tunnelling: potential returned an invalid Hessian-vector product");
    }
    return out;
  };
}

tunneling::RingHessianProducts
TunnelingJobContext::ringHessianProducts(const std::vector<VectorXd> &beads) {
  if (!potential_->supportsMassWeightedHessianVectorProduct()) {
    std::vector<MatrixXd> matrices;
    matrices.reserve(beads.size());
    for (const VectorXd &q : beads) {
      matrices.push_back(hessian(q, false));
    }
    return [matrices = std::move(matrices), count = products_](
               const std::vector<VectorXd> &v, std::vector<VectorXd> &out) {
      if (v.size() != matrices.size()) {
        throw std::invalid_argument(
            "tunnelling: ring vector does not match the beads");
      }
      out.resize(v.size());
      for (size_t j = 0; j < v.size(); ++j) {
        if (v[j].size() != matrices[j].cols()) {
          throw std::invalid_argument(
              "tunnelling: ring Hessian direction has the wrong size");
        }
        out[j] = matrices[j] * v[j];
      }
      count->fetch_add(static_cast<long>(v.size()), std::memory_order_relaxed);
    };
  }
  if (hessianProductThreads_ > 1) {
    const long threads =
        std::min(hessianProductThreads_, static_cast<long>(beads.size()));
    std::vector<std::shared_ptr<Potential>> potentials;
    potentials.reserve(static_cast<size_t>(threads));
    for (long thread = 0; thread < threads; ++thread) {
      std::shared_ptr<Potential> clone = potential_->clonePotential();
      if (!clone) {
        throw HessianVectorProductUnavailableError(
            "tunnelling: parallel Hessian products require a clonable "
            "potential");
      }
      if (!clone->supportsMassWeightedHessianVectorProduct()) {
        throw HessianVectorProductUnavailableError(
            "tunnelling: a potential clone has no exact mass-weighted "
            "Hessian-vector product");
      }
      potentials.push_back(std::move(clone));
    }
    std::vector<std::shared_ptr<Matter>> matters;
    matters.reserve(beads.size());
    for (const VectorXd &q : beads) {
      auto matter = std::make_shared<Matter>(reference_);
      coordinates_.place(q, *matter);
      matters.push_back(std::move(matter));
    }
    const VectorXi freeAtoms = coordinates_.freeAtoms();
    const long freedom = coordinates_.dimension();
    const auto count = products_;
    return [potentials = std::move(potentials), matters = std::move(matters),
            freeAtoms, freedom, count, threads](const std::vector<VectorXd> &v,
                                                std::vector<VectorXd> &out) {
      if (v.size() != matters.size()) {
        throw std::invalid_argument(
            "tunnelling: ring vector does not match the beads");
      }
      out.resize(v.size());
      std::vector<std::exception_ptr> failures(static_cast<size_t>(threads));
      std::vector<std::jthread> workers;
      workers.reserve(static_cast<size_t>(threads));
      for (long thread = 0; thread < threads; ++thread) {
        workers.emplace_back([&, thread]() {
          try {
            for (size_t j = static_cast<size_t>(thread); j < v.size();
                 j += static_cast<size_t>(threads)) {
              out[j] = potentials[static_cast<size_t>(thread)]
                           ->massWeightedHessianVectorProduct(*matters[j], v[j],
                                                              freeAtoms);
              count->fetch_add(1, std::memory_order_relaxed);
              if (out[j].size() != freedom || !out[j].allFinite()) {
                throw std::runtime_error(
                    "tunnelling: potential returned an invalid "
                    "Hessian-vector product");
              }
            }
          } catch (...) {
            failures[static_cast<size_t>(thread)] = std::current_exception();
          }
        });
      }
      workers.clear();
      for (const std::exception_ptr &failure : failures) {
        if (failure) {
          std::rethrow_exception(failure);
        }
      }
    };
  }
  std::vector<tunneling::HessianProduct> at;
  at.reserve(beads.size());
  for (const VectorXd &q : beads) {
    at.push_back(hessianProductAt(q));
  }
  return [at = std::move(at)](const std::vector<VectorXd> &v,
                              std::vector<VectorXd> &out) {
    if (v.size() != at.size()) {
      throw std::invalid_argument(
          "tunnelling: ring vector does not match the beads");
    }
    out.resize(v.size());
    for (size_t j = 0; j < v.size(); ++j) {
      out[j] = at[j](v[j]);
    }
  };
}

bool TunnelingJobContext::usesFiniteDifferenceHessianFallback() const noexcept {
  return !potential_->supportsMassWeightedHessianVectorProduct();
}

long TunnelingJobContext::hessianVectorProducts() const {
  return products_->load(std::memory_order_relaxed);
}

TunnelingJobContext::LowestMode
TunnelingJobContext::lowestModeAt(const VectorXd &q, const VectorXd &seed) {
  const long freedom = coordinates_.dimension();
  if (q.size() != freedom || seed.size() != freedom || !q.allFinite() ||
      !seed.allFinite() || !(seed.norm() > 0.0)) {
    throw std::invalid_argument(
        "tunnelling: matrix-free lowest mode needs finite matching vectors");
  }
  return lowestMode(hessianProductAt(q), seed, MatrixXd(freedom, 0));
}

TunnelingJobContext::LowestMode
TunnelingJobContext::lowestMode(const tunneling::HessianProduct &product,
                                const VectorXd &seed,
                                const MatrixXd &excluded) {
  const long freedom = seed.size();
  const long requested = parameters_.lanczos_options().max_iterations;
  if (requested < 1) {
    throw std::invalid_argument(
        "tunnelling: matrix-free lowest mode needs Lanczos iterations");
  }
  auto project = [&](VectorXd &v) {
    for (int pass = 0; pass < 2; ++pass) {
      if (excluded.cols() > 0) {
        v -= excluded * (excluded.transpose() * v);
      }
    }
  };
  VectorXd start = seed;
  project(start);
  if (!(start.norm() > 0.0)) {
    throw std::invalid_argument(
        "tunnelling: Lanczos seed lies in the excluded directions");
  }
  const long columns = std::min(freedom - excluded.cols(), requested);
  std::vector<VectorXd> basis{start.normalized()};
  std::vector<VectorXd> products;
  basis.reserve(static_cast<size_t>(columns));
  products.reserve(static_cast<size_t>(columns));

  LowestMode result;
  for (long i = 0; i < columns; ++i) {
    VectorXd image = product(basis[static_cast<size_t>(i)]);
    project(image);
    products.push_back(std::move(image));
    MatrixXd projected = MatrixXd::Zero(i + 1, i + 1);
    for (long row = 0; row <= i; ++row) {
      for (long column = 0; column <= i; ++column) {
        projected(row, column) = basis[static_cast<size_t>(row)].dot(
            products[static_cast<size_t>(column)]);
      }
    }
    const MatrixXd symmetric = 0.5 * (projected + projected.transpose());
    const Eigen::SelfAdjointEigenSolver<MatrixXd> eigensystem(symmetric);
    if (eigensystem.info() != Eigen::Success) {
      throw std::runtime_error(
          "tunnelling: matrix-free Lanczos eigensystem failed");
    }
    const VectorXd coefficients = eigensystem.eigenvectors().col(0);
    result.curvature = eigensystem.eigenvalues()(0);
    result.vector = VectorXd::Zero(freedom);
    VectorXd ritzImage = VectorXd::Zero(freedom);
    for (long column = 0; column <= i; ++column) {
      result.vector +=
          coefficients(column) * basis[static_cast<size_t>(column)];
      ritzImage += coefficients(column) * products[static_cast<size_t>(column)];
    }
    result.residual = (ritzImage - result.curvature * result.vector).norm();
    result.iterations = i + 1;
    const double scale = std::abs(result.curvature);
    // A Krylov space that fills the complement holds the exact spectrum.
    if (result.residual <= parameters_.lanczos_options().tolerance * scale ||
        i + 1 == freedom - excluded.cols()) {
      return result;
    }
    if (i + 1 == columns) {
      break;
    }
    VectorXd next = products.back();
    for (int pass = 0; pass < 2; ++pass) {
      for (long j = 0; j <= i; ++j) {
        next -= basis[static_cast<size_t>(j)].dot(next) *
                basis[static_cast<size_t>(j)];
      }
      project(next);
    }
    const double norm = next.norm();
    const double precision =
        std::numeric_limits<double>::epsilon() * products.back().norm();
    if (norm <= precision) {
      // An invariant subspace: its Ritz values are exact eigenvalues, but
      // the lowest of them need not be the lowest overall.
      throw std::runtime_error(
          "tunnelling: matrix-free Lanczos basis exhausted before "
          "convergence");
    }
    basis.push_back(next / norm);
  }
  throw std::runtime_error(
      "tunnelling: matrix-free Lanczos lowest mode did not converge");
}

void TunnelingJobContext::findRotationZeroModes(const VectorXd &q) {
  Matter matter(reference_);
  coordinates_.place(q, matter);
  const MatrixXd generators = coordinates_.rigidGenerators(matter);
  const tunneling::HessianProduct product = hessianProductAt(q);
  const Eigen::ColPivHouseholderQR<MatrixXd> qr(generators);
  const long rank = qr.rank();
  const long freedom = coordinates_.dimension();
  rotationZero_ = {false, false, false};
  rotationResidual_ = {0.0, 0.0, 0.0};
  rotationsKnown_ = true;
  if (rank >= freedom) {
    return;
  }
  const MatrixXd span =
      MatrixXd(qr.householderQ()) * MatrixXd::Identity(freedom, rank);
  VectorXd seed(freedom);
  for (long i = 0; i < freedom; ++i) {
    // Any start with a component off the generators reaches the softest
    // vibration; this one is deterministic.
    seed(i) = std::sin(static_cast<double>(i + 1));
  }
  const double softest = lowestMode(product, seed, span).curvature;
  for (int c = 0; c < 3; ++c) {
    const VectorXd r = generators.col(3 + c);
    const double rr = r.squaredNorm();
    if (!(rr > 0.0)) {
      continue;
    }
    const double rq = r.dot(product(r)) / rr;
    if (softest > 0.0) {
      rotationResidual_[static_cast<size_t>(c)] = rq / softest;
      rotationZero_[static_cast<size_t>(c)] =
          std::abs(rq) <= tunneling::kRotationZeroFraction * softest;
    } else {
      rotationResidual_[static_cast<size_t>(c)] =
          std::numeric_limits<double>::infinity();
    }
  }
}

const std::array<bool, 3> &TunnelingJobContext::rotationZero() const {
  return rotationZero_;
}

const std::array<double, 3> &TunnelingJobContext::rotationResidual() const {
  return rotationResidual_;
}

} // namespace eonc

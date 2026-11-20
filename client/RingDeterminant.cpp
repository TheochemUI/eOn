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
#include "eon/RingDeterminant.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>

namespace eonc::tunneling {
namespace {

constexpr double kPrecision = std::numeric_limits<double>::epsilon();

/// Column-major, unlike eonc::MatrixXd: ring vectors map onto it as one
/// column per bead, and the Krylov bases grow by columns.
using Dense = Eigen::MatrixXd;

/// Standard normal deviates by Box-Muller on mt19937_64, whose output the
/// standard fixes, so a seed gives the same probes on every platform.
class Gaussian {
public:
  explicit Gaussian(std::uint64_t seed)
      : engine_(seed) {}

  VectorXd vector(long n) {
    VectorXd v(n);
    for (long i = 0; i < n; ++i) {
      v(i) = next();
    }
    return v;
  }

private:
  double uniform() { return static_cast<double>(engine_() >> 11) * 0x1.0p-53; }

  double next() {
    if (cached_) {
      cached_ = false;
      return spare_;
    }
    const double radius = std::sqrt(-2.0 * std::log(1.0 - uniform()));
    const double angle = 2.0 * std::numbers::pi * uniform();
    spare_ = radius * std::sin(angle);
    cached_ = true;
    return radius * std::cos(angle);
  }

  std::mt19937_64 engine_;
  double spare_{0.0};
  bool cached_{false};
};

/// Removes the span of `basis` from v, twice for orthogonality to working
/// precision.
void project(const Dense &basis, VectorXd &v) {
  if (basis.cols() == 0) {
    return;
  }
  for (int pass = 0; pass < 2; ++pass) {
    v -= basis * (basis.transpose() * v);
  }
}

/// Orthonormalises the columns of `block` against `fixed`, `basis` and each
/// other. A column that vanishes to working precision against what
/// precedes it is dropped: the Krylov space is invariant along it.
Dense orthonormalBlock(const Dense &fixed, const Dense &basis,
                       const Dense &block) {
  Dense out(block.rows(), 0);
  for (long col = 0; col < block.cols(); ++col) {
    VectorXd v = block.col(col);
    const double before = v.norm();
    for (int pass = 0; pass < 2; ++pass) {
      project(fixed, v);
      project(basis, v);
      project(out, v);
    }
    const double norm = v.norm();
    if (!(norm >
          kPrecision * before * std::sqrt(static_cast<double>(block.rows())))) {
      continue;
    }
    out.conservativeResize(Eigen::NoChange, out.cols() + 1);
    out.col(out.cols() - 1) = v / norm;
  }
  return out;
}

/// A Gauss quadrature of z^T g(S) z from the Lanczos tridiagonal matrix of
/// S started at z / |z|.
struct Quadrature {
  double value{0.0};
  double squares{0.0};
  /// sum_i g(theta_i): the trace once the Krylov space is the whole space.
  double trace{0.0};
  double lowest{0.0};
  double nodes{0.0};
};

template <typename Nodes>
Quadrature gaussQuadrature(const std::vector<double> &alpha,
                           const std::vector<double> &beta, double weight,
                           const Nodes &g) {
  const long m = static_cast<long>(alpha.size());
  Dense t = Dense::Zero(m, m);
  for (long i = 0; i < m; ++i) {
    t(i, i) = alpha[static_cast<size_t>(i)];
    if (i + 1 < m) {
      t(i, i + 1) = t(i + 1, i) = beta[static_cast<size_t>(i)];
    }
  }
  const Eigen::SelfAdjointEigenSolver<Dense> es(t);
  Quadrature out;
  out.lowest = es.eigenvalues()(0);
  if (!(out.lowest > 0.0)) {
    return out;
  }
  for (long i = 0; i < m; ++i) {
    const double w = weight * es.eigenvectors()(0, i) * es.eigenvectors()(0, i);
    const double value = g(es.eigenvalues()(i));
    out.value += w * value;
    out.squares += w * value * value;
    out.trace += value;
    out.nodes += es.eigenvalues()(i);
  }
  return out;
}

/// One probe: Lanczos with full reorthogonalisation on the operator, from
/// z, until the quadrature's geometric tail falls to `tail` or the Krylov
/// space is exhausted. `orthogonal` is the space the operator acts on the
/// complement of. A Krylov space of the full `dimension` makes the Ritz
/// values the spectrum, so `trace` and `curvature` are then exact.
struct Probe {
  double value{0.0};
  double squares{0.0};
  double tail{0.0};
  double trace{0.0};
  double curvature{0.0};
  double rayleigh{0.0};
  double lowest{0.0};
  long steps{0};
  bool indefinite{false};
  bool filled{false};
};

template <typename Operator, typename Nodes>
Probe lanczosProbe(const Operator &apply, VectorXd z, const Dense &orthogonal,
                   long dimension, double tail, const Nodes &g) {
  Probe out;
  project(orthogonal, z);
  const double weight = z.squaredNorm();
  Dense v(z.size(), 0);
  std::vector<double> alpha;
  std::vector<double> beta;
  VectorXd current = z / std::sqrt(weight);
  double previous = 0.0;
  double previousDrop = 0.0;
  bool havePrevious = false;
  bool haveDrop = false;
  for (long step = 0; step < dimension; ++step) {
    v.conservativeResize(Eigen::NoChange, v.cols() + 1);
    v.col(v.cols() - 1) = current;
    VectorXd w = apply(current);
    const double scale = w.norm();
    alpha.push_back(current.dot(w));
    if (step == 0) {
      out.rayleigh = alpha.back();
    }
    for (int pass = 0; pass < 2; ++pass) {
      project(v, w);
      project(orthogonal, w);
    }
    const Quadrature q = gaussQuadrature(alpha, beta, weight, g);
    out.steps = step + 1;
    out.lowest = q.lowest;
    if (!(q.lowest > 0.0)) {
      out.indefinite = true;
      return out;
    }
    out.value = q.value;
    out.squares = q.squares;
    const double norm = w.norm();
    if (out.steps == dimension ||
        !(norm >
          kPrecision * scale * std::sqrt(static_cast<double>(z.size())))) {
      out.tail = 0.0;
      out.filled = out.steps == dimension;
      out.trace = q.trace;
      out.curvature = q.nodes / static_cast<double>(dimension);
      return out;
    }
    if (havePrevious) {
      // Once the extreme Ritz values settle, extrapolate the successive
      // Gauss-rule drops as a geometric series. This is a tail estimate,
      // not a certified Gauss--Radau bound.
      const double drop = previous - q.value;
      if (haveDrop && drop > 0.0 && drop < previousDrop) {
        const double ratio = drop / previousDrop;
        out.tail = drop * ratio / (1.0 - ratio);
        if (out.tail <= tail) {
          return out;
        }
      }
      previousDrop = drop;
      haveDrop = true;
    }
    previous = q.value;
    havePrevious = true;
    beta.push_back(norm);
    current = w / norm;
  }
  return out;
}

/// The probe-mean statistics, with the error split of the header.
struct ProbeMean {
  double sum{0.0};
  double squares{0.0};
  double tails{0.0};
  long steps{0};
  long count{0};

  void add(const Probe &p) {
    sum += p.value;
    squares += p.squares;
    tails += p.tail;
    steps += p.steps;
    ++count;
  }
  [[nodiscard]] double mean() const { return sum / static_cast<double>(count); }
  /// For an exact action, Var(z^T A z) = 2 |A|_F^2 for Gaussian z and
  /// z^T A^2 z estimates |A|_F^2 without bias. The Lanczos quadrature
  /// supplies the corresponding finite-Krylov estimate here.
  [[nodiscard]] double variance() const {
    return 2.0 * squares / static_cast<double>(count);
  }
  [[nodiscard]] double standardError() const {
    return std::sqrt(variance() / static_cast<double>(count));
  }
  [[nodiscard]] double quadratureError() const {
    return tails / static_cast<double>(count);
  }
  [[nodiscard]] double probesFor(double error) const {
    return variance() / (error * error);
  }
  [[nodiscard]] double meanSteps() const {
    return static_cast<double>(steps) / static_cast<double>(count);
  }
};

void requireTolerance(double tolerance) {
  if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
    throw std::invalid_argument(
        "log-determinant tolerance must be positive and finite");
  }
}

/// The ring in Fourier-preconditioned coordinates: J~ = M^{-1/2} J M^{-1/2}
/// on the complement of the rigid null space, applied matrix-free.
class PreconditionedRing {
public:
  PreconditionedRing(const RingHessianProducts &products, long beads,
                     long freedom, double c, double shift)
      : products_(products),
        beads_(beads),
        freedom_(freedom),
        c_(c),
        fourier_(Dense::Zero(beads, beads)),
        nu_(beads) {
    const double norm = std::sqrt(1.0 / static_cast<double>(beads));
    const double pair = std::sqrt(2.0 / static_cast<double>(beads));
    for (long j = 0; j < beads; ++j) {
      fourier_(j, 0) = norm;
    }
    nu_(0) = 0.0;
    long col = 1;
    for (long k = 1; 2 * k < beads; ++k) {
      const double s = std::sin(std::numbers::pi * static_cast<double>(k) /
                                static_cast<double>(beads));
      for (long j = 0; j < beads; ++j) {
        const double angle = 2.0 * std::numbers::pi *
                             static_cast<double>(k * j) /
                             static_cast<double>(beads);
        fourier_(j, col) = pair * std::cos(angle);
        fourier_(j, col + 1) = pair * std::sin(angle);
      }
      nu_(col) = nu_(col + 1) = 4.0 * c * s * s;
      col += 2;
    }
    if (beads % 2 == 0) {
      for (long j = 0; j < beads; ++j) {
        fourier_(j, col) = (j % 2 == 0 ? norm : -norm);
      }
      nu_(col) = 4.0 * c;
    }
    const VectorXd m = (nu_.array() + shift).matrix();
    logDetM_ = static_cast<double>(freedom) * m.array().log().sum();
    halfInverse_ = m.array().rsqrt().matrix();
    half_ = m.array().sqrt().matrix();
    inverse_ = m.cwiseInverse();
  }

  [[nodiscard]] long dimension() const { return beads_ * freedom_; }
  [[nodiscard]] double logDetM() const { return logDetM_; }
  [[nodiscard]] long products() const { return applications_; }

  [[nodiscard]] VectorXd halfInverse(const VectorXd &x) const {
    return fourierScale(x, halfInverse_);
  }
  [[nodiscard]] VectorXd half(const VectorXd &x) const {
    return fourierScale(x, half_);
  }
  [[nodiscard]] VectorXd inverse(const VectorXd &x) const {
    return fourierScale(x, inverse_);
  }

  /// J x in Euclidean coordinates.
  VectorXd ring(const VectorXd &x) {
    std::vector<VectorXd> in(static_cast<size_t>(beads_));
    for (long j = 0; j < beads_; ++j) {
      in[static_cast<size_t>(j)] = x.segment(j * freedom_, freedom_);
    }
    std::vector<VectorXd> out(static_cast<size_t>(beads_));
    products_(in, out);
    ++applications_;
    VectorXd y(x.size());
    for (long j = 0; j < beads_; ++j) {
      const auto &h = out[static_cast<size_t>(j)];
      if (h.size() != freedom_ || !h.allFinite()) {
        throw std::runtime_error(
            "ring Hessian product is missing or not finite at bead " +
            std::to_string(j));
      }
      const long prev = (j + beads_ - 1) % beads_;
      const long next = (j + 1) % beads_;
      y.segment(j * freedom_, freedom_) =
          h +
          c_ * (2.0 * in[static_cast<size_t>(j)] -
                in[static_cast<size_t>(prev)] - in[static_cast<size_t>(next)]);
    }
    return y;
  }

  VectorXd preconditioned(const VectorXd &x) {
    return halfInverse(ring(halfInverse(x)));
  }

private:
  [[nodiscard]] VectorXd fourierScale(const VectorXd &x,
                                      const VectorXd &scale) const {
    const Eigen::Map<const Dense> beadsByColumn(x.data(), freedom_, beads_);
    Dense modes = beadsByColumn * fourier_;
    modes *= scale.asDiagonal();
    VectorXd y(x.size());
    Eigen::Map<Dense>(y.data(), freedom_, beads_) =
        modes * fourier_.transpose();
    return y;
  }

  const RingHessianProducts &products_;
  long beads_;
  long freedom_;
  double c_;
  Dense fourier_;
  VectorXd nu_;
  VectorXd halfInverse_;
  VectorXd half_;
  VectorXd inverse_;
  double logDetM_{0.0};
  long applications_{0};
};

VectorXd flatten(const std::vector<VectorXd> &beads, long freedom) {
  VectorXd out(static_cast<long>(beads.size()) * freedom);
  for (size_t j = 0; j < beads.size(); ++j) {
    if (beads[j].size() != freedom) {
      throw std::invalid_argument("ring vector bead has the wrong size");
    }
    out.segment(static_cast<long>(j) * freedom, freedom) = beads[j];
  }
  return out;
}

} // namespace

StochasticLogDet
ringLogDetPrime(const RingHessianProducts &products, double c, double shift,
                const std::vector<VectorXd> &cycle,
                const std::vector<std::vector<VectorXd>> &rigid,
                double tolerance, std::uint64_t seed) {
  requireTolerance(tolerance);
  const long beads = static_cast<long>(cycle.size());
  if (beads < 3 || cycle.front().size() == 0) {
    throw std::invalid_argument("ring determinant needs at least three beads");
  }
  if (!(c > 0.0) || !(shift > 0.0)) {
    throw std::invalid_argument(
        "ring determinant needs positive spring and shift curvatures");
  }
  const long freedom = cycle.front().size();
  PreconditionedRing ring(products, beads, freedom, c, shift);
  const long n = ring.dimension();

  Dense euclideanRigid(n, static_cast<long>(rigid.size()));
  for (size_t r = 0; r < rigid.size(); ++r) {
    euclideanRigid.col(static_cast<long>(r)) = flatten(rigid[r], freedom);
  }
  // The rigid null space of J maps to M^{1/2} R in J~; its Gram matrix
  // R^T M R carries the Jacobian of that change of basis.
  Dense lifted(n, euclideanRigid.cols());
  for (long r = 0; r < euclideanRigid.cols(); ++r) {
    lifted.col(r) = ring.half(euclideanRigid.col(r));
  }
  double logDetRigid = 0.0;
  if (lifted.cols() > 0) {
    const Dense gram = lifted.transpose() * lifted;
    const Eigen::LLT<Dense> llt(gram);
    if (llt.info() != Eigen::Success) {
      throw std::invalid_argument("ring rigid vectors are linearly dependent");
    }
    logDetRigid = 2.0 * llt.matrixLLT().diagonal().array().log().sum();
  }
  const Dense rigidPre = orthonormalBlock(Dense(n, 0), Dense(n, 0), lifted);
  if (rigidPre.cols() != lifted.cols()) {
    throw std::invalid_argument("ring rigid vectors are linearly dependent");
  }
  const long space = n - rigidPre.cols();
  if (space < 2) {
    throw std::invalid_argument("ring has no directions beyond its rigid ones");
  }

  VectorXd tau = flatten(cycle, freedom);
  project(euclideanRigid, tau);
  const double tauNorm = tau.norm();
  if (!(tauNorm > 0.0)) {
    throw std::runtime_error("the ring cycle lies in its rigid null space");
  }
  tau /= tauNorm;

  Gaussian gaussian(seed);
  auto apply = [&](const VectorXd &x) {
    VectorXd y = x;
    project(rigidPre, y);
    y = ring.preconditioned(y);
    project(rigidPre, y);
    return y;
  };

  // The cycle and one Gaussian column enter the determinant exactly.
  constexpr long kBlock = 2;
  Dense start(n, kBlock);
  start.col(0) = ring.half(tau);
  start.col(1) = gaussian.vector(n);
  Dense basis = orthonormalBlock(rigidPre, Dense(n, 0), start);
  Dense images(n, basis.cols());
  for (long col = 0; col < basis.cols(); ++col) {
    images.col(col) = apply(basis.col(col));
  }
  long lastBlock = basis.cols();
  bool invariant = false;
  long targetBlocks = 1;
  long blocks = 1;
  StochasticLogDet result;

  for (;;) {
    while (blocks < targetBlocks && !invariant && basis.cols() < space) {
      const Dense fresh =
          orthonormalBlock(rigidPre, basis, images.rightCols(lastBlock));
      if (fresh.cols() == 0) {
        invariant = true;
        break;
      }
      const long old = basis.cols();
      const long take = std::min(fresh.cols(), space - old);
      basis.conservativeResize(Eigen::NoChange, old + take);
      images.conservativeResize(Eigen::NoChange, old + take);
      for (long col = 0; col < take; ++col) {
        basis.col(old + col) = fresh.col(col);
        images.col(old + col) = apply(fresh.col(col));
      }
      lastBlock = take;
      ++blocks;
    }
    const long d = basis.cols();
    const bool canGrow = !invariant && d < space;

    // The projected block T = U^T J~ U and its Euclidean Rayleigh-Ritz
    // pencil (T, X^T X) with X = M^{-1/2} U off the rigid space.
    Dense t = basis.transpose() * images;
    t = 0.5 * (t + t.transpose()).eval();
    const Eigen::SelfAdjointEigenSolver<Dense> blockEigen(t);
    const VectorXd theta = blockEigen.eigenvalues();
    Dense x(n, d);
    for (long col = 0; col < d; ++col) {
      VectorXd column = ring.halfInverse(basis.col(col));
      project(euclideanRigid, column);
      x.col(col) = column;
    }
    Dense gram = x.transpose() * x;
    gram = 0.5 * (gram + gram.transpose()).eval();
    const Eigen::GeneralizedSelfAdjointEigenSolver<Dense> pencil(
        t, gram, Eigen::ComputeEigenvectors | Eigen::Ax_lBx);
    if (pencil.info() != Eigen::Success) {
      throw std::runtime_error("ring Rayleigh-Ritz pencil failed");
    }
    const VectorXd lambda = pencil.eigenvalues();
    const Dense ritz = x * pencil.eigenvectors();
    long zeroIndex = 0;
    for (long i = 1; i < d; ++i) {
      if (std::abs(tau.dot(ritz.col(i))) >
          std::abs(tau.dot(ritz.col(zeroIndex)))) {
        zeroIndex = i;
      }
    }
    const double zero = lambda(zeroIndex);
    if (!(std::abs(zero) > 0.0)) {
      throw std::runtime_error("the ring's cycle eigenvalue is exactly zero");
    }
    VectorXd residual =
        ring.half(images * pencil.eigenvectors().col(zeroIndex));
    project(euclideanRigid, residual);
    residual -= zero * ritz.col(zeroIndex);
    double lowest = std::numeric_limits<double>::infinity();
    for (long i = 0; i < d; ++i) {
      if (i != zeroIndex) {
        lowest = std::min(lowest, lambda(i));
      }
    }
    // theta - lambda_0 = sum_i r_i^2 / (lambda_i - theta) to second order
    // in the residual r, over the other eigenpairs. The kinetic part
    // dominates the high-frequency lambda_i, so |J - theta| >= gamma M off
    // the cycle with gamma the gap of the pencil (J - theta, M), whose Ritz
    // values T - theta X^T X the block resolves; the estimate is then
    // r^T M^{-1} r / gamma. The Euclidean |r| alone weighs the
    // high-frequency residual by nu_max, which hardly moves theta.
    Dense shiftedPencil = t - zero * gram;
    shiftedPencil = 0.5 * (shiftedPencil + shiftedPencil.transpose()).eval();
    const VectorXd shiftedValues = Eigen::SelfAdjointEigenSolver<Dense>(
                                       shiftedPencil, Eigen::EigenvaluesOnly)
                                       .eigenvalues();
    long nearest = 0;
    for (long i = 1; i < d; ++i) {
      if (std::abs(shiftedValues(i)) < std::abs(shiftedValues(nearest))) {
        nearest = i;
      }
    }
    double gamma = std::numeric_limits<double>::infinity();
    for (long i = 0; i < d; ++i) {
      if (i != nearest) {
        gamma = std::min(gamma, std::abs(shiftedValues(i)));
      }
    }
    const double reach =
        residual.dot(ring.inverse(residual)) / (gamma * std::abs(zero));
    const double modeError = (d == space) ? 0.0
                             : reach < 1.0
                                 ? -std::log1p(-reach)
                                 : std::numeric_limits<double>::infinity();
    long negative = 0;
    double logDetT = 0.0;
    for (long i = 0; i < d; ++i) {
      negative += theta(i) < 0.0 ? 1 : 0;
      logDetT += std::log(std::abs(theta(i)));
    }

    result.deflated = d;
    result.zeroEigenvalue = zero;
    result.lowestEigenvalue = lowest;
    result.negativeModes = negative - (zero < 0.0 ? 1 : 0);
    result.modeError = modeError;
    result.exhausted = d == space;
    const double fixed =
        ring.logDetM() - logDetRigid + logDetT - std::log(std::abs(zero));
    if (result.exhausted) {
      result.value = fixed;
      result.standardError = 0.0;
      result.quadratureError = 0.0;
      result.modeError = 0.0;
      result.products = ring.products();
      return result;
    }
    if (!(modeError <= 0.25 * tolerance) && canGrow) {
      targetBlocks = 2 * blocks;
      continue;
    }

    // The Schur complement of the block on its orthogonal complement,
    // S = P J~ P - E T^{-1} E^T with E = (I - U U^T) J~ U. Its inertia is
    // J~'s less T's, so S is positive definite exactly when the block
    // holds every negative direction.
    const Dense coupling = (images - basis * t) * blockEigen.eigenvectors();
    const VectorXd inverseTheta = theta.cwiseInverse();
    Dense orthogonal(n, rigidPre.cols() + d);
    orthogonal << rigidPre, basis;
    auto schur = [&](const VectorXd &w) {
      VectorXd y = apply(w);
      project(basis, y);
      y -= coupling * (inverseTheta.asDiagonal() * (coupling.transpose() * w));
      return y;
    };
    auto logarithm = [](double value) { return std::log(value); };
    ProbeMean mean;
    bool indefinite = false;
    bool grow = false;
    for (;;) {
      const Probe probe = lanczosProbe(schur, gaussian.vector(n), orthogonal,
                                       space - d, 0.25 * tolerance, logarithm);
      if (probe.indefinite) {
        indefinite = true;
        break;
      }
      if (probe.filled) {
        result.value = fixed + probe.trace;
        result.standardError = 0.0;
        result.quadratureError = 0.0;
        result.modeError = 0.0;
        result.probes = mean.count + 1;
        result.products = ring.products();
        result.exhausted = true;
        return result;
      }
      mean.add(probe);
      const double needed = mean.probesFor(0.5 * tolerance);
      if (canGrow && needed * mean.meanSteps() > static_cast<double>(d)) {
        grow = true;
        break;
      }
      if (mean.standardError() <= 0.5 * tolerance) {
        break;
      }
    }
    if (indefinite && !canGrow) {
      throw std::runtime_error(
          "ring Schur complement is indefinite and the Krylov block cannot "
          "grow");
    }
    if (indefinite || grow) {
      targetBlocks = 2 * blocks;
      continue;
    }
    result.value = fixed + mean.mean();
    result.standardError = mean.standardError();
    result.quadratureError = mean.quadratureError();
    result.probes = mean.count;
    result.products = ring.products();
    return result;
  }
}

StochasticLogDet wellLogDet(const HessianProduct &product, long freedom,
                            const MatrixXd &rigid,
                            const std::vector<double> &shifts, double tolerance,
                            std::uint64_t seed) {
  requireTolerance(tolerance);
  if (freedom < 1 || rigid.rows() != freedom || shifts.empty()) {
    throw std::invalid_argument("well determinant inputs are inconsistent");
  }
  const Dense basis =
      orthonormalBlock(Dense(freedom, 0), Dense(freedom, 0), Dense(rigid));
  if (basis.cols() != rigid.cols()) {
    throw std::invalid_argument("well rigid vectors are linearly dependent");
  }
  const long space = freedom - basis.cols();
  if (space < 1) {
    throw std::invalid_argument("well has no directions beyond its rigid ones");
  }
  StochasticLogDet result;
  auto apply = [&](const VectorXd &v) {
    VectorXd y = product(v);
    if (y.size() != freedom || !y.allFinite()) {
      throw std::runtime_error(
          "Hessian product is missing or not finite at the well");
    }
    ++result.products;
    project(basis, y);
    return y;
  };
  auto shifted = [&](double value) {
    double sum = 0.0;
    for (const double s : shifts) {
      sum += std::log(value + s);
    }
    return sum;
  };
  Gaussian gaussian(seed);
  ProbeMean mean;
  double rayleigh = 0.0;
  bool fill = false;
  result.lowestEigenvalue = std::numeric_limits<double>::infinity();
  for (;;) {
    // A zero tail target runs the probe until its Krylov space is
    // invariant, which the whole complement is.
    const Probe probe =
        lanczosProbe(apply, gaussian.vector(freedom), basis, space,
                     fill ? 0.0 : 0.25 * tolerance, shifted);
    if (probe.indefinite) {
      throw std::runtime_error(
          "the reactant Hessian is not positive definite off its rigid modes");
    }
    result.lowestEigenvalue = std::min(result.lowestEigenvalue, probe.lowest);
    if (probe.filled) {
      result.value = probe.trace;
      result.meanCurvature = probe.curvature;
      result.probes = mean.count + 1;
      result.exhausted = true;
      return result;
    }
    mean.add(probe);
    // z^T H z / z^T z has mean tr(H) / dimension for an isotropic z.
    rayleigh += probe.rayleigh;
    if (mean.standardError() <= 0.5 * tolerance) {
      break;
    }
    // The probes still needed cost more than one that fills the space.
    fill = (mean.probesFor(0.5 * tolerance) - static_cast<double>(mean.count)) *
               mean.meanSteps() >
           static_cast<double>(space);
  }
  result.value = mean.mean();
  result.standardError = mean.standardError();
  result.quadratureError = mean.quadratureError();
  result.probes = mean.count;
  result.meanCurvature = rayleigh / static_cast<double>(mean.count);
  return result;
}

} // namespace eonc::tunneling

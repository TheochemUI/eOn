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
**
** The half-ring springs and the banded ring Hessian are adapted from i-PI
** under the MIT licence.
** i-PI Copyright (C) 2014-2015 i-PI developers
** Algorithms implemented by Yair Litman and Mariana Rossi, 2017.
*/
#include "eon/Tunneling.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace eonc::tunneling {

double massWeightedDistance(const Matter &a, const Matter &b) {
  if (a.numberOfAtoms() != b.numberOfAtoms()) {
    throw std::invalid_argument("the structures hold different atom counts");
  }
  const AtomMatrix dr = a.pbc(b.getPositions() - a.getPositions());
  const auto mass = a.getMasses();
  double sum = 0.0;
  for (long i = 0; i < a.numberOfAtoms(); ++i) {
    if (mass(i) <= 0.0) {
      throw std::invalid_argument(
          "every atom needs a positive mass for a mass-weighted path");
    }
    sum += mass(i) * dr.row(i).squaredNorm();
  }
  return std::sqrt(sum);
}

std::vector<double>
massWeightedPath(const std::vector<std::shared_ptr<Matter>> &band) {
  std::vector<double> s{0.0};
  s.reserve(band.size());
  for (size_t i = 1; i < band.size(); ++i) {
    s.push_back(s.back() + massWeightedDistance(*band[i - 1], *band[i]));
  }
  return s;
}

Profile::Profile(std::vector<double> s, std::vector<double> v)
    : s_(std::move(s)),
      v_(std::move(v)),
      m_(s_.size(), 0.0) {
  const size_t n = s_.size();
  if (n < 2 || v_.size() != n) {
    throw std::invalid_argument(
        "a profile needs matching s and V with two points or more");
  }
  for (size_t k = 1; k < n; ++k) {
    if (!(s_[k] > s_[k - 1])) {
      throw std::invalid_argument(
          "the path coordinate must increase along the band");
    }
  }
  for (size_t k = 1; k + 1 < n; ++k) {
    const double h0 = s_[k] - s_[k - 1];
    const double h1 = s_[k + 1] - s_[k];
    const double d0 = (v_[k] - v_[k - 1]) / h0;
    const double d1 = (v_[k + 1] - v_[k]) / h1;
    if (d0 * d1 <= 0.0) {
      m_[k] = 0.0;
    } else {
      const double w1 = 2.0 * h1 + h0;
      const double w2 = h1 + 2.0 * h0;
      m_[k] = (w1 + w2) / (w1 / d0 + w2 / d1);
    }
  }
}

double Profile::operator()(double x) const {
  x = std::clamp(x, s_.front(), s_.back());
  auto it = std::upper_bound(s_.begin(), s_.end(), x);
  size_t k = static_cast<size_t>(std::distance(s_.begin(), it));
  k = std::clamp<size_t>(k == 0 ? 0 : k - 1, 0, s_.size() - 2);
  const double h = s_[k + 1] - s_[k];
  const double t = (x - s_[k]) / h;
  const double t2 = t * t;
  const double t3 = t2 * t;
  return (2 * t3 - 3 * t2 + 1) * v_[k] + (t3 - 2 * t2 + t) * h * m_[k] +
         (-2 * t3 + 3 * t2) * v_[k + 1] + (t3 - t2) * h * m_[k + 1];
}

double wellCurvature(const Profile &p, bool leftEnd) {
  const auto &s = p.s();
  const auto &v = p.v();
  const size_t n = s.size();
  const double top = *std::max_element(v.begin(), v.end());
  const double floor = leftEnd ? v.front() : v.back();
  const double half = 0.5 * (top - floor);
  // Sums for the normal equations of y = a x^2 + b x^3.
  double s44 = 0, s45 = 0, s55 = 0, sy2 = 0, sy3 = 0;
  size_t used = 0;
  for (size_t j = 1; j < n; ++j) {
    const size_t i = leftEnd ? j : n - 1 - j;
    const double x = leftEnd ? s[i] - s.front() : s.back() - s[i];
    const double y = v[i] - floor;
    if (y > half) {
      break;
    }
    const double x2 = x * x;
    s44 += x2 * x2;
    s45 += x2 * x2 * x;
    s55 += x2 * x2 * x2;
    sy2 += y * x2;
    sy3 += y * x2 * x;
    ++used;
  }
  if (used == 0) {
    // The next image already stands above half the barrier: the parabola
    // through it is all the band says about this well.
    const size_t i = leftEnd ? 1 : n - 2;
    const double x = leftEnd ? s[i] - s.front() : s.back() - s[i];
    return 2.0 * (v[i] - floor) / (x * x);
  }
  if (used == 1) {
    return 2.0 * sy2 / s44;
  }
  const double det = s44 * s55 - s45 * s45;
  const double a = (sy2 * s55 - sy3 * s45) / det;
  return 2.0 * a;
}

double hbarOmega(double curvature) {
  if (!(curvature > 0.0)) {
    throw std::invalid_argument("a well needs a positive curvature");
  }
  return kHbar * std::sqrt(curvature);
}

double wkbAction(const Profile &p, double energy, int points) {
  const double a = p.s().front();
  const double b = p.s().back();
  const double h = (b - a) / (points - 1);
  double sum = 0.0;
  for (int i = 0; i < points; ++i) {
    const double gap = p(a + i * h) - energy;
    const double f = gap > 0.0 ? std::sqrt(2.0 * gap) : 0.0;
    sum += (i == 0 || i == points - 1) ? 0.5 * f : f;
  }
  return sum * h / kHbar;
}

double Splitting::tlsEnergy() const { return std::hypot(delta, delta0); }

Splitting wkbSplitting(const Profile &p, double hwReactant, double hwProduct) {
  const auto &v = p.v();
  Splitting out;
  const double top = *std::max_element(v.begin(), v.end());
  out.delta = v.back() - v.front();
  out.barrier = top - v.front();
  out.hwReactant = hwReactant;
  out.hwProduct = hwProduct;
  out.referenceEnergy =
      std::max(v.front() + 0.5 * hwReactant, v.back() + 0.5 * hwProduct);
  out.action = wkbAction(p, out.referenceEnergy);
  const double hw = std::sqrt(hwReactant * hwProduct);
  out.delta0 = hw / std::numbers::pi * std::exp(-out.action);
  out.deepWells =
      (top - v.front()) > hwReactant && (top - v.back()) > hwProduct;
  return out;
}

namespace {

std::vector<double> arcLengths(const std::vector<VectorXd> &path) {
  std::vector<double> s(path.size(), 0.0);
  for (size_t k = 1; k < path.size(); ++k) {
    s[k] = s[k - 1] + (path[k] - path[k - 1]).norm();
  }
  return s;
}

VectorXd atArcLength(const std::vector<VectorXd> &path,
                     const std::vector<double> &s, double target) {
  const auto it = std::upper_bound(s.begin(), s.end(), target);
  const size_t k = std::clamp<size_t>(
      static_cast<size_t>(std::distance(s.begin(), it)), 1, path.size() - 1);
  const double seg = s[k] - s[k - 1];
  const double t = seg > 0.0 ? (target - s[k - 1]) / seg : 0.0;
  return path[k - 1] + std::clamp(t, 0.0, 1.0) * (path[k] - path[k - 1]);
}

// Crossings V(s) = e nearest the barrier top, one on each side.
std::pair<double, double> turningPoints(const Profile &p, double sTop,
                                        double energy) {
  const double s0 = p.s().front();
  const double s1 = p.s().back();
  auto cross = [&](double from, double to) {
    const int steps = 2000;
    double a = from;
    double b = to;
    for (int k = 1; k <= steps; ++k) {
      const double sk = from + (to - from) * static_cast<double>(k) / steps;
      if (p(sk) <= energy) {
        a = from + (to - from) * static_cast<double>(k - 1) / steps;
        b = sk;
        break;
      }
      if (k == steps) {
        return to;
      }
    }
    for (int k = 0; k < 60; ++k) {
      const double m = 0.5 * (a + b);
      if (p(m) > energy) {
        a = m;
      } else {
        b = m;
      }
    }
    return 0.5 * (a + b);
  };
  return {cross(sTop, s0), cross(sTop, s1)};
}

// int_{s-}^{s+} ds / sqrt(2 (V - E)). The cosine substitution keeps the
// integrand finite at the turning points. Optional tables are the running
// integral and the arc length at the end of each panel.
double halfPeriod(const Profile &p, double sMinus, double sPlus, double energy,
                  std::vector<double> *cumulative = nullptr,
                  std::vector<double> *positions = nullptr) {
  const int panels = 4000;
  const double mid = 0.5 * (sMinus + sPlus);
  const double half = 0.5 * (sPlus - sMinus);
  double total = 0.0;
  if (cumulative != nullptr) {
    cumulative->assign(1, 0.0);
    positions->assign(1, sMinus);
  }
  for (int k = 0; k < panels; ++k) {
    const double phi =
        std::numbers::pi * (static_cast<double>(k) + 0.5) / panels;
    const double s = mid - half * std::cos(phi);
    const double under = 2.0 * (p(s) - energy);
    const double integrand =
        half * std::sin(phi) / std::sqrt(std::max(under, 1e-300));
    total += integrand * std::numbers::pi / panels;
    if (cumulative != nullptr) {
      cumulative->push_back(total);
      positions->push_back(
          mid - half * std::cos(std::numbers::pi * (k + 1.0) / panels));
    }
  }
  return total;
}

} // namespace

std::vector<VectorXd> ringFromPath(const std::vector<VectorXd> &path,
                                   const std::vector<double> &energies,
                                   double betaHbar, long beads) {
  if (path.size() < 3 || path.size() != energies.size() || beads < 4 ||
      !(betaHbar > 0.0)) {
    throw std::invalid_argument(
        "ringFromPath: a path of at least three points with energies, "
        "N >= 4 and beta hbar > 0");
  }
  const long width = path.front().size();
  for (const auto &q : path) {
    if (q.size() != width) {
      throw std::invalid_argument("ringFromPath: the path changes dimension");
    }
  }
  const std::vector<double> s = arcLengths(path);
  const Profile profile(s, energies);
  double sTop = s.front();
  double vTop = -std::numeric_limits<double>::infinity();
  const int grid = 4000;
  for (int k = 0; k <= grid; ++k) {
    const double sk =
        s.front() + (s.back() - s.front()) * static_cast<double>(k) / grid;
    if (profile(sk) > vTop) {
      vTop = profile(sk);
      sTop = sk;
    }
  }
  const double vLow = std::max(energies.front(), energies.back());
  if (!(vTop > vLow)) {
    throw std::invalid_argument("ringFromPath: the path has no barrier");
  }
  auto period = [&](double energy) {
    const auto [sMinus, sPlus] = turningPoints(profile, sTop, energy);
    return 2.0 * halfPeriod(profile, sMinus, sPlus, energy);
  };
  // The crossover along the path from a parabola through the three input
  // points around the barrier top, not from the interpolant, whose slope is
  // clamped to zero at the top node.
  {
    size_t top = 0;
    for (size_t k = 1; k < energies.size(); ++k) {
      if (energies[k] > energies[top]) {
        top = k;
      }
    }
    if (top == 0 || top + 1 == energies.size()) {
      throw std::invalid_argument(
          "ringFromPath: the barrier top is an end of the path");
    }
    const double h1 = s[top] - s[top - 1], h2 = s[top + 1] - s[top];
    const double curvature =
        2.0 *
        (h1 * energies[top + 1] - (h1 + h2) * energies[top] +
         h2 * energies[top - 1]) /
        (h1 * h2 * (h1 + h2));
    if (!(curvature < 0.0)) {
      throw std::invalid_argument("ringFromPath: no curvature at the top");
    }
    const double tc = kHbar * std::sqrt(-curvature) / (2.0 * std::numbers::pi);
    if (!(kHbar / betaHbar < tc)) {
      throw std::invalid_argument(
          "ringFromPath: the temperature is at or above the crossover along "
          "this path");
    }
  }
  // Bracket the orbit energy geometrically above the lower end: the period
  // grows only logarithmically as E approaches a well bottom.
  const double span = vTop - vLow;
  double eHi = vTop - 1e-9 * span;
  double eLo = vLow + 1e-14 * span;
  if (period(eLo) < betaHbar) {
    // The path does not reach a long enough orbit. The lowest one it
    // holds is the start.
    eHi = eLo;
  }
  for (int k = 0; k < 200 && eHi > eLo; ++k) {
    const double e = vLow + std::sqrt((eLo - vLow) * (eHi - vLow));
    if (period(e) > betaHbar) {
      eLo = e;
    } else {
      eHi = e;
    }
    if (eHi - eLo < 1e-15 * span) {
      break;
    }
  }
  const double energy = 0.5 * (eLo + eHi);
  const auto [sMinus, sPlus] = turningPoints(profile, sTop, energy);
  std::vector<double> tau;
  std::vector<double> pos;
  const double half = halfPeriod(profile, sMinus, sPlus, energy, &tau, &pos);
  if (tau.size() < 2 || tau.size() != pos.size()) {
    throw std::invalid_argument("ringFromPath: the orbit has no length");
  }
  // Bead j sits at imaginary time j * beta hbar / N on the way from the
  // reactant-side turning point to the other side. The return repeats it.
  std::vector<VectorXd> ring(static_cast<size_t>(beads), VectorXd::Zero(width));
  for (long j = 0; j <= beads / 2; ++j) {
    const double t =
        std::min(half, half * 2.0 * static_cast<double>(j) / beads);
    const auto it = std::upper_bound(tau.begin(), tau.end(), t);
    const size_t k = std::clamp<size_t>(
        static_cast<size_t>(std::distance(tau.begin(), it)), 1, tau.size() - 1);
    const double seg = tau[k] - tau[k - 1];
    const double w = seg > 0.0 ? (t - tau[k - 1]) / seg : 0.0;
    const double sj =
        pos[k - 1] + std::clamp(w, 0.0, 1.0) * (pos[k] - pos[k - 1]);
    ring[static_cast<size_t>(j)] = atArcLength(path, s, sj);
    if (j > 0 && j < beads - j) {
      ring[static_cast<size_t>(beads - j)] = ring[static_cast<size_t>(j)];
    }
  }
  return ring;
}

double wkbLogRateAlongPath(const Profile &profile, double beta,
                           double hwReactant) {
  if (!(beta > 0.0) || !(hwReactant > 0.0)) {
    throw std::invalid_argument(
        "wkbLogRateAlongPath: beta and hbar omega must be positive");
  }
  const double vReactant = profile.v().front();
  const double s0 = profile.s().front();
  const double s1 = profile.s().back();
  double vTop = vReactant;
  double sTop = s0;
  const int grid = 2000;
  for (int k = 0; k <= grid; ++k) {
    const double sk = s0 + (s1 - s0) * static_cast<double>(k) / grid;
    const double vk = profile(sk);
    if (vk >= vTop) {
      vTop = vk;
      sTop = sk;
    }
  }
  const double barrier = vTop - vReactant;
  if (!(barrier > 0.0)) {
    throw std::invalid_argument(
        "wkbLogRateAlongPath: no barrier above the reactant");
  }
  double ds = 1e-3 * (s1 - s0);
  ds = std::min(ds, std::min(sTop - s0, s1 - sTop));
  if (!(ds > 0.0)) {
    throw std::invalid_argument(
        "wkbLogRateAlongPath: the barrier top is at an end of the path");
  }
  const double curvature =
      std::max(1e-12, -(profile(sTop + ds) - 2.0 * vTop + profile(sTop - ds)) /
                          (ds * ds));
  const double hwBarrier = kHbar * std::sqrt(curvature);
  // int P(E) exp(-beta E) dE from the reactant up to where the Boltzmann
  // factor has died. E is measured from the reactant.
  const double eMax = barrier + 40.0 / beta;
  const int points = 600;
  auto logAdd = [](double a, double b) {
    if (a == -std::numeric_limits<double>::infinity()) {
      return b;
    }
    const double m = std::max(a, b);
    return m + std::log(std::exp(a - m) + std::exp(b - m));
  };
  double logTerms = -std::numeric_limits<double>::infinity();
  double prevLog = -std::numeric_limits<double>::infinity();
  double prevE = 0.0;
  for (int k = 0; k <= points; ++k) {
    const double energy = eMax * static_cast<double>(k) / points;
    double theta = 0.0;
    if (energy < barrier) {
      theta = wkbAction(profile, vReactant + energy);
    } else {
      theta = -std::numbers::pi * (energy - barrier) / hwBarrier;
    }
    const double logP =
        theta > 20.0 ? -2.0 * theta : -std::log1p(std::exp(2.0 * theta));
    const double logF = logP - beta * energy;
    if (k > 0) {
      const double segment =
          std::log(0.5 * (energy - prevE)) + logAdd(prevLog, logF);
      logTerms = logAdd(logTerms, segment);
    }
    prevLog = logF;
    prevE = energy;
  }
  const double logFlux = logTerms - std::log(2.0 * std::numbers::pi * kHbar);
  return logFlux + std::log(2.0 * std::sinh(0.5 * beta * hwReactant));
}

Splitting bandSplitting(const std::vector<std::shared_ptr<Matter>> &band,
                        double referenceEnergy) {
  std::vector<double> v;
  v.reserve(band.size());
  for (const auto &image : band) {
    v.push_back(image->getPotentialEnergy() - referenceEnergy);
  }
  const Profile p(massWeightedPath(band), std::move(v));
  return wkbSplitting(p, hbarOmega(wellCurvature(p, true)),
                      hbarOmega(wellCurvature(p, false)));
}

namespace {

/// Ring coordinates up to which the Newton step takes the dense spectrum
/// and a dense solve; beyond, the block chain and Lanczos.
constexpr long kDenseRing = 4096;

/// Lanczos steps a large ring's Newton view grows to at most; the basis
/// holds that many ring vectors.
// Overlap with the last climb below which the lowest mode takes over.
constexpr double kTrackOverlap = 0.3;
constexpr long kRitzCap = 400;

using ColMajorXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

// J = c T (x) I + blockdiag(A_k), T the tridiagonal (2, -1) spring matrix;
// the diagonal blocks hold 2c already, the off-diagonal blocks are -c I.
// Block LU: D_1 = A_1, D_k = A_k - c^2 D_{k-1}^-1.
class BlockChain {
public:
  BlockChain(double c, const std::vector<MatrixXd> &diag)
      : c_(c) {
    lu_.reserve(diag.size());
    for (size_t k = 0; k < diag.size(); ++k) {
      ColMajorXd d = diag[k];
      if (k > 0) {
        d -= c_ * c_ * lu_.back().inverse();
      }
      lu_.emplace_back(std::move(d));
      const auto &f = lu_.back();
      sign_ *= static_cast<int>(std::lround(f.permutationP().determinant()));
      for (long i = 0; i < d.rows(); ++i) {
        const double u = f.matrixLU()(i, i);
        if (u == 0.0) {
          throw std::runtime_error("instanton: singular chain Hessian block");
        }
        sign_ *= u < 0.0 ? -1 : 1;
        logAbsDet_ += std::log(std::abs(u));
      }
    }
  }
  double logAbsDet() const { return logAbsDet_; }
  int sign() const { return sign_; }
  // x = J^-1 b, b and x stacked by bead.
  std::vector<VectorXd> solve(const std::vector<VectorXd> &b) const {
    const size_t m = b.size();
    std::vector<VectorXd> y(m), x(m);
    y[0] = b[0];
    for (size_t k = 1; k < m; ++k) {
      y[k] = b[k] + c_ * lu_[k - 1].solve(y[k - 1]);
    }
    x[m - 1] = lu_[m - 1].solve(y[m - 1]);
    for (size_t k = m - 1; k-- > 0;) {
      x[k] = lu_[k].solve(y[k] + c_ * x[k + 1]);
    }
    return x;
  }
  // Same recurrence, several right-hand sides per bead.
  std::vector<MatrixXd> solve(const std::vector<MatrixXd> &b) const {
    const size_t m = b.size();
    std::vector<MatrixXd> y(m), x(m);
    y[0] = b[0];
    for (size_t k = 1; k < m; ++k) {
      y[k] = b[k] + c_ * lu_[k - 1].solve(y[k - 1]);
    }
    x[m - 1] = lu_[m - 1].solve(y[m - 1]);
    for (size_t k = m - 1; k-- > 0;) {
      x[k] = lu_[k].solve(y[k] + c_ * x[k + 1]);
    }
    return x;
  }

private:
  double c_;
  std::vector<Eigen::PartialPivLU<ColMajorXd>> lu_;
  double logAbsDet_ = 0.0;
  int sign_ = 1;
};

// Block LU of an open block-tridiagonal chain with given diagonal blocks
// and -c I between neighbours. Solves, and the inertia and log-determinant
// from the Schur complements (Haynsworth: the inertia of the chain is the
// sum over its Schur blocks).
class HaynsworthChain {
public:
  HaynsworthChain(double c, const std::vector<MatrixXd> &diag, bool spectrum)
      : c_(c) {
    lu_.reserve(diag.size());
    for (size_t k = 0; k < diag.size(); ++k) {
      const long f = diag[k].rows();
      MatrixXd d = 0.5 * (diag[k] + diag[k].transpose());
      if (k > 0) {
        const MatrixXd inv = lu_.back().inverse();
        d -= c * c * 0.5 * (inv + inv.transpose());
      }
      if (spectrum) {
        const ColMajorXd sym = d;
        const Eigen::SelfAdjointEigenSolver<ColMajorXd> es(
            sym, Eigen::EigenvaluesOnly);
        for (long i = 0; i < f; ++i) {
          const double lam = es.eigenvalues()(i);
          if (lam == 0.0) {
            throw std::runtime_error("instanton: singular chain Hessian block");
          }
          if (lam < 0.0) {
            ++negative_;
          }
          logAbsDet_ += std::log(std::abs(lam));
        }
      }
      lu_.emplace_back(ColMajorXd(d));
    }
  }
  double logAbsDet() const { return logAbsDet_; }
  long negative() const { return negative_; }
  std::vector<VectorXd> solve(const std::vector<VectorXd> &b) const {
    const size_t m = b.size();
    std::vector<VectorXd> y(m), x(m);
    y[0] = b[0];
    for (size_t k = 1; k < m; ++k) {
      y[k] = b[k] + c_ * lu_[k - 1].solve(y[k - 1]);
    }
    x[m - 1] = lu_[m - 1].solve(y[m - 1]);
    for (size_t k = m - 1; k-- > 0;) {
      x[k] = lu_[k].solve(y[k] + c_ * x[k + 1]);
    }
    return x;
  }

private:
  double c_;
  std::vector<Eigen::PartialPivLU<ColMajorXd>> lu_;
  double logAbsDet_ = 0.0;
  long negative_ = 0;
};

double dot(const std::vector<VectorXd> &a, const std::vector<VectorXd> &b) {
  double s = 0.0;
  for (size_t k = 0; k < a.size(); ++k) {
    s += a[k].dot(b[k]);
  }
  return s;
}

void scale(std::vector<VectorXd> &a, double f) {
  for (auto &v : a) {
    v *= f;
  }
}

// Point at arc-length fraction f along a polyline.
VectorXd alongPolyline(const std::vector<VectorXd> &pts,
                       const std::vector<double> &cum, double f) {
  const double target = f * cum.back();
  const auto it = std::upper_bound(cum.begin(), cum.end(), target);
  const size_t k = std::clamp<size_t>(
      static_cast<size_t>(std::distance(cum.begin(), it)), 1, pts.size() - 1);
  const double seg = cum[k] - cum[k - 1];
  const double t = seg > 0.0 ? (target - cum[k - 1]) / seg : 0.0;
  return pts[k - 1] + std::clamp(t, 0.0, 1.0) * (pts[k] - pts[k - 1]);
}

struct ActionEval {
  double action = 0.0;
  std::vector<double> v;         // interior beads
  std::vector<VectorXd> grad;    // dS/dq over interior beads
  std::vector<VectorXd> potGrad; // dV/dq over interior beads
};

ActionEval evaluateAction(const std::vector<VectorXd> &interior,
                          const VectorXd &start, const VectorXd &end,
                          double vStart, double vEnd, double dtau,
                          const BatchPotential &potential) {
  ActionEval out;
  potential(interior, out.v, out.potGrad);
  const size_t m = interior.size();
  if (out.v.size() != m || out.potGrad.size() != m) {
    throw std::runtime_error("instanton: potential returned the wrong count");
  }
  auto bead = [&](size_t j) -> const VectorXd & {
    return j == 0 ? start : (j == m + 1 ? end : interior[j - 1]);
  };
  double kinetic = 0.0;
  for (size_t j = 0; j <= m; ++j) {
    kinetic += (bead(j + 1) - bead(j)).squaredNorm();
  }
  double pot = 0.5 * (vStart + vEnd);
  for (double vj : out.v) {
    pot += vj;
  }
  out.action = 0.5 * kinetic / dtau + dtau * pot;
  out.grad.resize(m);
  for (size_t j = 1; j <= m; ++j) {
    out.grad[j - 1] = (2.0 * bead(j) - bead(j - 1) - bead(j + 1)) / dtau +
                      dtau * out.potGrad[j - 1];
  }
  return out;
}

double largestBeadNorm(const std::vector<VectorXd> &g) {
  double m = 0.0;
  for (const auto &v : g) {
    m = std::max(m, v.norm());
  }
  return m;
}

} // namespace

double pathOmega(const MatrixXd &hessStart, const MatrixXd &hessEnd,
                 const VectorXd &start, const VectorXd &end) {
  const VectorXd d = (end - start).normalized();
  const double k = std::max(d.dot(hessStart * d), d.dot(hessEnd * d));
  if (!(k > 0.0)) {
    throw std::invalid_argument(
        "pathOmega: no positive curvature along the path at either minimum");
  }
  return std::sqrt(k);
}

Instanton optimizeInstanton(const VectorXd &start, const VectorXd &end,
                            double betaHbar, std::vector<VectorXd> guess,
                            const BatchPotential &potential,
                            const InstantonOptions &options) {
  const long P = options.beads;
  if (P < 4 || !(betaHbar > 0.0) || start.size() != end.size()) {
    throw std::invalid_argument("optimizeInstanton: need P >= 4, beta hbar > 0 "
                                "and ends of one dimension");
  }
  Instanton inst;
  inst.betaHbar = betaHbar;
  inst.dtau = betaHbar / static_cast<double>(P);
  const double dtau = inst.dtau;

  std::vector<double> vEnds;
  std::vector<VectorXd> gEnds;
  potential({start, end}, vEnds, gEnds);
  if (vEnds.size() != 2) {
    throw std::runtime_error("instanton: potential returned the wrong count");
  }
  inst.asymmetry = vEnds[1] - vEnds[0];

  // Beads along the guess (or the straight line) on a tanh kink centred at
  // beta hbar / 2 whose width follows the harmonic decay of a well.
  if (guess.size() < 2) {
    guess = {start, end};
  }
  std::vector<double> cum(guess.size(), 0.0);
  for (size_t k = 1; k < guess.size(); ++k) {
    cum[k] = cum[k - 1] + (guess[k] - guess[k - 1]).norm();
  }
  if (!(cum.back() > 0.0)) {
    throw std::invalid_argument("optimizeInstanton: the two minima coincide");
  }
  const double width = betaHbar / (2.0 * options.betaHbarOmega);
  std::vector<VectorXd> x(static_cast<size_t>(P - 1));
  for (long j = 1; j < P; ++j) {
    const double tau = static_cast<double>(j) * dtau - 0.5 * betaHbar;
    const double f = 0.5 * (1.0 + std::tanh(tau / width));
    x[static_cast<size_t>(j - 1)] = alongPolyline(guess, cum, f);
  }

  // L-BFGS with a backtracking Armijo line search.
  ActionEval cur =
      evaluateAction(x, start, end, vEnds[0], vEnds[1], dtau, potential);
  std::deque<std::pair<std::vector<VectorXd>, std::vector<VectorXd>>> pairs;
  for (long it = 0; it < options.maxIterations; ++it) {
    inst.iterations = it;
    if (largestBeadNorm(cur.grad) / dtau < options.forceTolerance) {
      inst.converged = true;
      break;
    }
    std::vector<VectorXd> q = cur.grad;
    std::vector<double> alpha(pairs.size());
    for (size_t i = pairs.size(); i-- > 0;) {
      const double rho = 1.0 / dot(pairs[i].second, pairs[i].first);
      alpha[i] = rho * dot(pairs[i].first, q);
      for (size_t k = 0; k < q.size(); ++k) {
        q[k] -= alpha[i] * pairs[i].second[k];
      }
    }
    // Without history, half the inverse spring stiffness 2 / dtau.
    double gamma = dtau / 4.0;
    if (!pairs.empty()) {
      gamma = dot(pairs.back().first, pairs.back().second) /
              dot(pairs.back().second, pairs.back().second);
    }
    scale(q, gamma);
    for (size_t i = 0; i < pairs.size(); ++i) {
      const double rho = 1.0 / dot(pairs[i].second, pairs[i].first);
      const double beta = rho * dot(pairs[i].second, q);
      for (size_t k = 0; k < q.size(); ++k) {
        q[k] += (alpha[i] - beta) * pairs[i].first[k];
      }
    }
    // q is now the inverse-Hessian estimate times the gradient; step -q.
    double slope = -dot(cur.grad, q);
    if (!(slope < 0.0)) {
      pairs.clear();
      q = cur.grad;
      scale(q, dtau / 4.0);
      slope = -dot(cur.grad, q);
    }
    double step = 1.0;
    ActionEval next;
    std::vector<VectorXd> trial(x.size());
    bool accepted = false;
    for (int ls = 0; ls < 30; ++ls) {
      for (size_t k = 0; k < x.size(); ++k) {
        trial[k] = x[k] - step * q[k];
      }
      next = evaluateAction(trial, start, end, vEnds[0], vEnds[1], dtau,
                            potential);
      // Near the minimum the action changes by less than its round-off;
      // there a step that shrinks the gradient is progress too.
      const bool armijo = next.action <= cur.action + 1e-4 * step * slope;
      const bool flat = std::abs(next.action - cur.action) <=
                        1e-13 * std::max(1.0, std::abs(cur.action));
      if (armijo ||
          (flat && largestBeadNorm(next.grad) < largestBeadNorm(cur.grad))) {
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) {
      break;
    }
    std::vector<VectorXd> sk(x.size()), yk(x.size());
    for (size_t k = 0; k < x.size(); ++k) {
      sk[k] = trial[k] - x[k];
      yk[k] = next.grad[k] - cur.grad[k];
    }
    if (dot(sk, yk) > 0.0) {
      pairs.emplace_back(std::move(sk), std::move(yk));
      if (static_cast<long>(pairs.size()) > options.memory) {
        pairs.pop_front();
      }
    }
    x = std::move(trial);
    cur = std::move(next);
  }
  if (!inst.converged &&
      largestBeadNorm(cur.grad) / dtau < options.forceTolerance) {
    inst.converged = true;
  }

  inst.path.reserve(static_cast<size_t>(P + 1));
  inst.path.push_back(start);
  inst.path.insert(inst.path.end(), x.begin(), x.end());
  inst.path.push_back(end);
  inst.energies.reserve(static_cast<size_t>(P + 1));
  inst.energies.push_back(vEnds[0]);
  inst.energies.insert(inst.energies.end(), cur.v.begin(), cur.v.end());
  inst.energies.push_back(vEnds[1]);
  const double sWell = betaHbar * 0.5 * (vEnds[0] + vEnds[1]);
  inst.action = (cur.action - sWell) / kHbar;
  double s0 = 0.0;
  for (long j = 0; j < P; ++j) {
    s0 += (inst.path[static_cast<size_t>(j + 1)] -
           inst.path[static_cast<size_t>(j)])
              .squaredNorm();
  }
  inst.s0 = s0 / dtau;
  inst.symmetricEnough = std::abs(inst.asymmetry) * betaHbar / kHbar < 0.1;
  return inst;
}

void instantonSplitting(Instanton &inst, const BeadHessian &hessian,
                        const MatrixXd &hessStart, const MatrixXd &hessEnd) {
  const long P = static_cast<long>(inst.path.size()) - 1;
  if (P < 4 || !(inst.dtau > 0.0)) {
    throw std::invalid_argument("instantonSplitting: no optimised path");
  }
  const double dtau = inst.dtau;
  const double c = 1.0 / dtau;
  const long n = inst.path.front().size();
  const MatrixXd spring = 2.0 * c * MatrixXd::Identity(n, n);

  std::vector<MatrixXd> diag;
  diag.reserve(static_cast<size_t>(P - 1));
  for (long j = 1; j < P; ++j) {
    const MatrixXd h = hessian(j, inst.path[static_cast<size_t>(j)]);
    if (h.rows() != n || h.cols() != n) {
      throw std::runtime_error("instantonSplitting: bead Hessian size");
    }
    diag.push_back(spring + dtau * 0.5 * (h + h.transpose()));
  }
  // The inertia comes from the Schur blocks (Haynsworth), so a path with
  // two negative modes is not mistaken for a minimum by the determinant's
  // sign.
  const HaynsworthChain chain(c, diag, true);

  auto wellLogDet = [&](const MatrixXd &h) {
    const std::vector<MatrixXd> d(static_cast<size_t>(P - 1),
                                  spring + dtau * 0.5 * (h + h.transpose()));
    const HaynsworthChain well(c, d, true);
    if (well.negative() > 0) {
      throw std::runtime_error(
          "instantonSplitting: a well Hessian is not positive definite");
    }
    return well.logAbsDet();
  };
  const double logDetWell = 0.5 * (wellLogDet(hessStart) + wellLogDet(hessEnd));

  // The zero mode is the kink's translation in imaginary time, along the
  // discrete velocity v; det' J = det J (v^T J^-1 v) for v its eigenvector.
  std::vector<VectorXd> v(static_cast<size_t>(P - 1));
  for (long j = 1; j < P; ++j) {
    v[static_cast<size_t>(j - 1)] = inst.path[static_cast<size_t>(j + 1)] -
                                    inst.path[static_cast<size_t>(j - 1)];
  }
  scale(v, 1.0 / std::sqrt(dot(v, v)));
  const double vJv = dot(v, chain.solve(v));
  // A zero mode just below zero is one negative eigenvalue the prime
  // leaves out; any other negative eigenvalue is a second unstable mode.
  if (chain.negative() - (vJv < 0.0 ? 1 : 0) != 0) {
    throw std::runtime_error(
        "instantonSplitting: the path is not a minimum of the action "
        "(a negative mode besides the kink's translation)");
  }
  inst.zeroMode = 1.0 / vJv;
  const double logDetPrime = chain.logAbsDet() + std::log(std::abs(vJv));

  // Next eigenvalue: inverse iteration orthogonal to v.
  std::vector<VectorXd> w(v.size());
  for (size_t k = 0; k < w.size(); ++k) {
    w[k].resize(n);
    for (long i = 0; i < n; ++i) {
      w[k](i) = std::sin(0.7 * static_cast<double>(k) +
                         1.3 * static_cast<double>(i) + 0.1);
    }
  }
  double lambda1 = 0.0;
  for (int it = 0; it < 40; ++it) {
    const double proj = dot(v, w);
    for (size_t k = 0; k < w.size(); ++k) {
      w[k] -= proj * v[k];
    }
    scale(w, 1.0 / std::sqrt(dot(w, w)));
    std::vector<VectorXd> z = chain.solve(w);
    lambda1 = 1.0 / dot(w, z);
    w = std::move(z);
  }
  inst.modeSeparation = std::abs(lambda1 / inst.zeroMode);

  inst.delta0 = 2.0 * kHbar *
                std::sqrt(inst.s0 / (2.0 * std::numbers::pi * kHbar * dtau)) *
                std::exp(0.5 * (logDetWell - logDetPrime) - inst.action);
}

double crossoverTemperature(const MatrixXd &hessSaddle) {
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(
      0.5 * (hessSaddle + hessSaddle.transpose()));
  const double lambda = es.eigenvalues()(0);
  if (!(lambda < 0.0)) {
    throw std::invalid_argument(
        "crossoverTemperature: the saddle Hessian has no negative eigenvalue");
  }
  return kHbar * std::sqrt(-lambda) / (2.0 * std::numbers::pi * kBoltzmann);
}

namespace {

struct RingEval {
  double u = 0.0;             // U_N
  std::vector<double> v;      // V at each bead
  std::vector<VectorXd> grad; // dU_N/dq at each bead
};

RingEval evaluateRing(const std::vector<VectorXd> &x, double c,
                      const BatchPotential &potential,
                      double energyShift = 0.0) {
  RingEval out;
  std::vector<VectorXd> gv;
  potential(x, out.v, gv);
  for (double &v : out.v) {
    v -= energyShift;
  }
  const size_t n = x.size();
  if (out.v.size() != n || gv.size() != n) {
    throw std::runtime_error(
        "rate instanton: potential returned the wrong count");
  }
  out.grad.resize(n);
  double spring = 0.0;
  for (size_t j = 0; j < n; ++j) {
    const VectorXd &prev = x[(j + n - 1) % n];
    const VectorXd &next = x[(j + 1) % n];
    out.grad[j] = gv[j] + c * (2.0 * x[j] - prev - next);
    spring += (next - x[j]).squaredNorm();
    out.u += out.v[j];
  }
  out.u += 0.5 * c * spring;
  return out;
}

// Lowest eigenpair of the ring Hessian by Lanczos on finite-difference
// products, started from `start`; full reorthogonalisation.
double lowestMode(const std::vector<VectorXd> &x, const RingEval &here,
                  double c, const BatchPotential &potential,
                  std::vector<VectorXd> &mode, long steps, double eps,
                  bool mirror = false) {
  const size_t n = x.size();
  // The thermal instanton retraces, so the unstable mode is even. Probes
  // and products stay on that mirror and the potential call sees one half.
  const bool reflect = mirror && n % 2 == 0;
  auto snap = [&](std::vector<VectorXd> &q) {
    if (!reflect) {
      return;
    }
    const long m = static_cast<long>(n) / 2;
    for (long j = 1; j < m; ++j) {
      const size_t a = static_cast<size_t>(j);
      const size_t b = n - a;
      const VectorXd mid = 0.5 * (q[a] + q[b]);
      q[a] = mid;
      q[b] = mid;
    }
  };
  auto hv = [&](const std::vector<VectorXd> &u) {
    std::vector<VectorXd> xp(n);
    for (size_t j = 0; j < n; ++j) {
      xp[j] = x[j] + eps * u[j];
    }
    snap(xp);
    const RingEval e = evaluateRing(xp, c, potential);
    std::vector<VectorXd> out(n);
    for (size_t j = 0; j < n; ++j) {
      out[j] = (e.grad[j] - here.grad[j]) / eps;
    }
    snap(out);
    return out;
  };
  std::vector<std::vector<VectorXd>> basis;
  std::vector<double> alpha, beta;
  std::vector<VectorXd> q = mode;
  snap(q);
  scale(q, 1.0 / std::sqrt(dot(q, q)));
  for (long k = 0; k < steps; ++k) {
    basis.push_back(q);
    std::vector<VectorXd> w = hv(q);
    const double a = dot(w, q);
    alpha.push_back(a);
    for (const auto &b : basis) {
      const double p = dot(w, b);
      for (size_t j = 0; j < n; ++j) {
        w[j] -= p * b[j];
      }
    }
    snap(w);
    const double bnorm = std::sqrt(dot(w, w));
    if (!(bnorm > 1e-12) || k + 1 == steps) {
      break;
    }
    beta.push_back(bnorm);
    scale(w, 1.0 / bnorm);
    q = std::move(w);
  }
  const long m = static_cast<long>(alpha.size());
  MatrixXd t = MatrixXd::Zero(m, m);
  for (long i = 0; i < m; ++i) {
    t(i, i) = alpha[static_cast<size_t>(i)];
    if (i + 1 < m) {
      t(i, i + 1) = t(i + 1, i) = beta[static_cast<size_t>(i)];
    }
  }
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(t);
  const VectorXd y = es.eigenvectors().col(0);
  std::vector<VectorXd> ritz(n);
  for (size_t j = 0; j < n; ++j) {
    ritz[j] = VectorXd::Zero(x[j].size());
  }
  for (long i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      ritz[j] += y(i) * basis[static_cast<size_t>(i)][j];
    }
  }
  scale(ritz, 1.0 / std::sqrt(dot(ritz, ritz)));
  snap(ritz);
  const double rnorm = std::sqrt(dot(ritz, ritz));
  if (rnorm > 0.0) {
    scale(ritz, 1.0 / rnorm);
  }
  mode = std::move(ritz);
  return es.eigenvalues()(0);
}

// Lowest curvature of the ring Hessian over vectors odd under j -> N-j and
// orthogonal to the cycle mode x[j+1] - x[j-1], at a mirror-symmetric ring
// of even N. A search confined to the even sector cannot see these. The
// products move the ring off the mirror, so `potential` evaluates every bead.
// Returns +inf when no odd direction is left to probe.
double lowestOddMode(const std::vector<VectorXd> &x, const RingEval &here,
                     double c, const BatchPotential &potential,
                     std::vector<VectorXd> &mode, long steps, double eps) {
  const size_t n = x.size();
  const size_t m = n / 2;
  std::vector<VectorXd> cycle(n);
  for (size_t j = 0; j < n; ++j) {
    cycle[j] = x[(j + 1) % n] - x[(j + n - 1) % n];
  }
  auto odd = [&](std::vector<VectorXd> &q) {
    q[0].setZero();
    q[m].setZero();
    for (size_t a = 1; a < m; ++a) {
      const VectorXd half = 0.5 * (q[a] - q[n - a]);
      q[a] = half;
      q[n - a] = -half;
    }
  };
  odd(cycle);
  const double cnorm = std::sqrt(dot(cycle, cycle));
  if (cnorm > 0.0) {
    scale(cycle, 1.0 / cnorm);
  }
  auto project = [&](std::vector<VectorXd> &q) {
    odd(q);
    if (cnorm > 0.0) {
      const double p = dot(q, cycle);
      for (size_t j = 0; j < n; ++j) {
        q[j] -= p * cycle[j];
      }
    }
  };
  auto hv = [&](const std::vector<VectorXd> &u) {
    std::vector<VectorXd> xp(n);
    for (size_t j = 0; j < n; ++j) {
      xp[j] = x[j] + eps * u[j];
    }
    const RingEval e = evaluateRing(xp, c, potential);
    std::vector<VectorXd> out(n);
    for (size_t j = 0; j < n; ++j) {
      out[j] = (e.grad[j] - here.grad[j]) / eps;
    }
    project(out);
    return out;
  };
  // Start from the bead displacements with opposite signs on the two
  // halves: two copies of the turning region moving against each other.
  VectorXd mean = VectorXd::Zero(x[0].size());
  for (const auto &b : x) {
    mean += b;
  }
  mean /= static_cast<double>(n);
  std::vector<VectorXd> q(n);
  for (size_t j = 0; j < n; ++j) {
    q[j] = (j < m ? 1.0 : -1.0) * (x[j] - mean);
  }
  project(q);
  const double qnorm = std::sqrt(dot(q, q));
  if (!(qnorm > 1e-12)) {
    return std::numeric_limits<double>::infinity();
  }
  scale(q, 1.0 / qnorm);
  std::vector<std::vector<VectorXd>> basis;
  std::vector<double> alpha, beta;
  for (long k = 0; k < steps; ++k) {
    basis.push_back(q);
    std::vector<VectorXd> w = hv(q);
    alpha.push_back(dot(w, q));
    for (const auto &b : basis) {
      const double p = dot(w, b);
      for (size_t j = 0; j < n; ++j) {
        w[j] -= p * b[j];
      }
    }
    project(w);
    const double bnorm = std::sqrt(dot(w, w));
    if (!(bnorm > 1e-12) || k + 1 == steps) {
      break;
    }
    beta.push_back(bnorm);
    scale(w, 1.0 / bnorm);
    q = std::move(w);
  }
  const long dim = static_cast<long>(alpha.size());
  MatrixXd t = MatrixXd::Zero(dim, dim);
  for (long i = 0; i < dim; ++i) {
    t(i, i) = alpha[static_cast<size_t>(i)];
    if (i + 1 < dim) {
      t(i, i + 1) = t(i + 1, i) = beta[static_cast<size_t>(i)];
    }
  }
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(t);
  const VectorXd y = es.eigenvectors().col(0);
  mode.assign(n, VectorXd::Zero(x[0].size()));
  for (long i = 0; i < dim; ++i) {
    for (size_t j = 0; j < n; ++j) {
      mode[j] += y(i) * basis[static_cast<size_t>(i)][j];
    }
  }
  project(mode);
  scale(mode, 1.0 / std::sqrt(dot(mode, mode)));
  return es.eigenvalues()(0);
}

// Distance along +dir (sign = 1) or -dir (sign = -1) from the saddle at which
// V has dropped by `drop`, from a scan in steps of h; the lowest point found
// when V never drops that far.
double turningDistance(const VectorXd &saddle, const VectorXd &dir,
                       double vSaddle, double drop, double sign, double h,
                       long maxPoints, const BatchPotential &potential) {
  std::vector<VectorXd> pts;
  for (long k = 1; k <= maxPoints; ++k) {
    pts.push_back(saddle + sign * h * static_cast<double>(k) * dir);
  }
  std::vector<double> v;
  std::vector<VectorXd> g;
  potential(pts, v, g);
  double prevV = vSaddle, prevS = 0.0, best = 0.0, bestV = vSaddle;
  for (long k = 0; k < maxPoints; ++k) {
    const double sk = h * static_cast<double>(k + 1);
    const double vk = v[static_cast<size_t>(k)];
    if (vk <= vSaddle - drop) {
      const double t = (prevV - (vSaddle - drop)) / (prevV - vk);
      return prevS + t * (sk - prevS);
    }
    if (vk < bestV) {
      bestV = vk;
      best = sk;
    }
    if (vk > prevV + 1e-12 && k > 0) {
      break; // past the minimum on this side
    }
    prevV = vk;
    prevS = sk;
  }
  return best;
}

double sideDrop(const VectorXd &saddle, const VectorXd &dir, double vSaddle,
                double sign, double h, long maxPoints,
                const BatchPotential &potential) {
  std::vector<VectorXd> pts;
  for (long k = 1; k <= maxPoints; ++k) {
    pts.push_back(saddle + sign * h * static_cast<double>(k) * dir);
  }
  std::vector<double> v;
  std::vector<VectorXd> g;
  potential(pts, v, g);
  double lowest = vSaddle;
  for (long k = 0; k < maxPoints; ++k) {
    const double vk = v[static_cast<size_t>(k)];
    if (vk > lowest + 1e-12 && k > 0 && lowest < vSaddle) {
      return vSaddle - lowest; // found this side's minimum
    }
    lowest = std::min(lowest, vk);
  }
  return vSaddle - lowest; // still falling: the drop to the last point
}

struct RingMode {
  double theta = 0.0;
  double residual = 0.0;
  std::vector<VectorXd> vector;
};

// Lowest Ritz pairs of a symmetric ring operator, full reorthogonalisation.
// `steps` at the dimension is the whole spectrum.
std::vector<RingMode> lowestRingModes(
    const std::function<std::vector<VectorXd>(const std::vector<VectorXd> &)>
        &apply,
    std::vector<VectorXd> start, long steps) {
  const size_t n = start.size();
  const long f = start.empty() ? 0 : start.front().size();
  const double n0 = std::sqrt(dot(start, start));
  if (!(n0 > 0.0) || f < 1) {
    throw std::runtime_error("instantonRate: Lanczos was given a zero vector");
  }
  scale(start, 1.0 / n0);
  std::vector<std::vector<VectorXd>> basis;
  std::vector<double> alpha;
  std::vector<double> beta;
  std::vector<VectorXd> q = std::move(start);
  for (long k = 0; k < steps; ++k) {
    basis.push_back(q);
    std::vector<VectorXd> w = apply(q);
    alpha.push_back(dot(w, q));
    for (int pass = 0; pass < 2; ++pass) {
      for (const auto &b : basis) {
        const double p = dot(w, b);
        for (size_t j = 0; j < n; ++j) {
          w[j] -= p * b[j];
        }
      }
    }
    const double bnorm = std::sqrt(dot(w, w));
    if (!(bnorm > 1e-14) || k + 1 == steps) {
      break;
    }
    beta.push_back(bnorm);
    scale(w, 1.0 / bnorm);
    q = std::move(w);
  }
  const long m = static_cast<long>(alpha.size());
  MatrixXd tridiag = MatrixXd::Zero(m, m);
  for (long i = 0; i < m; ++i) {
    tridiag(i, i) = alpha[static_cast<size_t>(i)];
    if (i + 1 < m) {
      tridiag(i, i + 1) = tridiag(i + 1, i) = beta[static_cast<size_t>(i)];
    }
  }
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(tridiag);
  std::vector<RingMode> modes;
  modes.reserve(static_cast<size_t>(m));
  for (long i = 0; i < m; ++i) {
    RingMode mode;
    mode.theta = es.eigenvalues()(i);
    mode.vector.assign(n, VectorXd::Zero(f));
    const VectorXd y = es.eigenvectors().col(i);
    for (long s = 0; s < m; ++s) {
      for (size_t j = 0; j < n; ++j) {
        mode.vector[j] += y(s) * basis[static_cast<size_t>(s)][j];
      }
    }
    scale(mode.vector, 1.0 / std::sqrt(dot(mode.vector, mode.vector)));
    const std::vector<VectorXd> applied = apply(mode.vector);
    double residual = 0.0;
    for (size_t j = 0; j < n; ++j) {
      residual += (applied[j] - mode.theta * mode.vector[j]).squaredNorm();
    }
    mode.residual = std::sqrt(residual);
    modes.push_back(std::move(mode));
  }
  return modes;
}

// Factored cyclic ring Hessian: the open chain plus the corner coupling
// M between bead 0 and bead N-1.
// det(H_open + P M P^T) = det(H_open) det(M) det(M^{-1} + P^T H_open^{-1} P),
// and |det M| = c^{2f}. A singular ring makes the last factor zero.
struct CyclicFactor {
  BlockChain open;
  Eigen::PartialPivLU<ColMajorXd> cornerLu;
  double c;
  long f;
  long n;
  double logAbs;
  bool singular;

  CyclicFactor(double cIn, const std::vector<MatrixXd> &diag)
      : open(cIn, diag),
        c(cIn),
        f(diag.front().rows()),
        n(static_cast<long>(diag.size())),
        logAbs(0.0),
        singular(false) {
    const MatrixXd eye = MatrixXd::Identity(f, f);
    const MatrixXd zero = MatrixXd::Zero(f, f);
    std::vector<MatrixXd> rhs(static_cast<size_t>(n), zero);
    rhs.front() = eye;
    const std::vector<MatrixXd> fromFirst = open.solve(rhs);
    rhs.front() = zero;
    rhs.back() = eye;
    const std::vector<MatrixXd> fromLast = open.solve(rhs);
    MatrixXd corner(2 * f, 2 * f);
    corner.topLeftCorner(f, f) = fromFirst.front();
    corner.bottomLeftCorner(f, f) = fromFirst.back();
    corner.topRightCorner(f, f) = fromLast.front();
    corner.bottomRightCorner(f, f) = fromLast.back();
    corner.topRightCorner(f, f) -= eye / c;
    corner.bottomLeftCorner(f, f) -= eye / c;
    cornerLu.compute(ColMajorXd(corner));
    logAbs = open.logAbsDet() + 2.0 * static_cast<double>(f) * std::log(c);
    const MatrixXd &upper = cornerLu.matrixLU();
    for (long i = 0; i < upper.rows(); ++i) {
      const double pivot = upper(i, i);
      if (pivot == 0.0) {
        singular = true;
        logAbs = -std::numeric_limits<double>::infinity();
        return;
      }
      logAbs += std::log(std::abs(pivot));
    }
  }

  std::vector<VectorXd> solve(const std::vector<VectorXd> &rhs) const {
    if (singular || static_cast<long>(rhs.size()) != n) {
      throw std::runtime_error("cyclic ring: singular");
    }
    std::vector<VectorXd> y = open.solve(rhs);
    VectorXd g(2 * f);
    g.head(f) = y.front();
    g.tail(f) = y.back();
    const VectorXd z = cornerLu.solve(g);
    std::vector<VectorXd> bump(static_cast<size_t>(n), VectorXd::Zero(f));
    bump.front() = z.head(f);
    bump.back() = z.tail(f);
    const std::vector<VectorXd> corr = open.solve(bump);
    for (long j = 0; j < n; ++j) {
      y[static_cast<size_t>(j)] -= corr[static_cast<size_t>(j)];
    }
    return y;
  }
};

// J = T + G K G^T over the chain T: the closure blocks (-c I between the
// last bead and the first) when `closed`, and symmetric rank-one terms
// kappa_i u_i u_i^T. Woodbury gives J^{-1} b, the determinant lemma
// ln|det J| and Haynsworth the inertia, all O(N f^3), so the N f by N f
// matrix is never formed.
class WoodburyRing {
public:
  WoodburyRing(double c, const std::vector<MatrixXd> &diag, bool closed,
               const std::vector<std::vector<VectorXd>> &extras,
               const std::vector<double> &kappas, bool spectrum)
      : c_(c),
        n_(static_cast<long>(diag.size())),
        f_(diag.front().rows()),
        closed_(closed),
        chain_(c, diag, spectrum),
        extras_(extras) {
    const long base = closed ? 2 * f_ : 0;
    const long m = base + static_cast<long>(extras.size());
    kinv_ = MatrixXd::Zero(m, m);
    MatrixXd k = MatrixXd::Zero(m, m);
    if (closed) {
      kinv_.block(0, f_, f_, f_) = -MatrixXd::Identity(f_, f_) / c;
      kinv_.block(f_, 0, f_, f_) = -MatrixXd::Identity(f_, f_) / c;
      k.block(0, f_, f_, f_) = -c * MatrixXd::Identity(f_, f_);
      k.block(f_, 0, f_, f_) = -c * MatrixXd::Identity(f_, f_);
    }
    for (size_t i = 0; i < extras.size(); ++i) {
      const long r = base + static_cast<long>(i);
      kinv_(r, r) = 1.0 / kappas[i];
      k(r, r) = kappas[i];
    }
    if (m == 0) {
      ok_ = true;
      logAbsDet_ = chain_.logAbsDet();
      negative_ = chain_.negative();
      return;
    }
    gtg_ = MatrixXd::Zero(m, m);
    for (long col = 0; col < m; ++col) {
      gtg_.col(col) = pieces(chain_.solve(column(col)));
    }
    woodbury_.compute(ColMajorXd(kinv_ + gtg_));
    ok_ = gtg_.array().isFinite().all();
    if (spectrum) {
      const Eigen::PartialPivLU<ColMajorXd> lu(
          ColMajorXd(MatrixXd::Identity(m, m) + k * gtg_));
      double logDet = 0.0;
      for (long i = 0; i < m; ++i) {
        const double u = lu.matrixLU()(i, i);
        if (u == 0.0) {
          ok_ = false;
          logAbsDet_ = -std::numeric_limits<double>::infinity();
          return;
        }
        logDet += std::log(std::abs(u));
      }
      logAbsDet_ = chain_.logAbsDet() + logDet;
      // neg(J) = neg(T) + neg(S) - neg(-K^{-1}), S = -K^{-1} - G^T T^{-1} G.
      const ColMajorXd sMat = -kinv_ - gtg_;
      const ColMajorXd sSym = 0.5 * (sMat + sMat.transpose());
      const Eigen::SelfAdjointEigenSolver<ColMajorXd> es(
          sSym, Eigen::EigenvaluesOnly);
      const ColMajorXd kn = -kinv_;
      const Eigen::SelfAdjointEigenSolver<ColMajorXd> ek(
          kn, Eigen::EigenvaluesOnly);
      negative_ = chain_.negative() + (es.eigenvalues().array() < 0.0).count() -
                  (ek.eigenvalues().array() < 0.0).count();
    }
  }
  bool ok() const { return ok_; }
  std::vector<VectorXd> solve(const std::vector<VectorXd> &b) const {
    std::vector<VectorXd> tb = chain_.solve(b);
    if (kinv_.rows() == 0) {
      return tb;
    }
    const VectorXd y = woodbury_.solve(pieces(tb));
    const std::vector<VectorXd> gy = chain_.solve(expand(y));
    for (size_t j = 0; j < tb.size(); ++j) {
      tb[j] -= gy[j];
    }
    return tb;
  }
  double logAbsDet() const { return logAbsDet_; }
  long negative() const { return negative_; }

private:
  std::vector<VectorXd> column(long col) const {
    const long base = closed_ ? 2 * f_ : 0;
    if (col < base) {
      std::vector<VectorXd> g(static_cast<size_t>(n_), VectorXd::Zero(f_));
      g[static_cast<size_t>(col < f_ ? 0 : n_ - 1)](col % f_) = 1.0;
      return g;
    }
    return extras_[static_cast<size_t>(col - base)];
  }
  VectorXd pieces(const std::vector<VectorXd> &x) const {
    const long base = closed_ ? 2 * f_ : 0;
    VectorXd out(base + static_cast<long>(extras_.size()));
    if (closed_) {
      out.head(f_) = x[0];
      out.segment(f_, f_) = x[static_cast<size_t>(n_ - 1)];
    }
    for (size_t i = 0; i < extras_.size(); ++i) {
      out(base + static_cast<long>(i)) = dot(extras_[i], x);
    }
    return out;
  }
  std::vector<VectorXd> expand(const VectorXd &y) const {
    const long base = closed_ ? 2 * f_ : 0;
    std::vector<VectorXd> g(static_cast<size_t>(n_), VectorXd::Zero(f_));
    if (closed_) {
      g[0] += y.head(f_);
      g[static_cast<size_t>(n_ - 1)] += y.segment(f_, f_);
    }
    for (size_t i = 0; i < extras_.size(); ++i) {
      const double w = y(base + static_cast<long>(i));
      for (long j = 0; j < n_; ++j) {
        g[static_cast<size_t>(j)] += w * extras_[i][static_cast<size_t>(j)];
      }
    }
    return g;
  }

  double c_;
  long n_, f_;
  bool closed_;
  HaynsworthChain chain_;
  std::vector<std::vector<VectorXd>> extras_;
  MatrixXd kinv_, gtg_;
  Eigen::PartialPivLU<ColMajorXd> woodbury_;
  bool ok_ = false;
  double logAbsDet_ = 0.0;
  long negative_ = 0;
};

// Diagonal blocks of the closed ring (H + 2 c I) or of the half chain from
// one turning point to the other (end blocks H / 2 + c I).
std::vector<MatrixXd> ringDiagonal(const std::vector<MatrixXd> &physical,
                                   double c, bool half) {
  const long beads = static_cast<long>(physical.size());
  const long f = physical.front().rows();
  std::vector<MatrixXd> diag(static_cast<size_t>(beads));
  for (long j = 0; j < beads; ++j) {
    const bool end = half && (j == 0 || j + 1 == beads);
    MatrixXd block = 0.5 * (physical[static_cast<size_t>(j)] +
                            physical[static_cast<size_t>(j)].transpose());
    if (end) {
      block *= 0.5;
    }
    block += (end ? c : 2.0 * c) * MatrixXd::Identity(f, f);
    diag[static_cast<size_t>(j)] = block;
  }
  return diag;
}

// (diag - c neighbours) v for the closed ring or the open half chain.
std::vector<VectorXd> applyDiagonal(const std::vector<MatrixXd> &diag, double c,
                                    bool closed,
                                    const std::vector<VectorXd> &v) {
  const size_t n = v.size();
  std::vector<VectorXd> out(n);
  for (size_t j = 0; j < n; ++j) {
    out[j] = diag[j] * v[j];
    if (j > 0 || closed) {
      out[j] -= c * v[(j + n - 1) % n];
    }
    if (j + 1 < n || closed) {
      out[j] -= c * v[(j + 1) % n];
    }
  }
  return out;
}

} // namespace

namespace {

void requireCyclicBlocks(double c, const std::vector<MatrixXd> &diag,
                         const char *what) {
  if (!(c > 0.0) || diag.empty()) {
    throw std::invalid_argument(std::string(what) +
                                ": need a positive spring constant and blocks");
  }
  const long f = diag.front().rows();
  const long n = static_cast<long>(diag.size());
  if (f < 1 || n < 2) {
    throw std::invalid_argument(std::string(what) +
                                ": need at least two beads and one coordinate");
  }
  for (const auto &block : diag) {
    if (block.rows() != f || block.cols() != f) {
      throw std::invalid_argument(std::string(what) +
                                  ": the blocks differ in size");
    }
  }
}

} // namespace

RingSpectrum ringSpectrum(const std::vector<MatrixXd> &beadHessians, double c,
                          const std::vector<VectorXd> &tau) {
  if (beadHessians.empty() || beadHessians.size() != tau.size()) {
    throw std::invalid_argument(
        "ringSpectrum: N bead Hessians and N tau blocks");
  }
  const std::vector<MatrixXd> diag = ringDiagonal(beadHessians, c, false);
  const WoodburyRing ring(c, diag, true, {tau}, {1.0}, true);
  if (!ring.ok()) {
    throw std::runtime_error("ringSpectrum: singular ring");
  }
  RingSpectrum out;
  out.logDetPrime = ring.logAbsDet();
  out.negativeModes = ring.negative();
  out.zeroEigenvalue = dot(tau, applyDiagonal(diag, c, true, tau));
  return out;
}

double cyclicRingLogAbsDet(double c, const std::vector<MatrixXd> &diag) {
  requireCyclicBlocks(c, diag, "cyclicRingLogAbsDet");
  return CyclicFactor(c, diag).logAbs;
}

std::vector<VectorXd> cyclicRingSolve(double c,
                                      const std::vector<MatrixXd> &diag,
                                      const std::vector<VectorXd> &rhs) {
  requireCyclicBlocks(c, diag, "cyclicRingSolve");
  if (static_cast<long>(rhs.size()) != static_cast<long>(diag.size())) {
    throw std::invalid_argument(
        "cyclicRingSolve: one right-hand side per bead");
  }
  for (const auto &row : rhs) {
    if (row.size() != diag.front().rows()) {
      throw std::invalid_argument(
          "cyclicRingSolve: the right-hand side does not match the blocks");
    }
  }
  return CyclicFactor(c, diag).solve(rhs);
}

namespace {

// Central difference of dV/dq. One batch of 2 f displaced beads.
MatrixXd fdPhysicalHessian(const VectorXd &q, const BatchPotential &potential,
                           double eps) {
  const long f = q.size();
  std::vector<VectorXd> pts;
  pts.reserve(static_cast<size_t>(2 * f));
  for (long a = 0; a < f; ++a) {
    VectorXd qp = q;
    VectorXd qm = q;
    qp[a] += eps;
    qm[a] -= eps;
    pts.push_back(std::move(qp));
    pts.push_back(std::move(qm));
  }
  std::vector<double> v;
  std::vector<VectorXd> g;
  potential(pts, v, g);
  if (static_cast<long>(g.size()) != 2 * f) {
    throw std::runtime_error(
        "rate instanton: the Hessian sample returned the wrong count");
  }
  MatrixXd h(f, f);
  for (long a = 0; a < f; ++a) {
    h.col(a) =
        (g[static_cast<size_t>(2 * a)] - g[static_cast<size_t>(2 * a + 1)]) /
        (2.0 * eps);
  }
  return (0.5 * (h + h.transpose())).eval();
}

// Bofill mix of Powell's symmetric Broyden and SR1. H is d2V/dq2.
void bofillUpdate(MatrixXd &h, const VectorXd &dq, const VectorXd &dg) {
  const double dq2 = dq.squaredNorm();
  if (!(dq2 > 1e-24)) {
    return;
  }
  const VectorXd r = dg - h * dq;
  const double rq = r.dot(dq);
  const double r2 = r.squaredNorm();
  const double phi = r2 * dq2 > 1e-30 ? (rq * rq) / (r2 * dq2) : 0.0;
  const double inv = 1.0 / dq2;
  h.noalias() += (1.0 - phi) * inv *
                 (r * dq.transpose() + dq * r.transpose() -
                  (rq * inv) * (dq * dq.transpose()));
  if (std::abs(rq) > 1e-12 * std::sqrt(std::max(0.0, r2 * dq2))) {
    h.noalias() += (phi / rq) * (r * r.transpose());
  }
  h = (0.5 * (h + h.transpose())).eval();
}

VectorXd packBeads(const std::vector<VectorXd> &x) {
  const long f = x.front().size();
  VectorXd flat(static_cast<long>(x.size()) * f);
  for (long j = 0; j < static_cast<long>(x.size()); ++j) {
    flat.segment(j * f, f) = x[static_cast<size_t>(j)];
  }
  return flat;
}

void addPacked(std::vector<VectorXd> &x, const VectorXd &step) {
  const long f = x.front().size();
  for (long j = 0; j < static_cast<long>(x.size()); ++j) {
    x[static_cast<size_t>(j)] += step.segment(j * f, f);
  }
}

double packedBeadNorm(const VectorXd &step, long f) {
  double big = 0.0;
  const long n = step.size() / f;
  for (long j = 0; j < n; ++j) {
    big = std::max(big, step.segment(j * f, f).norm());
  }
  return big;
}

bool finiteBeads(const std::vector<VectorXd> &x) {
  for (const auto &q : x) {
    if (!q.array().isFinite().all()) {
      return false;
    }
  }
  return true;
}

template <typename Derived>
bool overlapsTau(const Eigen::MatrixBase<Derived> &mode, const VectorXd &tau) {
  return tau.size() == mode.size() && std::abs(mode.dot(tau)) > 0.5;
}

// Normalised bead velocity q_{j+1} - q_{j-1}. Empty when the beads coincide.
VectorXd timeTranslation(const std::vector<VectorXd> &x) {
  const long n = static_cast<long>(x.size());
  const long f = x.front().size();
  VectorXd tau(n * f);
  for (long j = 0; j < n; ++j) {
    const size_t prev = static_cast<size_t>((j + n - 1) % n);
    const size_t next = static_cast<size_t>((j + 1) % n);
    tau.segment(j * f, f) = 0.5 * (x[next] - x[prev]);
  }
  const double nrm = tau.norm();
  if (!(nrm > 0.0)) {
    return VectorXd();
  }
  tau /= nrm;
  return tau;
}

std::vector<VectorXd> physicalGradient(const std::vector<VectorXd> &x,
                                       const std::vector<VectorXd> &ringGrad,
                                       double c) {
  const size_t n = x.size();
  std::vector<VectorXd> g(n);
  for (size_t j = 0; j < n; ++j) {
    const size_t prev = (j + n - 1) % n;
    const size_t next = (j + 1) % n;
    g[j] = ringGrad[j] - c * (2.0 * x[j] - x[prev] - x[next]);
  }
  return g;
}

struct Climb {
  long index = -1;
  double curvature = 0.0;
  long negative = 0;
};

// Cosine between the turning points, opened by (1 - T / Tc) of the lower
// barrier. The far turning point sits on bead 0.
std::vector<VectorXd> cosineSeed(const VectorXd &saddle, const VectorXd &dir,
                                 double lambda0, double temperature,
                                 double crossover, long nBeads,
                                 const BatchPotential &potential) {
  std::vector<double> v0;
  std::vector<VectorXd> g0;
  potential({saddle}, v0, g0);
  const double vS = v0.at(0);
  const double h = 0.25 * std::sqrt(2.0 * kBoltzmann * crossover / -lambda0);
  const long pts = 200;
  const double dPlus = sideDrop(saddle, dir, vS, 1.0, h, pts, potential);
  const double dMinus = sideDrop(saddle, dir, vS, -1.0, h, pts, potential);
  const double dMin = std::min(dPlus, dMinus);
  const double drop = (1.0 - temperature / crossover) *
                      (dMin > 0.0 ? dMin : kBoltzmann * crossover);
  const double sPlus =
      turningDistance(saddle, dir, vS, drop, 1.0, h, pts, potential);
  const double sMinus =
      turningDistance(saddle, dir, vS, drop, -1.0, h, pts, potential);
  std::vector<VectorXd> guess(static_cast<size_t>(nBeads));
  for (long j = 0; j < nBeads; ++j) {
    const double ct = std::cos(2.0 * std::numbers::pi * static_cast<double>(j) /
                               static_cast<double>(nBeads));
    guess[static_cast<size_t>(j)] =
        saddle + dir * (ct >= 0.0 ? sPlus * ct : sMinus * ct);
  }
  return guess;
}

// Index-1 Newton step through the block chain: a negative climb eigenvalue
// stays and a raw Newton step climbs it; a positive one is flipped, as is
// every other negative Ritz value off the cycle; a tiny one is parked at a
// spring-sized curvature; the cycle itself is held with a spring-sized
// curvature and its component removed from the step. Each flip is a
// rank-one term in the Woodbury correction, so the solve stays O(N f^3).
VectorXd chainIndexOneStep(const std::vector<RingMode> &ritz,
                           const Climb &climb,
                           const std::vector<MatrixXd> &diag, double spring,
                           bool closed, const std::vector<VectorXd> &grad,
                           const VectorXd &tau,
                           const std::vector<std::vector<VectorXd>> &nullRing) {
  if (climb.index < 0 || ritz.empty() || grad.empty()) {
    return VectorXd();
  }
  const long f = grad.front().size();
  const long dim = static_cast<long>(grad.size()) * f;
  const double cut = -1e-8 * std::max(1.0, spring);
  const double tiny = 1e-8 * std::max(1.0, spring);
  const double parked = std::max(1.0, spring);
  std::vector<std::vector<VectorXd>> extras;
  std::vector<double> kappas;
  std::vector<VectorXd> tauRing;
  if (tau.size() == dim) {
    tauRing.assign(grad.size(), VectorXd::Zero(f));
    for (size_t j = 0; j < grad.size(); ++j) {
      tauRing[j] = tau.segment(static_cast<long>(j) * f, f);
    }
    // On a discrete ring the time shift has a small curvature of its own,
    // and the stationary ring is where Newton on that curvature leads; the
    // lift is only for a cycle too flat to solve with.
    double cycleCurvature = 0.0;
    for (const auto &m : ritz) {
      if (std::abs(dot(m.vector, tauRing)) > 0.5) {
        cycleCurvature = m.theta;
        break;
      }
    }
    if (std::abs(cycleCurvature) <= 1e-6 * std::max(1.0, spring)) {
      extras.push_back(tauRing);
      kappas.push_back(spring);
    }
  }
  // The rigid ring motions are lifted like the cycle.
  auto onNull = [&](const std::vector<VectorXd> &m) {
    for (const auto &r : nullRing) {
      if (std::abs(dot(m, r)) > 0.5) {
        return true;
      }
    }
    return false;
  };
  for (const auto &r : nullRing) {
    extras.push_back(r);
    kappas.push_back(spring);
  }
  for (size_t i = 0; i < ritz.size(); ++i) {
    if (!tauRing.empty() && std::abs(dot(ritz[i].vector, tauRing)) > 0.5) {
      continue;
    }
    if (onNull(ritz[i].vector)) {
      continue;
    }
    const double li = ritz[i].theta;
    const bool isClimb = static_cast<long>(i) == climb.index;
    if ((isClimb && li > 0.0) || (!isClimb && li < cut)) {
      extras.push_back(ritz[i].vector);
      kappas.push_back(-2.0 * li);
    } else if (!isClimb && std::abs(li) <= tiny) {
      extras.push_back(ritz[i].vector);
      kappas.push_back(parked - li);
    }
  }
  std::vector<VectorXd> rhs = grad;
  scale(rhs, -1.0);
  // The operator the step solves with: the ring plus every rank-one term.
  auto applyShifted = [&](const std::vector<VectorXd> &x) {
    std::vector<VectorXd> out = applyDiagonal(diag, spring, closed, x);
    for (size_t i = 0; i < extras.size(); ++i) {
      const double w = kappas[i] * dot(extras[i], x);
      for (size_t j = 0; j < out.size(); ++j) {
        out[j] += w * extras[i][j];
      }
    }
    return out;
  };
  const double rhsNorm = std::sqrt(dot(rhs, rhs));
  std::vector<VectorXd> stepRing;
  bool solved = false;
  try {
    if (dim <= kDenseRing) {
      throw std::runtime_error("small ring: dense solve");
    }
    const WoodburyRing ring(spring, diag, closed, extras, kappas, false);
    if (ring.ok()) {
      stepRing = ring.solve(rhs);
      // Cutting the ring open can leave a Schur pivot near zero next to
      // the barrier, and the Woodbury solve then loses digits; iterative
      // refinement against the exact operator recovers them.
      for (int pass = 0; pass < 4; ++pass) {
        std::vector<VectorXd> r = applyShifted(stepRing);
        for (size_t j = 0; j < r.size(); ++j) {
          r[j] = rhs[j] - r[j];
        }
        const double rn = std::sqrt(dot(r, r));
        if (!std::isfinite(rn)) {
          break;
        }
        if (rn <= 1e-10 * std::max(1.0, rhsNorm)) {
          solved = true;
          break;
        }
        const std::vector<VectorXd> dx = ring.solve(r);
        for (size_t j = 0; j < stepRing.size(); ++j) {
          stepRing[j] += dx[j];
        }
      }
    }
  } catch (const std::runtime_error &) {
    solved = false;
  }
  if (!solved && dim <= kDenseRing) {
    // A small ring is cheap to solve densely when the chain cannot.
    const long n = static_cast<long>(grad.size());
    ColMajorXd jt = ColMajorXd::Zero(dim, dim);
    for (long col = 0; col < dim; ++col) {
      std::vector<VectorXd> e(static_cast<size_t>(n), VectorXd::Zero(f));
      e[static_cast<size_t>(col / f)](col % f) = 1.0;
      jt.col(col) = packBeads(applyShifted(e));
    }
    const Eigen::PartialPivLU<ColMajorXd> lu(jt);
    const VectorXd x = lu.solve(packBeads(rhs));
    stepRing.assign(static_cast<size_t>(n), VectorXd::Zero(f));
    for (long j = 0; j < n; ++j) {
      stepRing[static_cast<size_t>(j)] = x.segment(j * f, f);
    }
    solved = x.array().isFinite().all();
  }
  if (!solved && stepRing.empty()) {
    return VectorXd();
  }
  // The cycle component stays: a gradient along it has to be stepped out.
  VectorXd step = packBeads(stepRing);
  for (const auto &r : nullRing) {
    const VectorXd rf = packBeads(r);
    if (rf.size() == step.size()) {
      step -= step.dot(rf) * rf;
    }
  }
  if (!step.array().isFinite().all()) {
    return VectorXd();
  }
  return step;
}

// Non-increasing isotonic regression (pool adjacent violators).
std::vector<double> pavaNonIncreasing(const std::vector<double> &values) {
  const long n = static_cast<long>(values.size());
  struct Block {
    long start;
    long end;
    double average;
    double weight;
  };
  std::vector<Block> blocks;
  blocks.reserve(static_cast<size_t>(n));
  for (long i = 0; i < n; ++i) {
    blocks.push_back(Block{i, i, -values[static_cast<size_t>(i)], 1.0});
    while (blocks.size() >= 2 &&
           blocks[blocks.size() - 2].average > blocks.back().average) {
      const Block b = blocks.back();
      blocks.pop_back();
      Block &a = blocks.back();
      const double w = a.weight + b.weight;
      a.average = (a.average * a.weight + b.average * b.weight) / w;
      a.end = b.end;
      a.weight = w;
    }
  }
  std::vector<double> out(static_cast<size_t>(n));
  for (const Block &b : blocks) {
    for (long i = b.start; i <= b.end; ++i) {
      out[static_cast<size_t>(i)] = -b.average;
    }
  }
  return out;
}

// Beads 0..last move along the unit vector dir until that coordinate is
// monotone.
// The sense follows the two ends, so an arbitrary eigenvector sign is harmless.
void projectMonotonePrefix(std::vector<VectorXd> &q, long last,
                           const VectorXd &dir) {
  if (last < 1 || q.empty() || dir.size() != q.front().size()) {
    return;
  }
  const long m = last;
  std::vector<double> s(static_cast<size_t>(m + 1));
  for (long j = 0; j <= m; ++j) {
    s[static_cast<size_t>(j)] = dir.dot(q[static_cast<size_t>(j)]);
  }
  std::vector<double> target;
  if (s.front() >= s.back()) {
    target = pavaNonIncreasing(s);
  } else {
    std::vector<double> flipped(s.size());
    for (size_t j = 0; j < s.size(); ++j) {
      flipped[j] = -s[j];
    }
    target = pavaNonIncreasing(flipped);
    for (double &v : target) {
      v = -v;
    }
  }
  for (long j = 0; j <= m; ++j) {
    q[static_cast<size_t>(j)] +=
        (target[static_cast<size_t>(j)] - s[static_cast<size_t>(j)]) * dir;
  }
}

/// The rigid motions of a whole ring: the three translations and each free
/// rotation about the ring's centre of mass, built from the beads
/// themselves and orthonormalised (rank-revealing, so a linear ring keeps
/// two rotations). A rotation moves bead j along sqrt(m) e x (r_j - centre),
/// which differs from bead to bead, so a single structure's generator copied
/// to every bead is not a null vector of the ring Hessian.
std::vector<std::vector<VectorXd>>
ringRigidBasis(const std::vector<VectorXd> &q,
               const std::vector<double> &sqrtMasses, const VectorXd &reference,
               const std::array<bool, 3> &rotations) {
  std::vector<std::vector<VectorXd>> out;
  const long nAtoms = static_cast<long>(sqrtMasses.size());
  if (nAtoms == 0 || q.empty() || reference.size() != 3 * nAtoms ||
      q.front().size() != 3 * nAtoms) {
    return out;
  }
  const long nb = static_cast<long>(q.size());
  // Cartesian positions and the ring's centre of mass.
  Eigen::Vector3d centre = Eigen::Vector3d::Zero();
  double total = 0.0;
  for (long j = 0; j < nb; ++j) {
    for (long k = 0; k < nAtoms; ++k) {
      const double sm = sqrtMasses[static_cast<size_t>(k)];
      const Eigen::Vector3d r =
          reference.segment<3>(3 * k) +
          q[static_cast<size_t>(j)].segment<3>(3 * k) / sm;
      centre += sm * sm * r;
      total += sm * sm;
    }
  }
  centre /= total;
  std::vector<int> kinds{0, 1, 2};
  for (int c = 0; c < 3; ++c) {
    if (rotations[static_cast<size_t>(c)]) {
      kinds.push_back(3 + c);
    }
  }
  const long dim = nb * 3 * nAtoms;
  MatrixXd g = MatrixXd::Zero(dim, static_cast<long>(kinds.size()));
  for (size_t col = 0; col < kinds.size(); ++col) {
    const int kind = kinds[col];
    for (long j = 0; j < nb; ++j) {
      for (long k = 0; k < nAtoms; ++k) {
        const double sm = sqrtMasses[static_cast<size_t>(k)];
        Eigen::Vector3d d = Eigen::Vector3d::Zero();
        if (kind < 3) {
          d(kind) = sm;
        } else {
          const Eigen::Vector3d r =
              reference.segment<3>(3 * k) +
              q[static_cast<size_t>(j)].segment<3>(3 * k) / sm;
          Eigen::Vector3d e = Eigen::Vector3d::Zero();
          e(kind - 3) = 1.0;
          d = sm * e.cross(r - centre);
        }
        g.block(j * 3 * nAtoms + 3 * k, static_cast<long>(col), 3, 1) = d;
      }
    }
  }
  const ColMajorXd gc = g;
  const Eigen::ColPivHouseholderQR<ColMajorXd> qr(gc);
  const long rank = qr.rank();
  const ColMajorXd basis = qr.householderQ() * ColMajorXd::Identity(dim, rank);
  for (long r = 0; r < rank; ++r) {
    std::vector<VectorXd> u(static_cast<size_t>(nb));
    for (long j = 0; j < nb; ++j) {
      u[static_cast<size_t>(j)] =
          basis.col(r).segment(j * 3 * nAtoms, 3 * nAtoms);
    }
    out.push_back(std::move(u));
  }
  return out;
}

struct NewtonOut {
  std::vector<VectorXd> beads;
  std::vector<double> energies;
  double ringPotential = 0.0;
  double bN = 0.0;
  long iterations = 0;
  bool converged = false;
  // Half ring, stationary, and not index 1. The odd-mode probe decides
  // whether the search continues on the whole ring.
  bool stalledHalf = false;
};

// Index-1 Newton on one ring. `x` holds N beads. A half ring optimises beads
// 0..N/2 and mirrors them. Trust is the largest bead displacement.
NewtonOut newtonInstanton(std::vector<VectorXd> guess, double c,
                          const MatrixXd &hessSaddle,
                          const RateInstantonOptions &options,
                          const BatchPotential &potential) {
  const long nBeads = options.beads;
  // Fold only an even ring that already matches under j -> N - j. An empty
  // guess is seeded into that shape before this call.
  bool half = options.halfRing && nBeads % 2 == 0 &&
              static_cast<long>(guess.size()) == nBeads;
  if (half) {
    const long m = nBeads / 2;
    for (long j = 1; j < m; ++j) {
      if ((guess[static_cast<size_t>(j)] -
           guess[static_cast<size_t>(nBeads - j)])
              .norm() > 1e-8) {
        half = false;
        break;
      }
    }
  }
  std::vector<VectorXd> x;
  if (half) {
    const long m = nBeads / 2;
    x.resize(static_cast<size_t>(m + 1));
    for (long j = 0; j <= m; ++j) {
      x[static_cast<size_t>(j)] = guess[static_cast<size_t>(j)];
    }
  } else {
    x = std::move(guess);
  }
  const long f = x.front().size();
  const MatrixXd hS = (0.5 * (hessSaddle + hessSaddle.transpose())).eval();
  // The rigid quotient: translations and free rotations of the whole ring
  // about its centre of mass, from the current beads, orthonormalised.
  const long nAtoms = static_cast<long>(options.rigidSqrtMasses.size());
  const bool quotient = nAtoms > 0 && 3 * nAtoms == hS.rows() &&
                        options.rigidReference.size() == 3 * nAtoms;
  auto ringRigid = [&](const std::vector<VectorXd> &q) {
    if (!quotient) {
      return std::vector<std::vector<VectorXd>>{};
    }
    return ringRigidBasis(q, options.rigidSqrtMasses, options.rigidReference,
                          options.rigidRotations);
  };
  std::vector<std::vector<VectorXd>> nullRing;
  // A one-dimensional well is already the saddle curvature. In more
  // dimensions the turning points are not the saddle, so each bead starts
  // from its own curvature and the Bofill update carries it.
  std::vector<MatrixXd> physical;
  if (f == 1 || options.initialHessians != "finite_difference") {
    physical.assign(x.size(), hS);
  } else {
    const double eps = options.lanczosStep > 0.0 ? options.lanczosStep : 1e-4;
    for (size_t j = 0; j < x.size(); ++j) {
      physical[j] = fdPhysicalHessian(x[j], potential, eps);
    }
  }
  // Exact zeros of the saddle Hessian are rigid displacements. The same
  // vector on every bead is a null vector of the chain, and the block LU
  // then returns an arbitrary step. Hold those directions at a
  // spring-sized curvature for the solve only.
  std::vector<VectorXd> rigid;
  {
    const Eigen::SelfAdjointEigenSolver<ColMajorXd> esRigid(hS);
    const auto ev = esRigid.eigenvalues();
    const double span = std::max(std::abs(ev(0)), std::abs(ev(ev.size() - 1)));
    const double lim = 1e-6 * std::max(1.0, span);
    for (long i = 0; i < ev.size(); ++i) {
      if (std::abs(ev(i)) <= lim) {
        rigid.push_back(VectorXd(esRigid.eigenvectors().col(i)));
      }
    }
  }

  struct Obj {
    std::vector<VectorXd> grad;
    std::vector<VectorXd> gradPot;
    std::vector<double> energies;
    double u = 0.0;
  };
  auto objective = [&](const std::vector<VectorXd> &q) {
    Obj out;
    if (half) {
      // The potential is evaluated on beads 0..N/2. The closed ring is that
      // chain plus its mirror, and the folded gradient is half the derivative
      // of the closed-ring energy.
      const long m = static_cast<long>(q.size()) - 1;
      const long n = 2 * m;
      std::vector<double> vu;
      std::vector<VectorXd> gu;
      potential(q, vu, gu);
      for (double &vj : vu) {
        vj -= options.energyShift;
      }
      if (vu.size() != q.size() || gu.size() != q.size()) {
        throw std::runtime_error(
            "rate instanton: potential returned the wrong count");
      }
      std::vector<VectorXd> full(static_cast<size_t>(n));
      std::vector<double> vFull(static_cast<size_t>(n));
      std::vector<VectorXd> gFull(static_cast<size_t>(n));
      for (long j = 0; j <= m; ++j) {
        full[static_cast<size_t>(j)] = q[static_cast<size_t>(j)];
        vFull[static_cast<size_t>(j)] = vu[static_cast<size_t>(j)];
        gFull[static_cast<size_t>(j)] = gu[static_cast<size_t>(j)];
      }
      for (long j = 1; j < m; ++j) {
        full[static_cast<size_t>(n - j)] = q[static_cast<size_t>(j)];
        vFull[static_cast<size_t>(n - j)] = vu[static_cast<size_t>(j)];
        gFull[static_cast<size_t>(n - j)] = gu[static_cast<size_t>(j)];
      }
      std::vector<VectorXd> gRing(static_cast<size_t>(n));
      double uFull = 0.0;
      double springE = 0.0;
      for (long j = 0; j < n; ++j) {
        const long prev = (j + n - 1) % n;
        const long next = (j + 1) % n;
        gRing[static_cast<size_t>(j)] =
            gFull[static_cast<size_t>(j)] +
            c * (2.0 * full[static_cast<size_t>(j)] -
                 full[static_cast<size_t>(prev)] -
                 full[static_cast<size_t>(next)]);
        springE +=
            (full[static_cast<size_t>(next)] - full[static_cast<size_t>(j)])
                .squaredNorm();
        uFull += vFull[static_cast<size_t>(j)];
      }
      uFull += 0.5 * c * springE;
      out.u = 0.5 * uFull;
      out.grad.resize(q.size());
      out.grad.front() = 0.5 * gRing.front();
      out.grad.back() = 0.5 * gRing[static_cast<size_t>(m)];
      for (long j = 1; j < m; ++j) {
        out.grad[static_cast<size_t>(j)] =
            0.5 *
            (gRing[static_cast<size_t>(j)] + gRing[static_cast<size_t>(n - j)]);
      }
      out.gradPot = std::move(gu);
      out.energies = std::move(vu);
    } else {
      const RingEval ev = evaluateRing(q, c, potential, options.energyShift);
      out.u = ev.u;
      out.grad = ev.grad;
      out.gradPot = physicalGradient(q, ev.grad, c);
      out.energies = ev.v;
    }
    return out;
  };

  struct View {
    bool ok = false;
    std::vector<RingMode> ritz; // lowest Ritz pairs, ascending
    std::vector<MatrixXd> diag; // ring blocks the step solves with
    VectorXd tau;
    Climb climb;
    double gmax = 0.0;
  };
  // Turning-point blocks store half the closed-ring derivative. The
  // residual that stops the climb is the closed-ring residual.
  auto closedGmax = [&](const Obj &ev) {
    if (!half || ev.grad.size() < 2) {
      return largestBeadNorm(ev.grad);
    }
    double big = 2.0 * ev.grad.front().norm();
    big = std::max(big, 2.0 * ev.grad.back().norm());
    const long last = static_cast<long>(ev.grad.size()) - 1;
    for (long j = 1; j < last; ++j) {
      big = std::max(big, ev.grad[static_cast<size_t>(j)].norm());
    }
    return big;
  };
  // Lowest modes of the ring Hessian from matrix-vector products alone,
  // started along the saddle's unstable direction on every bead.
  std::vector<VectorXd> ritzStart(x.size(), VectorXd::Zero(f));
  VectorXd climbDir(f);
  {
    const ColMajorXd hs0 = hS;
    const Eigen::SelfAdjointEigenSolver<ColMajorXd> es0(hs0);
    climbDir = es0.eigenvectors().col(0);
    for (auto &v : ritzStart) {
      v = climbDir;
    }
  }
  // The climb follows the mode that overlaps the last one, as in dimer and
  // minimum-mode following: the lowest curvature can switch to another
  // channel, and climbing it walks the ring to a neighbouring saddle.
  std::vector<VectorXd> track = ritzStart;
  {
    const double n0 = std::sqrt(dot(track, track));
    if (n0 > 0.0) {
      scale(track, 1.0 / n0);
    }
  }
  // A half chain that folds back is a second bounce. The single instanton
  // is monotone in the saddle's unstable direction.
  if (half) {
    projectMonotonePrefix(x, static_cast<long>(x.size()) - 1, climbDir);
  }
  auto viewOf = [&](const Obj &ev) {
    nullRing = ringRigid(x);
    View v;
    v.tau = half ? VectorXd() : timeTranslation(x);
    std::vector<MatrixXd> pinned = physical;
    if (!rigid.empty()) {
      const double pin = std::max(1.0, c);
      for (auto &block : pinned) {
        for (const auto &r : rigid) {
          block.noalias() += pin * (r * r.transpose());
        }
      }
    }
    v.diag = ringDiagonal(pinned, c, half);
    auto apply = [&](const std::vector<VectorXd> &vec) {
      return applyDiagonal(v.diag, c, !half, vec);
    };
    const long dim = static_cast<long>(x.size()) * f;
    if (dim <= kDenseRing) {
      // A small ring takes its whole spectrum densely: every negative
      // curvature in every symmetry sector, exactly.
      const long nb = static_cast<long>(x.size());
      ColMajorXd big = ColMajorXd::Zero(dim, dim);
      for (long col = 0; col < dim; ++col) {
        std::vector<VectorXd> e(static_cast<size_t>(nb), VectorXd::Zero(f));
        e[static_cast<size_t>(col / f)](col % f) = 1.0;
        big.col(col) = packBeads(apply(e));
      }
      const ColMajorXd sym = 0.5 * (big + big.transpose());
      const Eigen::SelfAdjointEigenSolver<ColMajorXd> es(sym);
      v.ritz.clear();
      for (long i = 0; i < dim; ++i) {
        RingMode m;
        m.theta = es.eigenvalues()(i);
        m.residual = 0.0;
        m.vector.assign(static_cast<size_t>(nb), VectorXd::Zero(f));
        for (long j = 0; j < nb; ++j) {
          m.vector[static_cast<size_t>(j)] =
              es.eigenvectors().col(i).segment(j * f, f);
        }
        v.ritz.push_back(std::move(m));
      }
    } else {
      // The number of negative curvatures off the cycle, exactly, from the
      // inertia of the same chain with the cycle lifted by c.
      long negatives = 0;
      {
        std::vector<std::vector<VectorXd>> lift;
        std::vector<double> kap;
        if (v.tau.size() == dim) {
          std::vector<VectorXd> tauRing(x.size(), VectorXd::Zero(f));
          for (size_t j = 0; j < x.size(); ++j) {
            tauRing[j] = v.tau.segment(static_cast<long>(j) * f, f);
          }
          lift.push_back(std::move(tauRing));
          kap.push_back(c);
        }
        try {
          negatives =
              WoodburyRing(c, v.diag, !half, lift, kap, true).negative();
        } catch (const std::runtime_error &) {
          negatives = 0;
        }
      }
      // Lanczos deepens until every one of those negative curvatures is a
      // resolved Ritz pair; a flip on an unresolved vector corrupts the step.
      // The ring spectrum reaches 4 c, so residuals scale with c.
      const double resolvedTol = 1e-8 * std::max(1.0, 4.0 * c);
      auto resolvedNegatives = [&](const std::vector<RingMode> &modes) {
        long count = 0;
        for (const auto &m : modes) {
          if (m.theta < 0.0 && m.residual <= resolvedTol) {
            ++count;
          }
        }
        return count;
      };
      long steps = std::min(dim, std::max(60L, 4 * (negatives + 2)));
      for (;;) {
        // J commutes with the ring's mirror about any bead, so a start that
        // is mirror-symmetric spans no antisymmetric mode at all; a fixed
        // pseudo-random admixture reaches every symmetry sector.
        std::vector<VectorXd> start = ritzStart;
        std::uint64_t h = 0x9E3779B97F4A7C15ULL;
        for (auto &bead : start) {
          for (long a2 = 0; a2 < bead.size(); ++a2) {
            h ^= h << 13;
            h ^= h >> 7;
            h ^= h << 17;
            bead(a2) += 1e-2 * (static_cast<double>(h >> 11) * 0x1.0p-53 - 0.5);
          }
        }
        v.ritz = lowestRingModes(apply, std::move(start), steps);
        if (steps >= std::min(dim, kRitzCap) ||
            resolvedNegatives(v.ritz) >= negatives) {
          break;
        }
        steps = std::min(std::min(dim, kRitzCap), 2 * steps);
      }
      // Only resolved pairs enter the classification and the flips.
      v.ritz.erase(std::remove_if(v.ritz.begin(), v.ritz.end(),
                                  [](const RingMode &m) {
                                    return m.residual >
                                           1e-6 *
                                               std::max(1.0, std::abs(m.theta));
                                  }),
                   v.ritz.end());
    }
    v.ok = !v.ritz.empty();
    if (v.ok) {
      // Classify from the Ritz values: the lowest mode off the cycle is
      // the climb; further negatives below a thousandth of it count.
      const double cut = -1e-8 * std::max(1.0, c);
      v.climb = Climb{};
      std::vector<VectorXd> tauRing;
      if (v.tau.size() == dim) {
        tauRing.assign(x.size(), VectorXd::Zero(f));
        for (size_t j = 0; j < x.size(); ++j) {
          tauRing[j] = v.tau.segment(static_cast<long>(j) * f, f);
        }
      }
      auto onCycle = [&](const std::vector<VectorXd> &m) {
        if (!tauRing.empty() && std::abs(dot(m, tauRing)) > 0.5) {
          return true;
        }
        for (const auto &r : nullRing) {
          if (std::abs(dot(m, r)) > 0.5) {
            return true;
          }
        }
        return false;
      };
      double best = 0.0;
      long lowest = -1;
      for (size_t i = 0; i < v.ritz.size(); ++i) {
        if (onCycle(v.ritz[i].vector)) {
          continue;
        }
        if (lowest < 0 ||
            v.ritz[i].theta < v.ritz[static_cast<size_t>(lowest)].theta) {
          lowest = static_cast<long>(i);
        }
        if (v.ritz[i].theta < 0.0 && track.size() == v.ritz[i].vector.size()) {
          const double o = std::abs(dot(v.ritz[i].vector, track));
          if (o > best) {
            best = o;
            v.climb.index = static_cast<long>(i);
          }
        }
      }
      if (v.climb.index < 0 || best < kTrackOverlap) {
        v.climb.index = lowest;
      }
      if (v.climb.index >= 0) {
        v.climb.curvature = v.ritz[static_cast<size_t>(v.climb.index)].theta;
      }
      if (v.climb.curvature < 0.0) {
        for (size_t i = 0; i < v.ritz.size(); ++i) {
          if (onCycle(v.ritz[i].vector)) {
            continue;
          }
          if (v.ritz[i].theta < cut &&
              v.ritz[i].theta <= 1e-3 * v.climb.curvature) {
            ++v.climb.negative;
          }
        }
      }
      v.gmax = closedGmax(ev);
      if (v.climb.index >= 0) {
        ritzStart = v.ritz[static_cast<size_t>(v.climb.index)].vector;
        if (v.climb.curvature < 0.0) {
          track = ritzStart;
          if (track.size() == ritzStart.size() && dot(track, track) > 0.0) {
            scale(track, 1.0 / std::sqrt(dot(track, track)));
          }
        }
      }
    }
    return v;
  };
  auto done = [&](const View &v) {
    return v.ok && v.gmax < options.forceTolerance && v.climb.negative == 1 &&
           v.climb.curvature < 0.0;
  };

  Obj cur = objective(x);
  double trust = options.maxStep;
  long entries = 0;
  bool converged = false;
  bool stalledHalf = false;
  const double trustFloor = std::min(1e-4, options.maxStep);
  // Finite-difference rebuilds of the bead blocks when the search stalls.
  constexpr int kMaxHessianRefreshes = 3;
  int refreshes = 0;
  bool exactAtX = false;
  for (long it = 0; it < options.maxIterations; ++it) {
    ++entries;
    if (f == 1) {
      const double eps = options.lanczosStep > 0.0 ? options.lanczosStep : 1e-4;
      for (size_t j = 0; j < x.size(); ++j) {
        physical[j] = fdPhysicalHessian(x[j], potential, eps);
      }
    }
    const View v = viewOf(cur);
    if (done(v)) {
      converged = true;
      break;
    }
    // A shorter Newton step when the quadratic model does not match. A step
    // that only reduces the residual is a walk into a well.
    auto accept = [&](VectorXd dir) {
      if (dir.size() == 0 || !dir.array().isFinite().all()) {
        return false;
      }
      const double big = packedBeadNorm(dir, f);
      if (!(big > 0.0)) {
        return false;
      }
      std::vector<VectorXd> trial = x;
      addPacked(trial, dir);
      if (half) {
        projectMonotonePrefix(trial, static_cast<long>(trial.size()) - 1,
                              climbDir);
      }
      if (!finiteBeads(trial)) {
        return false;
      }
      VectorXd taken(dir.size());
      for (size_t k = 0; k < x.size(); ++k) {
        taken.segment(static_cast<long>(k) * f, f) =
            trial[k] - x[static_cast<size_t>(k)];
      }
      if (!(packedBeadNorm(taken, f) > 0.0) ||
          !taken.array().isFinite().all()) {
        return false;
      }
      dir = taken;
      Obj next = objective(trial);
      if (!finiteBeads(next.grad) || !std::isfinite(next.u)) {
        return false;
      }
      bool ratioOk = false;
      double ratio = 0.0;
      if (v.ok && static_cast<long>(x.size()) * f == dir.size()) {
        const VectorXd gflat = packBeads(cur.grad);
        std::vector<VectorXd> dirRing(x.size(), VectorXd::Zero(f));
        for (size_t k = 0; k < x.size(); ++k) {
          dirRing[k] = dir.segment(static_cast<long>(k) * f, f);
        }
        const VectorXd jd = packBeads(applyDiagonal(v.diag, c, !half, dirRing));
        const double pred = gflat.dot(dir) + 0.5 * dir.dot(jd);
        const double actual = next.u - cur.u;
        // Predicted and actual changes below the energy resolution are
        // roundoff. Their ratio is not a disagreement with the model,
        // and rejecting it repeats one micro-step until the budget ends.
        const double unresolved = 1e-12 * std::max(1.0, std::abs(cur.u));
        if (std::abs(pred) <= unresolved) {
          const bool quiet = std::abs(actual) <= unresolved &&
                             packedBeadNorm(dir, f) <= trustFloor;
          ratio = quiet ? 1.0 : 0.0;
        } else {
          ratio = actual / pred;
        }
        ratioOk = std::isfinite(ratio) && ratio >= 0.1 && ratio <= 3.0;
        // Near the saddle the predicted change sits at the round-off of
        // U_N and the ratio is noise; there the step stands on the residual.
        double magnitude = std::abs(cur.u);
        for (const double vj : cur.energies) {
          magnitude += std::abs(vj);
        }
        const double roundoff =
            1e3 * std::numeric_limits<double>::epsilon() * magnitude;
        if (!ratioOk && std::abs(pred) < roundoff &&
            closedGmax(next) <= closedGmax(cur)) {
          ratioOk = true;
          ratio = 1.0;
        }
      }
      if (!ratioOk) {
        return false;
      }
      for (size_t k = 0; k < x.size(); ++k) {
        bofillUpdate(physical[k], trial[k] - x[k],
                     next.gradPot[k] - cur.gradPot[k]);
      }
      x = std::move(trial);
      cur = std::move(next);
      exactAtX = false;
      if (ratioOk && ratio > 0.75 && ratio < 1.25 &&
          packedBeadNorm(dir, f) >= 0.99 * trust) {
        trust = std::min(2.0 * trust, options.maxStep);
      }
      return true;
    };
    // A converged gradient is classified with exact bead Hessians, as
    // i-PI's hessian_final does: the Bofill blocks can carry negative
    // curvatures the surface does not have. A second negative curvature
    // that survives the rebuild is a higher-index stationary ring, where
    // the flipped Newton step vanishes; a trust-sized displacement down
    // that mode leaves it.
    if (v.ok && v.gmax < options.forceTolerance && v.climb.negative != 1 &&
        !exactAtX) {
      const double eps = options.lanczosStep > 0.0 ? options.lanczosStep : 1e-4;
      for (size_t j = 0; j < x.size(); ++j) {
        physical[j] = fdPhysicalHessian(x[j], potential, eps);
      }
      exactAtX = true;
      continue;
    }
    if (v.ok && v.gmax < options.forceTolerance && v.climb.negative > 1) {
      long down = -1;
      for (size_t i = 0; i < v.ritz.size(); ++i) {
        const auto &m = v.ritz[i];
        if (static_cast<long>(i) == v.climb.index || !(m.theta < 0.0)) {
          continue;
        }
        bool held = false;
        if (v.tau.size() == static_cast<long>(x.size()) * f) {
          held = std::abs(packBeads(m.vector).dot(v.tau)) > 0.5;
        }
        for (const auto &r : nullRing) {
          held = held || std::abs(dot(m.vector, r)) > 0.5;
        }
        if (!held &&
            (down < 0 || m.theta < v.ritz[static_cast<size_t>(down)].theta)) {
          down = static_cast<long>(i);
        }
      }
      if (down >= 0) {
        trust = options.maxStep;
        const VectorXd mode =
            packBeads(v.ritz[static_cast<size_t>(down)].vector);
        if (accept(mode) || accept(-mode)) {
          continue;
        }
      }
    }
    // A half ring that is stationary at the wrong index cannot see a mode
    // odd under the mirror. Leave it for the odd-mode probe instead of
    // spending the remaining iterations on this point.
    if (half && options.checkOddSector && exactAtX && v.ok &&
        v.gmax < options.forceTolerance && v.climb.negative != 1) {
      stalledHalf = true;
      break;
    }
    const VectorXd step = chainIndexOneStep(v.ritz, v.climb, v.diag, c, !half,
                                            cur.grad, v.tau, nullRing);
    bool moved = false;
    VectorXd dir = step;
    // Clip to the trust radius once, then halve that capped step.
    if (dir.size() > 0) {
      const double big0 = packedBeadNorm(dir, f);
      if (big0 > trust) {
        dir *= trust / big0;
      }
    }
    for (int bt = 0; bt < 4 && !moved; ++bt) {
      moved = accept(dir);
      dir *= 0.5;
    }
    if (!moved) {
      trust = std::max(0.5 * trust, trustFloor);
      // At the trust floor the Bofill blocks no longer model the ring and
      // every step is refused; rebuild them from finite differences, 2 f
      // gradient calls per bead, a few times at most, and start the trust
      // region again.
      if (trust <= trustFloor && refreshes < kMaxHessianRefreshes) {
        ++refreshes;
        const double eps =
            options.lanczosStep > 0.0 ? options.lanczosStep : 1e-4;
        for (size_t j = 0; j < x.size(); ++j) {
          physical[j] = fdPhysicalHessian(x[j], potential, eps);
        }
        trust = options.maxStep;
      }
    }
  }
  if (!converged && done(viewOf(cur))) {
    converged = true;
  }

  NewtonOut out;
  out.iterations = entries;
  out.converged = converged;
  out.stalledHalf = stalledHalf;
  if (half) {
    const long m = static_cast<long>(x.size()) - 1;
    const long n = 2 * m;
    out.beads.resize(static_cast<size_t>(n));
    out.energies.assign(static_cast<size_t>(n), 0.0);
    for (long j = 0; j <= m; ++j) {
      out.beads[static_cast<size_t>(j)] = x[static_cast<size_t>(j)];
      out.energies[static_cast<size_t>(j)] =
          cur.energies[static_cast<size_t>(j)];
    }
    for (long j = 1; j < m; ++j) {
      out.beads[static_cast<size_t>(n - j)] = x[static_cast<size_t>(j)];
      out.energies[static_cast<size_t>(n - j)] =
          cur.energies[static_cast<size_t>(j)];
    }
    out.ringPotential = 2.0 * cur.u;
  } else {
    out.beads = std::move(x);
    out.ringPotential = cur.u;
    out.energies = std::move(cur.energies);
  }
  out.bN = 0.0;
  for (size_t j = 0; j < out.beads.size(); ++j) {
    const size_t next = (j + 1) % out.beads.size();
    out.bN += (out.beads[next] - out.beads[j]).squaredNorm();
  }
  return out;
}

// A one-dimensional bead carries its own curvature, so the requested
// temperature is solved first. Below 0.75 Tc an empty guess that does not
// converge walks down from 0.85 Tc, and that walk replaces the first try.
RateInstanton optimizeRateByNewton(const VectorXd &saddle,
                                   const MatrixXd &hessSaddle, double beta,
                                   std::vector<VectorXd> guess,
                                   const BatchPotential &potential,
                                   const RateInstantonOptions &options) {
  const long nBeads = options.beads;
  if (nBeads < 4 || !(beta > 0.0) || hessSaddle.rows() != saddle.size()) {
    throw std::invalid_argument(
        "optimizeRateInstanton: need N >= 4, beta > 0 and a saddle Hessian "
        "of the saddle's dimension");
  }
  RateInstanton inst;
  inst.beta = beta;
  inst.betaN = beta / static_cast<double>(nBeads);
  inst.temperature = 1.0 / (kBoltzmann * beta);
  inst.crossover = crossoverTemperature(hessSaddle);
  if (!(inst.temperature < inst.crossover)) {
    throw std::invalid_argument(
        "optimizeRateInstanton: T is at or above the crossover temperature; "
        "the ring collapses onto the saddle and steepest descent needs the "
        "parabolic barrier correction, of which classical transition-state "
        "theory is only the one-bead limit");
  }

  const bool cool = static_cast<long>(guess.size()) != nBeads &&
                    inst.temperature < 0.75 * inst.crossover;
  // One dimension solves the requested temperature before any walk. A
  // higher-dimensional empty guess below 0.75 Tc keeps the walk, which
  // starts from a curvature copied off the saddle.
  if (saddle.size() == 1 || !cool) {
    if (static_cast<long>(guess.size()) != nBeads) {
      const ColMajorXd hS = 0.5 * (hessSaddle + hessSaddle.transpose());
      const Eigen::SelfAdjointEigenSolver<ColMajorXd> es(hS);
      guess = cosineSeed(saddle, es.eigenvectors().col(0), es.eigenvalues()(0),
                         inst.temperature, inst.crossover, nBeads, potential);
    }
    const double bnh = inst.betaN * kHbar;
    const double spring = 1.0 / (bnh * bnh);
    const NewtonOut got = newtonInstanton(std::move(guess), spring, hessSaddle,
                                          options, potential);
    inst.beads = got.beads;
    inst.energies = got.energies;
    inst.ringPotential = got.ringPotential;
    inst.bN = got.bN;
    inst.iterations = got.iterations;
    inst.converged = got.converged;
    if (inst.converged || !cool) {
      return inst;
    }
  }
  if (cool) {
    std::vector<double> temps;
    for (double t = 0.85 * inst.crossover; t > inst.temperature * 1.05;
         t *= 0.75) {
      temps.push_back(t);
    }
    temps.push_back(inst.temperature);
    std::vector<VectorXd> beads;
    RateInstanton last;
    long used = 0;
    bool targetRan = false;
    const ColMajorXd hS = 0.5 * (hessSaddle + hessSaddle.transpose());
    const Eigen::SelfAdjointEigenSolver<ColMajorXd> es(hS);
    for (size_t s = 0; s < temps.size(); ++s) {
      const long remain = options.maxIterations - used;
      if (remain <= 0) {
        break;
      }
      RateInstantonOptions opt = options;
      opt.maxIterations = remain;
      const bool target = s + 1 == temps.size();
      opt.checkOddSector = options.checkOddSector && target;
      const double betaStage = target ? beta : 1.0 / (kBoltzmann * temps[s]);
      std::vector<VectorXd> stageGuess = beads;
      if (static_cast<long>(stageGuess.size()) != nBeads) {
        stageGuess =
            cosineSeed(saddle, es.eigenvectors().col(0), es.eigenvalues()(0),
                       temps[s], inst.crossover, nBeads, potential);
      }
      last = optimizeRateInstanton(saddle, hessSaddle, betaStage,
                                   std::move(stageGuess), potential, opt);
      used += last.iterations;
      beads = last.beads;
      if (target) {
        targetRan = true;
      }
    }
    last.iterations = used;
    if (targetRan) {
      last.beta = beta;
      last.betaN = beta / static_cast<double>(nBeads);
      last.temperature = inst.temperature;
      last.crossover = inst.crossover;
    } else {
      last.converged = false;
    }
    return last;
  }

  if (static_cast<long>(guess.size()) != nBeads) {
    const ColMajorXd hS = 0.5 * (hessSaddle + hessSaddle.transpose());
    const Eigen::SelfAdjointEigenSolver<ColMajorXd> es(hS);
    guess = cosineSeed(saddle, es.eigenvectors().col(0), es.eigenvalues()(0),
                       inst.temperature, inst.crossover, nBeads, potential);
  }
  const double bnh = inst.betaN * kHbar;
  const double spring = 1.0 / (bnh * bnh);
  NewtonOut got =
      newtonInstanton(std::move(guess), spring, hessSaddle, options, potential);
  // A half ring cannot see modes odd under j -> N-j, so two copies of the
  // instanton are a stationary point it can accept or sit on. An unstable
  // odd mode there sends the search onto the whole ring from a kick along
  // that mode.
  const long nGot = static_cast<long>(got.beads.size());
  bool mirrored = (got.converged || got.stalledHalf) && options.halfRing &&
                  options.checkOddSector && nGot == nBeads && nGot % 2 == 0;
  for (long j = 1; mirrored && j < nGot / 2; ++j) {
    mirrored = (got.beads[static_cast<size_t>(j)] -
                got.beads[static_cast<size_t>(nGot - j)])
                   .norm() <= 1e-8;
  }
  if (mirrored) {
    const Eigen::SelfAdjointEigenSolver<MatrixXd> es0(
        0.5 * (hessSaddle + hessSaddle.transpose()));
    const double barrierCurvature = std::abs(es0.eigenvalues()(0));
    const RingEval here =
        evaluateRing(got.beads, spring, potential, options.energyShift);
    std::vector<VectorXd> oddMode;
    const double oddCurv =
        lowestOddMode(got.beads, here, spring, potential, oddMode,
                      options.lanczosFirst, options.lanczosStep);
    if (oddCurv < -1e-3 * barrierCurvature &&
        oddMode.size() == got.beads.size()) {
      std::vector<VectorXd> kicked = got.beads;
      const double kick = std::sqrt(2.0 / (inst.betaN * -oddCurv));
      for (size_t k = 0; k < kicked.size(); ++k) {
        kicked[k] += kick * oddMode[k];
      }
      RateInstantonOptions whole = options;
      whole.halfRing = false;
      const long before = got.iterations;
      got = newtonInstanton(std::move(kicked), spring, hessSaddle, whole,
                            potential);
      got.iterations += before;
    }
  }
  inst.beads = got.beads;
  inst.energies = got.energies;
  inst.ringPotential = got.ringPotential;
  inst.bN = got.bN;
  inst.iterations = got.iterations;
  inst.converged = got.converged;
  return inst;
}

} // namespace

RateInstanton optimizeRateInstanton(const VectorXd &saddle,
                                    const MatrixXd &hessSaddle, double beta,
                                    std::vector<VectorXd> guess,
                                    const BatchPotential &potential,
                                    const RateInstantonOptions &options) {
  const long N = options.beads;
  if (N < 4 || !(beta > 0.0) || hessSaddle.rows() != saddle.size()) {
    throw std::invalid_argument(
        "optimizeRateInstanton: need N >= 4, beta > 0 and a saddle Hessian "
        "of the saddle's dimension");
  }
  RateInstanton inst;
  inst.beta = beta;
  inst.betaN = beta / static_cast<double>(N);
  inst.temperature = 1.0 / (kBoltzmann * beta);
  inst.crossover = crossoverTemperature(hessSaddle);
  if (!(inst.temperature < inst.crossover)) {
    throw std::invalid_argument(
        "optimizeRateInstanton: T is at or above the crossover temperature; "
        "the ring collapses onto the saddle and steepest descent needs the "
        "parabolic barrier correction, of which classical transition-state "
        "theory is only the one-bead limit");
  }
  bool mirror = options.halfRing && N % 2 == 0;
  if (mirror && static_cast<long>(guess.size()) == N) {
    const long mid = N / 2;
    for (long j = 1; j < mid && mirror; ++j) {
      if ((guess[static_cast<size_t>(j)] - guess[static_cast<size_t>(N - j)])
              .norm() > 1e-8) {
        mirror = false;
      }
    }
  }
  const long active = mirror ? (N / 2 + 1) : N;
  if (options.newtonLimit > 0 &&
      active * saddle.size() <= options.newtonLimit) {
    return optimizeRateByNewton(saddle, hessSaddle, beta, std::move(guess),
                                potential, options);
  }
  const double bnh = inst.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);

  const ColMajorXd saddleCurvature =
      0.5 * (hessSaddle + hessSaddle.transpose());
  const Eigen::SelfAdjointEigenSolver<ColMajorXd> es(saddleCurvature);
  const VectorXd dir = es.eigenvectors().col(0);

  if (static_cast<long>(guess.size()) != N) {
    std::vector<double> v0;
    std::vector<VectorXd> g0;
    potential({saddle}, v0, g0);
    const double vS = v0.at(0);
    // Scan in steps of a quarter of the length over which the barrier's
    // curvature drops V by kB T_c.
    const double h = 0.25 * std::sqrt(2.0 * kBoltzmann * inst.crossover /
                                      -es.eigenvalues()(0));
    const long pts = 200;
    const double dPlus = sideDrop(saddle, dir, vS, 1.0, h, pts, potential);
    const double dMinus = sideDrop(saddle, dir, vS, -1.0, h, pts, potential);
    const double dMin = std::min(dPlus, dMinus);
    const double drop =
        (1.0 - inst.temperature / inst.crossover) *
        (std::isfinite(dMin) ? dMin : kBoltzmann * inst.crossover);
    const double sPlus =
        turningDistance(saddle, dir, vS, drop, 1.0, h, pts, potential);
    const double sMinus =
        turningDistance(saddle, dir, vS, drop, -1.0, h, pts, potential);
    guess.resize(static_cast<size_t>(N));
    for (long j = 0; j < N; ++j) {
      const double ct =
          std::cos(2.0 * std::numbers::pi * static_cast<double>(j) /
                   static_cast<double>(N));
      guess[static_cast<size_t>(j)] =
          saddle + dir * (ct >= 0.0 ? sPlus * ct : sMinus * ct);
    }
  }

  // An even count whose beads already match under j -> N-j is the
  // out-and-back instanton. The images are assigned equal, so the
  // potential on one half is copied, and the step stays on the closed ring.
  const bool wantMirror = options.halfRing && N % 2 == 0;
  std::vector<VectorXd> x = std::move(guess);
  bool fold = false;
  if (wantMirror) {
    fold = true;
    const long m = N / 2;
    for (long j = 1; j < m; ++j) {
      if ((x[static_cast<size_t>(j)] - x[static_cast<size_t>(N - j)]).norm() >
          1e-8) {
        fold = false;
        break;
      }
    }
  }
  BatchPotential evalPot = potential;
  if (fold) {
    evalPot = [&](const std::vector<VectorXd> &q, std::vector<double> &v,
                  std::vector<VectorXd> &g) {
      const long m = N / 2;
      bool sym = true;
      for (long j = 1; j < m; ++j) {
        if ((q[static_cast<size_t>(j)] - q[static_cast<size_t>(N - j)])
                .squaredNorm() != 0.0) {
          sym = false;
          break;
        }
      }
      if (!sym) {
        potential(q, v, g);
        return;
      }
      std::vector<VectorXd> uniq(static_cast<size_t>(m + 1));
      for (long j = 0; j <= m; ++j) {
        uniq[static_cast<size_t>(j)] = q[static_cast<size_t>(j)];
      }
      std::vector<double> vu;
      std::vector<VectorXd> gu;
      potential(uniq, vu, gu);
      v.assign(static_cast<size_t>(N), 0.0);
      g.assign(static_cast<size_t>(N), VectorXd());
      for (long j = 0; j <= m; ++j) {
        v[static_cast<size_t>(j)] = vu[static_cast<size_t>(j)];
        g[static_cast<size_t>(j)] = gu[static_cast<size_t>(j)];
      }
      for (long j = 1; j < m; ++j) {
        v[static_cast<size_t>(N - j)] = vu[static_cast<size_t>(j)];
        g[static_cast<size_t>(N - j)] = gu[static_cast<size_t>(j)];
      }
    };
  }
  auto symmetrize = [&](std::vector<VectorXd> &q) {
    if (!fold) {
      return;
    }
    const long m = N / 2;
    for (long j = 1; j < m; ++j) {
      const size_t a = static_cast<size_t>(j);
      const size_t b = static_cast<size_t>(N - j);
      const VectorXd mid = 0.5 * (q[a] + q[b]);
      q[a] = mid;
      q[b] = mid;
    }
  };
  // The mirror of an out-and-back bounce is symmetric too. A monotone
  // reaction coordinate on beads 0..N/2 leaves that bounce out of the step.
  auto projectFold = [&](std::vector<VectorXd> &q) {
    if (!fold) {
      return;
    }
    projectMonotonePrefix(q, N / 2, dir);
    const long m = N / 2;
    for (long j = 1; j < m; ++j) {
      q[static_cast<size_t>(N - j)] = q[static_cast<size_t>(j)];
    }
  };
  symmetrize(x);
  projectFold(x);

  RingEval cur = evaluateRing(x, c, evalPot, options.energyShift);
  // The unstable mode of the ring starts as every bead moving along the
  // saddle's unstable direction.
  std::vector<VectorXd> mode(x.size(), dir);
  double curvature = lowestMode(x, cur, c, evalPot, mode, options.lanczosFirst,
                                options.lanczosStep, fold);
  std::deque<std::pair<std::vector<VectorXd>, std::vector<VectorXd>>> pairs;
  auto effective = [&](const std::vector<VectorXd> &g) {
    const double par = dot(g, mode);
    std::vector<VectorXd> e = g;
    const double f = curvature < 0.0 ? 2.0 : 1.0;
    for (size_t j = 0; j < e.size(); ++j) {
      e[j] -= f * par * mode[j];
      if (!(curvature < 0.0)) {
        e[j] = -par * mode[j]; // climb along the mode only
      }
    }
    return e;
  };
  std::vector<VectorXd> geff = effective(cur.grad);
  long limit = options.maxIterations;
  for (long it = 0; it < limit; ++it) {
    inst.iterations = it;
    if (curvature < 0.0 && largestBeadNorm(cur.grad) < options.forceTolerance) {
      std::vector<VectorXd> oddMode;
      const double oddCurv =
          fold ? lowestOddMode(x, cur, c, potential, oddMode,
                               options.lanczosFirst, options.lanczosStep)
               : 0.0;
      if (!(oddCurv < -1e-3 * std::abs(curvature))) {
        inst.converged = true;
        break;
      }
      // A mirror-symmetric stationary point with a second unstable mode
      // odd under the mirror, such as two copies of the instanton on one
      // ring. The rest of the search runs on the whole ring from a kick
      // along that mode which lowers beta_N U_N by one in the quadratic
      // model.
      fold = false;
      evalPot = potential;
      const double kick = std::sqrt(2.0 / (inst.betaN * -oddCurv));
      for (size_t k = 0; k < x.size(); ++k) {
        x[k] += kick * oddMode[k];
      }
      cur = evaluateRing(x, c, evalPot, options.energyShift);
      curvature = lowestMode(x, cur, c, evalPot, mode, options.lanczosFirst,
                             options.lanczosStep, fold);
      pairs.clear();
      geff = effective(cur.grad);
      limit = it + options.maxIterations;
    }
    std::vector<VectorXd> trial(x.size());
    std::vector<VectorXd> d = geff;
    std::vector<double> alpha(pairs.size());
    for (size_t i = pairs.size(); i-- > 0;) {
      const double rho = 1.0 / dot(pairs[i].second, pairs[i].first);
      alpha[i] = rho * dot(pairs[i].first, d);
      for (size_t k = 0; k < d.size(); ++k) {
        d[k] -= alpha[i] * pairs[i].second[k];
      }
    }
    double gamma = 1.0 / (4.0 * c); // spring stiffness sets the first scale
    if (!pairs.empty()) {
      gamma = dot(pairs.back().first, pairs.back().second) /
              dot(pairs.back().second, pairs.back().second);
    }
    scale(d, gamma);
    for (size_t i = 0; i < pairs.size(); ++i) {
      const double rho = 1.0 / dot(pairs[i].second, pairs[i].first);
      const double b = rho * dot(pairs[i].second, d);
      for (size_t k = 0; k < d.size(); ++k) {
        d[k] += (alpha[i] - b) * pairs[i].first[k];
      }
    }
    if (!(dot(d, geff) > 0.0)) { // not a descent direction on geff
      pairs.clear();
      d = geff;
      scale(d, 1.0 / (4.0 * c));
    }
    const double big = largestBeadNorm(d);
    if (big > options.maxStep) {
      scale(d, options.maxStep / big);
    }
    for (size_t k = 0; k < x.size(); ++k) {
      trial[k] = x[k] - d[k];
    }
    symmetrize(trial);
    projectFold(trial);
    RingEval next = evaluateRing(trial, c, evalPot, options.energyShift);
    const double prevCurv = curvature;
    const long restart = options.lanczosRestart;
    curvature = lowestMode(trial, next, c, evalPot, mode, restart,
                           options.lanczosStep, fold);
    std::vector<VectorXd> geffNext = effective(next.grad);
    if ((prevCurv < 0.0) != (curvature < 0.0)) {
      pairs.clear();
    } else {
      std::vector<VectorXd> sk(x.size()), yk(x.size());
      for (size_t k = 0; k < x.size(); ++k) {
        sk[k] = trial[k] - x[k];
        yk[k] = geffNext[k] - geff[k];
      }
      if (dot(sk, yk) > 0.0) {
        pairs.emplace_back(std::move(sk), std::move(yk));
        if (static_cast<long>(pairs.size()) > options.memory) {
          pairs.pop_front();
        }
      }
    }
    x = std::move(trial);
    cur = std::move(next);
    geff = std::move(geffNext);
  }
  if (!inst.converged && curvature < 0.0 &&
      largestBeadNorm(cur.grad) < options.forceTolerance) {
    inst.converged = true;
  }
  inst.beads = x;
  inst.energies = cur.v;
  inst.ringPotential = cur.u;
  inst.bN = 0.0;
  for (size_t j = 0; j < x.size(); ++j) {
    inst.bN += (x[(j + 1) % x.size()] - x[j]).squaredNorm();
  }
  return inst;
}

namespace {

/// Flags the count entries of lam nearest zero, among the entries from
/// index `first` on. A saddle passes first = 1 so that its unstable mode,
/// eigenvalue 0 in ascending order, is never taken for a rigid mode however
/// soft the barrier is.
std::vector<bool> nearestZero(const VectorXd &lam, long count, long first = 0) {
  std::vector<long> order(static_cast<size_t>(lam.size() - first));
  std::iota(order.begin(), order.end(), first);
  std::sort(order.begin(), order.end(), [&](long a, long b) {
    return std::abs(lam(a)) < std::abs(lam(b));
  });
  std::vector<bool> out(static_cast<size_t>(lam.size()), false);
  for (long k = 0; k < std::min<long>(count, static_cast<long>(order.size()));
       ++k) {
    out[static_cast<size_t>(order[static_cast<size_t>(k)])] = true;
  }
  return out;
}

} // namespace

void instantonRate(RateInstanton &inst, const RingBeadHessian &hessian,
                   const MatrixXd &hessReactant, double vReactant,
                   const MatrixXd &hessSaddle, double vSaddle, long rigidModes,
                   long denseLimit, const RingRigidBodies &rigidBodies) {
  const long N = static_cast<long>(inst.beads.size());
  if (N < 4 || !(inst.betaN > 0.0)) {
    throw std::invalid_argument("instantonRate: no optimised ring");
  }
  const long f = inst.beads.front().size();
  if (rigidModes < 0 || rigidModes > f) {
    throw std::invalid_argument(
        "instantonRate: rigidModes exceeds the degrees of freedom");
  }
  const double bnh = inst.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);
  const MatrixXd eye = MatrixXd::Identity(f, f);

  std::vector<MatrixXd> hBead(static_cast<size_t>(N));
  std::vector<MatrixXd> diag(static_cast<size_t>(N));
  for (long j = 0; j < N; ++j) {
    const MatrixXd h = hessian(j, inst.beads[static_cast<size_t>(j)]);
    if (h.rows() != f || h.cols() != f) {
      throw std::runtime_error("instantonRate: bead Hessian size");
    }
    hBead[static_cast<size_t>(j)] = 0.5 * (h + h.transpose());
    diag[static_cast<size_t>(j)] =
        hBead[static_cast<size_t>(j)] + 2.0 * c * eye;
  }

  const Eigen::SelfAdjointEigenSolver<MatrixXd> er(
      0.5 * (hessReactant + hessReactant.transpose()));
  const VectorXd &lr = er.eigenvalues();
  const std::vector<bool> rigidR = nearestZero(lr, rigidModes);
  MatrixXd nullBasis(f, 0);
  // The rigid vectors leave the product whether or not every bead Hessian
  // annihilates them exactly; a finite-difference Hessian never does, and
  // the reactant side drops the same modes at k = 0.
  for (long m = 0; m < lr.size(); ++m) {
    if (!rigidR[static_cast<size_t>(m)]) {
      continue;
    }
    nullBasis.conservativeResize(f, nullBasis.cols() + 1);
    nullBasis.col(nullBasis.cols() - 1) = er.eigenvectors().col(m);
  }
  if (rigidModes >= f) {
    throw std::runtime_error("instantonRate: every direction is a rigid mode");
  }
  // det' through the block chain: the cyclic zero mode and the rigid null
  // vectors leave the product by the determinant lemma, the inertia comes
  // from the Schur complements, and the N f by N f matrix is never formed.
  std::vector<VectorXd> cycle(static_cast<size_t>(N));
  double cycleNorm = 0.0;
  for (long j = 0; j < N; ++j) {
    cycle[static_cast<size_t>(j)] =
        0.5 * (inst.beads[static_cast<size_t>((j + 1) % N)] -
               inst.beads[static_cast<size_t>((j + N - 1) % N)]);
    cycleNorm += cycle[static_cast<size_t>(j)].squaredNorm();
  }
  if (!(cycleNorm > 0.0)) {
    throw std::runtime_error(
        "instantonRate: the beads coincide, so the ring has collapsed");
  }
  scale(cycle, 1.0 / std::sqrt(cycleNorm));
  // Each omitted direction is lifted by a spring-sized curvature c, far
  // above any physical near-zero eigenvalue, and the lift comes off the
  // log-determinant again: det(J + c u u^T) = c det' J when J u = 0.
  std::vector<std::vector<VectorXd>> dropped{cycle};
  std::vector<double> kappas{c};
  if (!rigidBodies.sqrtMasses.empty()) {
    // The ring's own rigid motions, orthonormalised against the cycle; the
    // lemma needs an orthonormal basis of the null space, not the
    // generators themselves.
    std::vector<std::vector<VectorXd>> rigid =
        ringRigidBasis(inst.beads, rigidBodies.sqrtMasses,
                       rigidBodies.reference, rigidBodies.rotations);
    if (static_cast<long>(rigid.size()) != nullBasis.cols()) {
      throw std::runtime_error(
          "instantonRate: the ring has " + std::to_string(rigid.size()) +
          " rigid motions and the reactant " +
          std::to_string(nullBasis.cols()));
    }
    for (auto &u : rigid) {
      for (const auto &prev : dropped) {
        const double o = dot(prev, u);
        for (long j = 0; j < N; ++j) {
          u[static_cast<size_t>(j)] -= o * prev[static_cast<size_t>(j)];
        }
      }
      const double un = std::sqrt(dot(u, u));
      if (!(un > 1e-8)) {
        throw std::runtime_error(
            "instantonRate: a rigid motion of the ring lies along its cycle");
      }
      scale(u, 1.0 / un);
      dropped.push_back(std::move(u));
      kappas.push_back(c);
    }
  } else {
    for (long r = 0; r < nullBasis.cols(); ++r) {
      dropped.emplace_back(
          static_cast<size_t>(N),
          (nullBasis.col(r) / std::sqrt(static_cast<double>(N))).eval());
      kappas.push_back(c);
    }
  }
  const WoodburyRing ring(c, diag, true, dropped, kappas, true);
  if (!ring.ok() || !std::isfinite(ring.logAbsDet())) {
    throw std::runtime_error(
        "instantonRate: the ring Hessian is singular and the zero mode was "
        "not removed with the rigid modes");
  }
  inst.zeroEigenvalue = dot(cycle, applyDiagonal(diag, c, true, cycle));
  inst.negativeModes = ring.negative();
  // The lowest ring eigenvalue, for the report, from products alone.
  {
    auto applyFull = [&](const std::vector<VectorXd> &vec) {
      return applyDiagonal(diag, c, true, vec);
    };
    const long dim = N * f;
    const long steps = std::min(dim, static_cast<long>(60));
    std::vector<VectorXd> start = cycle;
    std::uint64_t h = 0x9E3779B97F4A7C15ULL;
    for (auto &bead : start) {
      for (long a2 = 0; a2 < bead.size(); ++a2) {
        h ^= h << 13;
        h ^= h >> 7;
        h ^= h << 17;
        bead(a2) += 0.1 * (static_cast<double>(h >> 11) * 0x1.0p-53 - 0.5);
      }
    }
    const std::vector<RingMode> modes =
        lowestRingModes(applyFull, std::move(start), steps);
    inst.negativeEigenvalue = 0.0;
    for (const auto &mode : modes) {
      if (mode.theta < inst.negativeEigenvalue &&
          std::abs(dot(mode.vector, cycle)) < 0.5) {
        inst.negativeEigenvalue = mode.theta;
      }
    }
    // A numerical null eigenvalue can sit just below zero off the cycle,
    // where the lift along tau does not reach it; it is not a second
    // unstable mode when it is tiny next to the barrier curvature.
    if (inst.negativeModes > 1 && inst.negativeEigenvalue < 0.0) {
      long tiny = 0;
      for (const auto &mode : modes) {
        if (mode.theta < 0.0 && mode.theta > 1e-3 * inst.negativeEigenvalue &&
            std::abs(dot(mode.vector, cycle)) < 0.5) {
          ++tiny;
        }
      }
      inst.negativeModes = std::max(1L, inst.negativeModes - tiny);
    }
  }
  const long nDrop = 1 + nullBasis.cols();
  const double logDetPrime =
      ring.logAbsDet() - static_cast<double>(nDrop) * std::log(c);
  const double logProd =
      static_cast<double>(N * f - nDrop) * std::log(bnh) + 0.5 * logDetPrime;

  inst.logRateTimesZr = -std::log(bnh) +
                        0.5 * std::log(inst.bN / (2.0 * std::numbers::pi *
                                                  inst.betaN * kHbar * kHbar)) -
                        logProd - inst.betaN * inst.ringPotential;

  for (long m = 0; m < lr.size(); ++m) {
    if (!rigidR[static_cast<size_t>(m)] && !(lr(m) > 0.0)) {
      throw std::runtime_error(
          "instantonRate: the reactant Hessian is not positive definite");
    }
  }
  // The rigid modes leave the centroid (k = 0) factor, as the ring's own
  // rigid modes leave its product; both carry them as free particles for
  // k > 0.
  double logZr = -inst.beta * vReactant;
  for (long k = 0; k < N; ++k) {
    const double sk = std::sin(std::numbers::pi * static_cast<double>(k) /
                               static_cast<double>(N));
    for (long m = 0; m < lr.size(); ++m) {
      if (k == 0 && rigidR[static_cast<size_t>(m)]) {
        continue;
      }
      const double l = rigidR[static_cast<size_t>(m)] ? 0.0 : lr(m);
      logZr -= std::log(bnh) + 0.5 * std::log(l + 4.0 * c * sk * sk);
    }
  }
  inst.logZr = logZr;
  inst.logRate = inst.logRateTimesZr - logZr;
  inst.rate = std::exp(inst.logRate) / kTimeUnitSeconds;
  inst.effectiveBarrier =
      -std::log(2.0 * std::numbers::pi * kHbar * inst.beta) / inst.beta -
      inst.logRate / inst.beta;

  if (hessSaddle.size() > 0) {
    inst.classicalLogRate = harmonicTstLogRate(
        hessReactant, hessSaddle, inst.beta, vSaddle - vReactant, rigidModes);
    inst.classicalRate = std::exp(inst.classicalLogRate) / kTimeUnitSeconds;
  }
}

double parabolicFactor(double temperature, double crossover) {
  if (!(temperature > 0.0) || !(crossover > 0.0)) {
    throw std::invalid_argument(
        "parabolicFactor: temperature and crossover must be positive");
  }
  if (!(temperature > crossover)) {
    throw std::invalid_argument(
        "parabolicFactor: T is at or below the crossover; the factor "
        "diverges there");
  }
  // beta hbar omega_b / 2 = pi T_c / T, since T_c = hbar omega_b / (2 pi kB).
  const double phase = std::numbers::pi * crossover / temperature;
  const double s = std::sin(phase);
  if (!(s > 0.0)) {
    throw std::invalid_argument(
        "parabolicFactor: the sine of the barrier phase is not positive");
  }
  return phase / s;
}

double harmonicTstLogRate(const MatrixXd &hessReactant,
                          const MatrixXd &hessSaddle, double beta,
                          double barrier, long rigidModes) {
  if (!(beta > 0.0) || hessReactant.size() == 0 || hessSaddle.size() == 0 ||
      hessReactant.rows() != hessReactant.cols() ||
      hessSaddle.rows() != hessSaddle.cols() ||
      hessReactant.rows() != hessSaddle.rows() || rigidModes < 0) {
    throw std::invalid_argument(
        "harmonicTstLogRate: need beta > 0, matching square Hessians and a "
        "non-negative rigid-mode count");
  }
  const Eigen::SelfAdjointEigenSolver<MatrixXd> er(
      0.5 * (hessReactant + hessReactant.transpose()), Eigen::EigenvaluesOnly);
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(
      0.5 * (hessSaddle + hessSaddle.transpose()), Eigen::EigenvaluesOnly);
  const VectorXd &lr = er.eigenvalues();
  const VectorXd &ls = es.eigenvalues();
  if (!(ls(0) < 0.0)) {
    throw std::invalid_argument(
        "harmonicTstLogRate: the saddle Hessian has no negative eigenvalue");
  }
  const std::vector<bool> rigidR = nearestZero(lr, rigidModes);
  const std::vector<bool> rigidS = nearestZero(ls, rigidModes, 1);
  double logRatio = 0.0;
  for (long m = 0; m < lr.size(); ++m) {
    if (!rigidR[static_cast<size_t>(m)]) {
      logRatio += 0.5 * std::log(lr(m));
    }
  }
  // Eigenvalue 0 is the unstable mode, the most negative.
  for (long m = 1; m < ls.size(); ++m) {
    if (!rigidS[static_cast<size_t>(m)]) {
      logRatio -= 0.5 * std::log(std::abs(ls(m)));
    }
  }
  return logRatio - std::log(2.0 * std::numbers::pi) - beta * barrier;
}

double quantumHarmonicTstLogRate(const MatrixXd &hessReactant,
                                 const MatrixXd &hessSaddle, double beta,
                                 double barrier, long rigidModes) {
  if (!(beta > 0.0) || hessReactant.size() == 0 || hessSaddle.size() == 0 ||
      hessReactant.rows() != hessReactant.cols() ||
      hessSaddle.rows() != hessSaddle.cols() ||
      hessReactant.rows() != hessSaddle.rows() || rigidModes < 0) {
    throw std::invalid_argument(
        "quantumHarmonicTstLogRate: need beta > 0, matching square Hessians "
        "and a non-negative rigid-mode count");
  }
  const ColMajorXd hr = 0.5 * (hessReactant + hessReactant.transpose());
  const ColMajorXd hsd = 0.5 * (hessSaddle + hessSaddle.transpose());
  const Eigen::SelfAdjointEigenSolver<ColMajorXd> er(hr,
                                                     Eigen::EigenvaluesOnly);
  const Eigen::SelfAdjointEigenSolver<ColMajorXd> es(hsd,
                                                     Eigen::EigenvaluesOnly);
  const VectorXd lr = er.eigenvalues();
  const VectorXd ls = es.eigenvalues();
  if (!(ls(0) < 0.0)) {
    throw std::invalid_argument("quantumHarmonicTstLogRate: the saddle "
                                "Hessian has no negative eigenvalue");
  }
  const std::vector<bool> rigidR = nearestZero(lr, rigidModes);
  const std::vector<bool> rigidS = nearestZero(ls, rigidModes, 1);
  const double bh = beta * kHbar;
  // ln(2 sinh(x / 2)) without overflow for large x.
  auto logTwoSinhHalf = [](double x) {
    return 0.5 * x + std::log1p(-std::exp(-x));
  };
  double logRatio = 0.0;
  for (long m = 0; m < lr.size(); ++m) {
    if (!rigidR[static_cast<size_t>(m)]) {
      logRatio += logTwoSinhHalf(bh * std::sqrt(lr(m)));
    }
  }
  for (long m = 1; m < ls.size(); ++m) {
    if (!rigidS[static_cast<size_t>(m)]) {
      logRatio -= logTwoSinhHalf(bh * std::sqrt(std::abs(ls(m))));
    }
  }
  return logRatio - std::log(2.0 * std::numbers::pi * bh) - beta * barrier;
}

} // namespace eonc::tunneling

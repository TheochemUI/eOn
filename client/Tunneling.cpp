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
#include <deque>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>

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
  const BlockChain chain(c, diag);

  auto wellLogDet = [&](const MatrixXd &h) {
    const std::vector<MatrixXd> d(static_cast<size_t>(P - 1),
                                  spring + dtau * 0.5 * (h + h.transpose()));
    const BlockChain well(c, d);
    if (well.sign() < 0) {
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
  const int signPrime = chain.sign() * (vJv < 0.0 ? -1 : 1);
  if (signPrime < 0) {
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
  (void)mirror;
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
    return out;
  };
  std::vector<std::vector<VectorXd>> basis;
  std::vector<double> alpha, beta;
  std::vector<VectorXd> q = mode;
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
  mode = std::move(ritz);
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
  return std::numeric_limits<double>::infinity(); // still falling
}

// sum_{k=1}^{N-1} log(4 c sin^2(pi k / N)). The product of those sines is N,
// so the sum is 2 log N + (N - 1) log c.
double flatSpringLog(long n, double c) {
  return 2.0 * std::log(static_cast<double>(n)) +
         static_cast<double>(n - 1) * std::log(c);
}

// Orthonormal complement of a full-rank thin basis. The leading columns of
// the Householder Q span the basis, and the rest are orthogonal to it.
MatrixXd complementOf(const MatrixXd &nullBasis) {
  const long f = nullBasis.rows();
  const long k = nullBasis.cols();
  const Eigen::HouseholderQR<MatrixXd> qr(nullBasis);
  const MatrixXd q = qr.householderQ() * MatrixXd::Identity(f, f);
  return q.rightCols(f - k);
}

std::vector<MatrixXd> congruences(const std::vector<MatrixXd> &diag,
                                  const MatrixXd &cred) {
  std::vector<MatrixXd> out;
  out.reserve(diag.size());
  const MatrixXd ct = cred.transpose();
  for (const auto &block : diag) {
    out.push_back(ct * block * cred);
  }
  return out;
}

bool logAbsAgrees(double a, double b) {
  if (!std::isfinite(a) || !std::isfinite(b)) {
    return false;
  }
  const double d = std::abs(a - b);
  return d <= 1e-6 || d <= 1e-6 * std::max(1.0, std::abs(a));
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

struct OmittedMode {
  double theta = 0.0;
  std::vector<VectorXd> vector;
};

void projectOut(std::vector<VectorXd> &v,
                const std::vector<std::vector<VectorXd>> &locked) {
  for (const auto &p : locked) {
    const double s = dot(v, p);
    for (size_t j = 0; j < v.size(); ++j) {
      v[j] -= s * p[j];
    }
  }
}

// Eigenvalues closest to zero, by inverse iteration on a factored ring.
// `apply` is that same operator, used for the Rayleigh quotient.
std::vector<OmittedMode> modesClosestToZero(
    const CyclicFactor &fac,
    const std::function<std::vector<VectorXd>(const std::vector<VectorXd> &)>
        &apply,
    const std::vector<VectorXd> &seed, long count) {
  const long n = static_cast<long>(seed.size());
  const long f = seed.empty() ? 0 : seed.front().size();
  if (count < 1 || f < 1 || fac.singular) {
    throw std::runtime_error(
        "instantonRate: the cyclic zero mode is not resolved");
  }
  std::vector<std::vector<VectorXd>> locked;
  std::vector<OmittedMode> out;
  out.reserve(static_cast<size_t>(count));
  for (long m = 0; m < count; ++m) {
    std::vector<VectorXd> v;
    if (m == 0) {
      v = seed;
    } else {
      v.assign(static_cast<size_t>(n), VectorXd::Zero(f));
      v[static_cast<size_t>(m % n)](m % f) = 1.0;
    }
    double theta = 0.0;
    for (int it = 0; it < 40; ++it) {
      projectOut(v, locked);
      const double nv = std::sqrt(dot(v, v));
      if (!(nv > 1e-14)) {
        throw std::runtime_error(
            "instantonRate: the cyclic zero mode is not resolved");
      }
      scale(v, 1.0 / nv);
      std::vector<VectorXd> y = fac.solve(v);
      projectOut(y, locked);
      const double ny = std::sqrt(dot(y, y));
      if (!(ny > 0.0) || !std::isfinite(ny)) {
        throw std::runtime_error(
            "instantonRate: the cyclic zero mode is not resolved");
      }
      scale(y, 1.0 / ny);
      const std::vector<VectorXd> hy = apply(y);
      theta = dot(y, hy);
      v = std::move(y);
    }
    if (!std::isfinite(theta)) {
      throw std::runtime_error(
          "instantonRate: the cyclic zero mode is not resolved");
    }
    OmittedMode mode;
    mode.theta = theta;
    mode.vector = std::move(v);
    locked.push_back(mode.vector);
    out.push_back(std::move(mode));
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

// Closed ring: diagonal H + 2 c I, neighbour coupling -c I, including the
// corner. Each edge is written once.
ColMajorXd closedRingMatrix(const std::vector<MatrixXd> &physical, double c) {
  const long nBeads = static_cast<long>(physical.size());
  const long f = physical.front().rows();
  ColMajorXd big = ColMajorXd::Zero(nBeads * f, nBeads * f);
  const ColMajorXd eye = ColMajorXd::Identity(f, f);
  for (long j = 0; j < nBeads; ++j) {
    ColMajorXd block = 0.5 * (physical[static_cast<size_t>(j)] +
                              physical[static_cast<size_t>(j)].transpose());
    block += 2.0 * c * eye;
    big.block(j * f, j * f, f, f) = block;
  }
  for (long j = 0; j < nBeads; ++j) {
    const long k = (j + 1) % nBeads;
    if (j < k) {
      big.block(j * f, k * f, f, f) = -c * eye;
      big.block(k * f, j * f, f, f) = -c * eye;
    }
  }
  if (nBeads > 1) {
    const long last = nBeads - 1;
    big.block(last * f, 0, f, f) = -c * eye;
    big.block(0, last * f, f, f) = -c * eye;
  }
  return big;
}

// Half chain from one turning point to the other. End blocks hold half the
// physical Hessian, interior blocks the whole of it. Neighbour coupling is
// -c I. On a symmetric ring this is the Hessian of half the closed-ring energy.
ColMajorXd halfRingMatrix(const std::vector<MatrixXd> &physical, double c) {
  const long beads = static_cast<long>(physical.size());
  const long f = physical.front().rows();
  ColMajorXd big = ColMajorXd::Zero(beads * f, beads * f);
  const ColMajorXd eye = ColMajorXd::Identity(f, f);
  for (long j = 0; j < beads; ++j) {
    const bool end = j == 0 || j + 1 == beads;
    ColMajorXd block = 0.5 * (physical[static_cast<size_t>(j)] +
                              physical[static_cast<size_t>(j)].transpose());
    if (end) {
      block *= 0.5;
    }
    block += (end ? c : 2.0 * c) * eye;
    big.block(j * f, j * f, f, f) = block;
    if (j + 1 < beads) {
      big.block(j * f, (j + 1) * f, f, f) = -c * eye;
      big.block((j + 1) * f, j * f, f, f) = -c * eye;
    }
  }
  return big;
}

struct RingSpectrum {
  VectorXd values;
  ColMajorXd vectors;
  ColMajorXd ring;
  bool ok = false;
};

RingSpectrum spectrumOf(ColMajorXd ring) {
  RingSpectrum out;
  const Eigen::SelfAdjointEigenSolver<ColMajorXd> es(ring);
  out.ok = es.info() == Eigen::Success && ring.array().isFinite().all();
  if (out.ok) {
    out.values = es.eigenvalues();
    out.vectors = es.eigenvectors();
    out.ring = std::move(ring);
  }
  return out;
}

struct Climb {
  long index = -1;
  double curvature = 0.0;
  long negative = 0;
};

Climb classify(const RingSpectrum &sp, const VectorXd &tau, double spring) {
  Climb out;
  const double cut = -1e-8 * std::max(1.0, spring);
  double lowest = std::numeric_limits<double>::infinity();
  for (long i = 0; i < sp.values.size(); ++i) {
    if (overlapsTau(sp.vectors.col(i), tau)) {
      continue;
    }
    if (sp.values(i) < cut) {
      ++out.negative;
    }
    if (sp.values(i) < lowest) {
      lowest = sp.values(i);
      out.index = i;
      out.curvature = sp.values(i);
    }
  }
  return out;
}

// Leave a negative climb eigenvalue in place: a raw Newton step climbs it.
// Flip every other negative eigenvalue. The imaginary-time cycle is left out,
// and a shift along it keeps the solve off that direction.
VectorXd indexOneStep(const RingSpectrum &sp, const Climb &climb,
                      const VectorXd &gflat, const VectorXd &tau,
                      double spring) {
  if (!sp.ok || climb.index < 0 || gflat.size() != sp.values.size()) {
    return VectorXd();
  }
  const double cut = -1e-8 * std::max(1.0, spring);
  ColMajorXd jt = sp.ring;
  if (tau.size() == gflat.size()) {
    jt.noalias() += spring * (tau * tau.transpose());
  }
  for (long i = 0; i < sp.values.size(); ++i) {
    if (overlapsTau(sp.vectors.col(i), tau)) {
      continue;
    }
    const double li = sp.values(i);
    const bool flip =
        (i == climb.index && li > 0.0) || (i != climb.index && li < cut);
    if (!flip) {
      continue;
    }
    const VectorXd v = sp.vectors.col(i);
    jt.noalias() -= (2.0 * li) * (v * v.transpose());
  }
  jt = (0.5 * (jt + jt.transpose())).eval();
  const Eigen::PartialPivLU<ColMajorXd> lu(jt);
  VectorXd step = lu.solve(-gflat);
  const double rhs = std::max(1.0, gflat.norm());
  if (!step.array().isFinite().all() ||
      (jt * step + gflat).norm() > 1e-6 * rhs) {
    return VectorXd();
  }
  if (tau.size() == step.size()) {
    step -= step.dot(tau) * tau;
  }
  if (!step.array().isFinite().all()) {
    return VectorXd();
  }
  return step;
}

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
                      (std::isfinite(dMin) ? dMin : kBoltzmann * crossover);
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

struct NewtonOut {
  std::vector<VectorXd> beads;
  std::vector<double> energies;
  double ringPotential = 0.0;
  double bN = 0.0;
  long iterations = 0;
  bool converged = false;
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
  std::vector<MatrixXd> physical(x.size(), hS);

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
    RingSpectrum sp;
    VectorXd tau;
    Climb climb;
    double gmax = 0.0;
  };
  auto viewOf = [&](const Obj &ev) {
    View v;
    v.tau = half ? VectorXd() : timeTranslation(x);
    v.sp = spectrumOf(half ? halfRingMatrix(physical, c)
                           : closedRingMatrix(physical, c));
    if (v.sp.ok) {
      v.climb = classify(v.sp, v.tau, c);
      v.gmax = largestBeadNorm(ev.grad);
    }
    return v;
  };
  auto done = [&](const View &v) {
    return v.sp.ok && v.gmax < options.forceTolerance &&
           v.climb.negative == 1 && v.climb.curvature < 0.0;
  };

  Obj cur = objective(x);
  double trust = options.maxStep;
  long entries = 0;
  bool converged = false;
  const double trustFloor = std::min(1e-4, options.maxStep);
  for (long it = 0; it < options.maxIterations; ++it) {
    ++entries;
    const View v = viewOf(cur);
    if (done(v)) {
      converged = true;
      break;
    }
    VectorXd step = indexOneStep(v.sp, v.climb, packBeads(cur.grad), v.tau, c);
    if (step.size() == 0) {
      trust = std::max(0.5 * trust, trustFloor);
      continue;
    }
    const double big = packedBeadNorm(step, f);
    if (big > trust) {
      step *= trust / big;
    }
    const VectorXd gflat = packBeads(cur.grad);
    const double pred = gflat.dot(step) + 0.5 * step.dot(v.sp.ring * step);
    std::vector<VectorXd> trial = x;
    addPacked(trial, step);
    if (!finiteBeads(trial)) {
      trust = std::max(0.5 * trust, trustFloor);
      continue;
    }
    Obj next = objective(trial);
    const double ratio = std::abs(pred) > 1e-30 ? (next.u - cur.u) / pred : 1.0;
    if (!std::isfinite(ratio) || ratio < 0.1 || ratio > 3.0 ||
        !finiteBeads(next.grad) || !std::isfinite(next.u)) {
      trust = std::max(0.5 * trust, trustFloor);
      continue;
    }
    for (size_t k = 0; k < x.size(); ++k) {
      bofillUpdate(physical[k], trial[k] - x[k],
                   next.gradPot[k] - cur.gradPot[k]);
    }
    x = std::move(trial);
    cur = std::move(next);
    if (ratio > 0.75 && ratio < 1.25 && big >= 0.99 * trust) {
      trust = std::min(2.0 * trust, options.maxStep);
    }
  }
  if (!converged && done(viewOf(cur))) {
    converged = true;
  }

  NewtonOut out;
  out.iterations = entries;
  out.converged = converged;
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

// Below 0.75 Tc an empty guess walks down from 0.85 Tc. Each stage passes a
// full-length ring onward, so the walk is not repeated on the way back in.
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
  const NewtonOut got =
      newtonInstanton(std::move(guess), spring, hessSaddle, options, potential);
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
  if (active * saddle.size() <= 4096) {
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
  // out-and-back instanton. The potential on one half is copied onto the
  // other. The step stays in the closed-ring product, because the
  // derivative on an interior bead is the sum of both images.
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
        // A few ulps: the image force matches, and the log rate stays put.
        // A looser threshold of 1e-10 moves the LJ13 rate.
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
  for (long it = 0; it < options.maxIterations; ++it) {
    inst.iterations = it;
    if (curvature < 0.0 && largestBeadNorm(cur.grad) < options.forceTolerance) {
      inst.converged = true;
      break;
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
    RingEval next = evaluateRing(trial, c, evalPot, options.energyShift);
    const double prevCurv = curvature;
    const long restart = fold ? 4L : options.lanczosRestart;
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

/// Flags the count entries of lam nearest zero.
std::vector<bool> nearestZero(const VectorXd &lam, long count) {
  std::vector<long> order(static_cast<size_t>(lam.size()));
  std::iota(order.begin(), order.end(), 0L);
  std::sort(order.begin(), order.end(), [&](long a, long b) {
    return std::abs(lam(a)) < std::abs(lam(b));
  });
  std::vector<bool> out(static_cast<size_t>(lam.size()), false);
  for (long k = 0; k < std::min<long>(count, lam.size()); ++k) {
    out[static_cast<size_t>(order[static_cast<size_t>(k)])] = true;
  }
  return out;
}

} // namespace

void instantonRate(RateInstanton &inst, const RingBeadHessian &hessian,
                   const MatrixXd &hessReactant, double vReactant,
                   const MatrixXd &hessSaddle, double vSaddle, long rigidModes,
                   long denseLimit) {
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
  const double zeroCut = 1e-8 * c;
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
  bool exactRigid = rigidModes > 0;
  for (long m = 0; m < lr.size(); ++m) {
    if (!rigidR[static_cast<size_t>(m)]) {
      continue;
    }
    if (!(std::abs(lr(m)) < zeroCut)) {
      exactRigid = false;
    }
    nullBasis.conservativeResize(f, nullBasis.cols() + 1);
    nullBasis.col(nullBasis.cols() - 1) = er.eigenvectors().col(m);
  }
  if (exactRigid && rigidModes > 0) {
    for (long j = 0; j < N; ++j) {
      const double residual =
          (hBead[static_cast<size_t>(j)] * nullBasis).norm();
      if (residual > std::sqrt(zeroCut) * static_cast<double>(rigidModes)) {
        exactRigid = false;
        break;
      }
    }
  }
  if (!exactRigid) {
    nullBasis.resize(f, 0);
  } else if (rigidModes >= f) {
    throw std::runtime_error("instantonRate: every direction is a rigid mode");
  }

  auto flatFraction = [&](const std::vector<VectorXd> &vec) {
    if (nullBasis.cols() == 0) {
      return 0.0;
    }
    double flat = 0.0;
    double total = 0.0;
    for (const auto &bead : vec) {
      total += bead.squaredNorm();
      flat += (nullBasis.transpose() * bead).squaredNorm();
    }
    return total > 0.0 ? flat / total : 0.0;
  };

  double logProd = 0.0;
  const bool useDense =
      denseLimit < 0 || (denseLimit > 0 && N * f <= denseLimit);
  if (useDense) {
    MatrixXd big = MatrixXd::Zero(N * f, N * f);
    for (long j = 0; j < N; ++j) {
      big.block(j * f, j * f, f, f) = diag[static_cast<size_t>(j)];
      const long k = (j + 1) % N;
      big.block(j * f, k * f, f, f) -= c * eye;
      big.block(k * f, j * f, f, f) -= c * eye;
    }
    // MatrixXd is row-major. The self-adjoint solver reads a column-major
    // triangle, so the ring matrix is copied before the decomposition.
    const ColMajorXd ring = big;
    const Eigen::SelfAdjointEigenSolver<ColMajorXd> es(ring,
                                                       Eigen::EigenvaluesOnly);
    const VectorXd lam = es.eigenvalues();
    if (exactRigid) {
      std::vector<long> order(static_cast<size_t>(lam.size()));
      std::iota(order.begin(), order.end(), 0L);
      std::sort(order.begin(), order.end(), [&](long a, long b) {
        return std::abs(lam(a)) < std::abs(lam(b));
      });
      double denseKept = 0.0;
      for (long k = rigidModes; k < lam.size(); ++k) {
        denseKept += std::log(std::abs(lam(order[static_cast<size_t>(k)])));
      }
      const double blockLog =
          cyclicRingLogAbsDet(c, congruences(diag, complementOf(nullBasis))) +
          static_cast<double>(rigidModes) * flatSpringLog(N, c);
      if (!logAbsAgrees(denseKept, blockLog)) {
        throw std::runtime_error(
            "instantonRate: the reduced ring determinant disagrees with the "
            "dense product (" +
            std::to_string(blockLog) + " against " + std::to_string(denseKept) +
            ")");
      }
    } else {
      double denseLog = 0.0;
      for (long k = 0; k < lam.size(); ++k) {
        denseLog += std::log(std::abs(lam(k)));
      }
      // A numerical null mode makes log|det| disagree by tens of nats while
      // the modes kept in the rate still match. Check only a definite ring.
      const double blockLog = cyclicRingLogAbsDet(c, diag);
      if (lam.cwiseAbs().minCoeff() > zeroCut &&
          !logAbsAgrees(blockLog, denseLog)) {
        throw std::runtime_error(
            "instantonRate: the cyclic block determinant disagrees with the "
            "dense ring Hessian (" +
            std::to_string(blockLog) + " against " + std::to_string(denseLog) +
            ")");
      }
    }
    // The zero mode (the ring's translation in imaginary time) and the rigid
    // modes are the 1 + rigidModes eigenvalues nearest zero.
    const std::vector<bool> dropped = nearestZero(lam, 1 + rigidModes);
    inst.zeroEigenvalue = 0.0;
    for (long k = 0; k < lam.size(); ++k) {
      if (dropped[static_cast<size_t>(k)] &&
          std::abs(lam(k)) > std::abs(inst.zeroEigenvalue)) {
        inst.zeroEigenvalue = lam(k);
      }
    }
    inst.negativeModes = 0;
    inst.negativeEigenvalue = 0.0;
    for (long k = 0; k < lam.size(); ++k) {
      if (dropped[static_cast<size_t>(k)]) {
        continue;
      }
      if (lam(k) < 0.0) {
        ++inst.negativeModes;
        inst.negativeEigenvalue = std::min(inst.negativeEigenvalue, lam(k));
      }
      logProd += std::log(bnh) + 0.5 * std::log(std::abs(lam(k)));
    }
  } else {
    const long nDrop = 1 + rigidModes;
    const long nMore = exactRigid ? 1 : nDrop;
    const MatrixXd cred = exactRigid ? complementOf(nullBasis) : MatrixXd();
    const std::vector<MatrixXd> blocks =
        exactRigid ? congruences(diag, cred) : diag;
    const CyclicFactor fac(c, blocks);
    double logAbs = fac.logAbs;
    if (exactRigid) {
      logAbs += static_cast<double>(rigidModes) * flatSpringLog(N, c);
    }
    if (!std::isfinite(logAbs)) {
      throw std::runtime_error(
          "instantonRate: the ring Hessian is singular and the zero mode was "
          "not removed with the rigid modes");
    }
    std::vector<VectorXd> cycleFull(static_cast<size_t>(N));
    double cycleNorm = 0.0;
    for (long j = 0; j < N; ++j) {
      cycleFull[static_cast<size_t>(j)] =
          inst.beads[static_cast<size_t>((j + 1) % N)] -
          inst.beads[static_cast<size_t>((j + N - 1) % N)];
      cycleNorm += cycleFull[static_cast<size_t>(j)].squaredNorm();
    }
    if (!(cycleNorm > 0.0)) {
      throw std::runtime_error(
          "instantonRate: the beads coincide, so the ring has collapsed");
    }
    scale(cycleFull, 1.0 / std::sqrt(cycleNorm));
    std::vector<VectorXd> cycle = cycleFull;
    if (exactRigid) {
      for (long j = 0; j < N; ++j) {
        cycle[static_cast<size_t>(j)] =
            cred.transpose() * cycleFull[static_cast<size_t>(j)];
      }
      const double reducedNorm = std::sqrt(dot(cycle, cycle));
      if (!(reducedNorm > 1e-8)) {
        throw std::runtime_error(
            "instantonRate: the cyclic zero mode is not resolved");
      }
      scale(cycle, 1.0 / reducedNorm);
    }
    auto applyFull = [&](const std::vector<VectorXd> &vec) {
      std::vector<VectorXd> out(static_cast<size_t>(N));
      for (long j = 0; j < N; ++j) {
        const long prev = (j + N - 1) % N;
        const long next = (j + 1) % N;
        out[static_cast<size_t>(j)] =
            hBead[static_cast<size_t>(j)] * vec[static_cast<size_t>(j)] +
            c * (2.0 * vec[static_cast<size_t>(j)] -
                 vec[static_cast<size_t>(prev)] -
                 vec[static_cast<size_t>(next)]);
      }
      return out;
    };
    auto applyBlocks = [&](const std::vector<VectorXd> &vec) {
      if (!exactRigid) {
        return applyFull(vec);
      }
      std::vector<VectorXd> full(static_cast<size_t>(N));
      for (long j = 0; j < N; ++j) {
        full[static_cast<size_t>(j)] = cred * vec[static_cast<size_t>(j)];
      }
      const std::vector<VectorXd> acted = applyFull(full);
      std::vector<VectorXd> out(static_cast<size_t>(N));
      for (long j = 0; j < N; ++j) {
        out[static_cast<size_t>(j)] =
            cred.transpose() * acted[static_cast<size_t>(j)];
      }
      return out;
    };
    // The product divides by these eigenvalues. Inverse iteration on the
    // factored ring resolves a near-zero mode; Lanczos of an interior
    // eigenvalue does not, once the ring is long.
    const std::vector<OmittedMode> omitted =
        modesClosestToZero(fac, applyBlocks, cycle, nMore);
    inst.zeroEigenvalue = 0.0;
    double bestOverlap = 0.0;
    for (const auto &mode : omitted) {
      const double overlap = std::abs(dot(mode.vector, cycle));
      if (overlap >= bestOverlap) {
        bestOverlap = overlap;
        inst.zeroEigenvalue = mode.theta;
      }
      logAbs -= std::log(std::abs(mode.theta));
    }
    if (!(bestOverlap > 0.5) || !std::isfinite(logAbs)) {
      throw std::runtime_error(
          "instantonRate: the cyclic zero mode is not resolved (overlap " +
          std::to_string(bestOverlap) + ")");
    }
    // A negative Ritz value along the bead velocity is the cyclic zero,
    // already omitted. Any other resolved negative value is an extra
    // unstable mode.
    const double residualCut = 1e-4 * c;
    const long dim = N * f;
    const long steps = dim <= 1024 ? dim : std::min(dim, static_cast<long>(80));
    std::vector<VectorXd> start = cycleFull;
    start.front()(0) += 0.1;
    const std::vector<RingMode> modes =
        lowestRingModes(applyFull, std::move(start), steps);
    inst.negativeModes = 0;
    inst.negativeEigenvalue = 0.0;
    double unresolved = 0.0;
    std::vector<double> negative;
    for (const auto &mode : modes) {
      if (!(mode.theta < 0.0)) {
        continue;
      }
      const double overlap = std::abs(dot(mode.vector, cycleFull));
      if (overlap > 0.5 || flatFraction(mode.vector) > 0.5) {
        continue;
      }
      if (mode.residual > residualCut) {
        unresolved = std::max(unresolved, mode.residual);
        continue;
      }
      negative.push_back(mode.theta);
      inst.negativeEigenvalue = std::min(inst.negativeEigenvalue, mode.theta);
    }
    // A numerical null eigenvalue can sit just below zero. It is not a
    // second unstable mode when it is tiny next to the barrier curvature.
    for (const double theta : negative) {
      if (theta <= 1e-3 * inst.negativeEigenvalue) {
        ++inst.negativeModes;
      }
    }
    if (inst.negativeModes == 0 && unresolved > 0.0) {
      throw std::runtime_error(
          "instantonRate: the negative ring mode is not resolved (residual " +
          std::to_string(unresolved) + ")");
    }
    const long nKept = N * f - nDrop;
    logProd = static_cast<double>(nKept) * std::log(bnh) + 0.5 * logAbs;
  }

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
    const Eigen::SelfAdjointEigenSolver<MatrixXd> ets(
        0.5 * (hessSaddle + hessSaddle.transpose()), Eigen::EigenvaluesOnly);
    const VectorXd &ls = ets.eigenvalues();
    const std::vector<bool> rigidS = nearestZero(ls, rigidModes);
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
    inst.classicalLogRate = logRatio - std::log(2.0 * std::numbers::pi) -
                            inst.beta * (vSaddle - vReactant);
    inst.classicalRate = std::exp(inst.classicalLogRate) / kTimeUnitSeconds;
  }
}

} // namespace eonc::tunneling

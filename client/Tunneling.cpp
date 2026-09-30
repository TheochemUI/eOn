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
#include "eon/Tunneling.h"

#include <Eigen/Eigenvalues>

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

// J = c T (x) I + blockdiag(A_k), T the tridiagonal (2, -1) spring matrix;
// the diagonal blocks hold 2c already, the off-diagonal blocks are -c I.
// Block LU: D_1 = A_1, D_k = A_k - c^2 D_{k-1}^-1.
class BlockChain {
public:
  BlockChain(double c, const std::vector<MatrixXd> &diag)
      : c_(c) {
    lu_.reserve(diag.size());
    for (size_t k = 0; k < diag.size(); ++k) {
      MatrixXd d = diag[k];
      if (k > 0) {
        d -= c_ * c_ * lu_.back().inverse();
      }
      lu_.emplace_back(d);
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

private:
  double c_;
  std::vector<Eigen::PartialPivLU<MatrixXd>> lu_;
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
                      const BatchPotential &potential) {
  RingEval out;
  std::vector<VectorXd> gv;
  potential(x, out.v, gv);
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
                  std::vector<VectorXd> &mode, long steps, double eps) {
  const size_t n = x.size();
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
        "the ring collapses onto the saddle and classical transition-state "
        "theory applies");
  }
  const double bnh = inst.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);

  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(
      0.5 * (hessSaddle + hessSaddle.transpose()));
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

  std::vector<VectorXd> x = std::move(guess);
  RingEval cur = evaluateRing(x, c, potential);
  // The unstable mode of the ring starts as every bead moving along the
  // saddle's unstable direction.
  std::vector<VectorXd> mode(x.size(), dir);
  double curvature = lowestMode(x, cur, c, potential, mode,
                                options.lanczosFirst, options.lanczosStep);
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
    std::vector<VectorXd> trial(x.size());
    for (size_t k = 0; k < x.size(); ++k) {
      trial[k] = x[k] - d[k];
    }
    RingEval next = evaluateRing(trial, c, potential);
    const double prevCurv = curvature;
    curvature = lowestMode(trial, next, c, potential, mode,
                           options.lanczosRestart, options.lanczosStep);
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
                   const MatrixXd &hessSaddle, double vSaddle,
                   long rigidModes) {
  const long N = static_cast<long>(inst.beads.size());
  if (N < 4 || !(inst.betaN > 0.0)) {
    throw std::invalid_argument("instantonRate: no optimised ring");
  }
  const long f = inst.beads.front().size();
  const double bnh = inst.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);
  MatrixXd big = MatrixXd::Zero(N * f, N * f);
  for (long j = 0; j < N; ++j) {
    const MatrixXd h = hessian(j, inst.beads[static_cast<size_t>(j)]);
    if (h.rows() != f || h.cols() != f) {
      throw std::runtime_error("instantonRate: bead Hessian size");
    }
    big.block(j * f, j * f, f, f) =
        0.5 * (h + h.transpose()) + 2.0 * c * MatrixXd::Identity(f, f);
    const long k = (j + 1) % N;
    big.block(j * f, k * f, f, f) -= c * MatrixXd::Identity(f, f);
    big.block(k * f, j * f, f, f) -= c * MatrixXd::Identity(f, f);
  }
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(big, Eigen::EigenvaluesOnly);
  const VectorXd &lam = es.eigenvalues();
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
  double logProd = 0.0;
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
  inst.logRateTimesZr = -std::log(bnh) +
                        0.5 * std::log(inst.bN / (2.0 * std::numbers::pi *
                                                  inst.betaN * kHbar * kHbar)) -
                        logProd - inst.betaN * inst.ringPotential;

  const Eigen::SelfAdjointEigenSolver<MatrixXd> er(
      0.5 * (hessReactant + hessReactant.transpose()), Eigen::EigenvaluesOnly);
  const VectorXd &lr = er.eigenvalues();
  const std::vector<bool> rigidR = nearestZero(lr, rigidModes);
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

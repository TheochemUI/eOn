/*
** This file is part of eOn.
**
** i-PI Copyright (C) 2014-2015 i-PI developers
**
** Permission is hereby granted, free of charge, to any person obtaining
** a copy of this software and associated documentation files (the
** "Software"), to deal in the Software without restriction, including
** without limitation the rights to use, copy, modify, merge, publish,
** distribute, sublicense, and/or sell copies of the Software, and to
** permit persons to whom the Software is furnished to do so, subject to
** the following conditions:
**
** The above copyright notice and this permission notice shall be
** included in all copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
** MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
** NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
** BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
** ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
** CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
** SOFTWARE.
**
** SPDX-License-Identifier: MIT
*/

#include "eon/PathIntegral.h"

#include <unsupported/Eigen/MatrixFunctions>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace eonc::pathintegral {
namespace {

constexpr double kPi = std::numbers::pi;

std::string lower(std::string s) {
  for (char &c : s) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return s;
}

double ecoResponse(double x) {
  const double z = 0.5 * x;
  if (z < 0.25) {
    const double z2 = z * z;
    const double denom =
        1.0 / 3.0 - z2 / 45.0 + 2.0 * z2 * z2 / 945.0 - z2 * z2 * z2 / 4725.0;
    return 4.0 / denom;
  }
  return (x * x) / (z / std::tanh(z) - 1.0);
}

VectorXd ecoFit(long nBeads, double xmax) {
  const long nfree = nBeads / 2;
  VectorXd mult(nfree);
  mult.setConstant(2.0);
  if (nBeads % 2 == 0) {
    mult[nfree - 1] = 1.0;
  }
  const long m = std::max(std::lround(10.0 * xmax), 100L);
  VectorXd x(m);
  VectorXd f(m);
  for (long i = 0; i < m; ++i) {
    x[i] = (static_cast<double>(i) + 0.5) * (xmax / static_cast<double>(m));
    f[i] = ecoResponse(x[i]);
  }

  auto objective = [&](const VectorXd &y, double &s, VectorXd &g, MatrixXd &h) {
    MatrixXd d(m, nfree);
    MatrixXd e(m, nfree);
    MatrixXd dg(m, nfree);
    MatrixXd d2(m, nfree);
    VectorXd r(m);
    for (long i = 0; i < m; ++i) {
      double row = 0.0;
      for (long k = 0; k < nfree; ++k) {
        d(i, k) = 1.0 / (y[k] * y[k] + x[i] * x[i]);
        e(i, k) = f[i] * mult[k] * d(i, k);
        row += e(i, k);
        dg(i, k) = -2.0 * d(i, k) * e(i, k) * y[k];
        d2(i, k) = 2.0 * d(i, k) * d(i, k) * e(i, k) *
                   (3.0 * y[k] * y[k] - x[i] * x[i]);
      }
      r[i] = row - 1.0;
    }
    s = 0.5 * r.squaredNorm() / static_cast<double>(m);
    g = VectorXd::Zero(nfree);
    for (long i = 0; i < m; ++i) {
      for (long k = 0; k < nfree; ++k) {
        g[k] += r[i] * dg(i, k);
      }
    }
    g /= static_cast<double>(m);
    h = dg.transpose() * dg;
    h += (d2.transpose() * r).asDiagonal();
    h /= static_cast<double>(m);
  };

  VectorXd y(nfree);
  for (long k = 0; k < nfree; ++k) {
    y[k] = 2.0 * kPi * static_cast<double>(k + 1);
  }
  double s = 0.0;
  VectorXd g(nfree);
  MatrixXd h(nfree, nfree);
  objective(y, s, g, h);
  bool exhausted = true;
  // The near-degenerate high-frequency group leaves a valley whose Hessian
  // eigenvalues span about twelve decades; the shifted Newton step crosses
  // it in several hundred iterations, so the cap follows the paper's 10000.
  for (int iter = 0; iter < 10000; ++iter) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(h);
    if (es.info() != Eigen::Success) {
      throw std::runtime_error("economised spring fit: Hessian eigensolve");
    }
    const VectorXd eva = es.eigenvalues();
    const MatrixXd vec = es.eigenvectors();
    const double delta = std::max(1e-16 * eva[nfree - 1], -2.0 * eva[0]);
    const VectorXd rhs = vec.transpose() * g;
    VectorXd scaled(nfree);
    for (long k = 0; k < nfree; ++k) {
      scaled[k] = -rhs[k] / (eva[k] + delta);
    }
    const VectorXd dy = vec * scaled;
    const double previous = s;
    double c = 1.0;
    bool accepted = false;
    for (int cut = 0; cut < 60; ++cut) {
      const VectorXd z = y + c * dy;
      bool ordered = z[0] >= 0.0;
      for (long k = 1; ordered && k < nfree; ++k) {
        ordered = z[k] >= z[k - 1];
      }
      if (ordered) {
        double sn = 0.0;
        VectorXd gn(nfree);
        MatrixXd hn(nfree, nfree);
        objective(z, sn, gn, hn);
        if (sn <= s) {
          y = z;
          s = sn;
          g = std::move(gn);
          h = std::move(hn);
          accepted = true;
          break;
        }
      }
      c *= 0.5;
    }
    if (!accepted) {
      exhausted = false;
      break;
    }
    if (previous - s <= 1e-12 * previous) {
      exhausted = false;
      break;
    }
  }
  if (exhausted) {
    throw std::runtime_error("economised spring fit did not converge");
  }
  return y;
}

MatrixXd factorCovariance(const MatrixXd &cov) {
  const MatrixXd sym = 0.5 * (cov + cov.transpose());
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(sym);
  if (es.info() != Eigen::Success) {
    throw std::runtime_error("GLE covariance factorisation failed");
  }
  return es.eigenvectors() *
         es.eigenvalues().cwiseMax(0.0).cwiseSqrt().asDiagonal();
}

struct GleFile {
  long nModes{0};
  long dim{0};
  std::vector<MatrixXd> drift;
  std::vector<MatrixXd> covariance;
};

GleFile readGle(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot read GLE matrices from " + path);
  }
  std::vector<double> nums;
  std::string line;
  while (std::getline(in, line)) {
    const auto hash = line.find('#');
    if (hash != std::string::npos) {
      line.resize(hash);
    }
    std::istringstream ss(line);
    double v = 0.0;
    while (ss >> v) {
      nums.push_back(v);
    }
  }
  if (nums.size() < 2) {
    throw std::runtime_error("GLE matrix file " + path + " is empty");
  }
  GleFile out;
  out.nModes = std::lround(nums[0]);
  out.dim = std::lround(nums[1]);
  if (out.nModes < 1 || out.dim < 1) {
    throw std::runtime_error("GLE matrix file " + path +
                             " has no modes or no matrix");
  }
  const long block = out.dim * out.dim;
  const long need = 2 + out.nModes * 2 * block;
  if (static_cast<long>(nums.size()) != need) {
    throw std::runtime_error("GLE matrix file " + path +
                             " does not match its header");
  }
  long cursor = 2;
  out.drift.resize(static_cast<size_t>(out.nModes));
  out.covariance.resize(static_cast<size_t>(out.nModes));
  for (long mode = 0; mode < out.nModes; ++mode) {
    MatrixXd a(out.dim, out.dim);
    MatrixXd c(out.dim, out.dim);
    for (long i = 0; i < out.dim; ++i) {
      for (long j = 0; j < out.dim; ++j) {
        a(i, j) = nums[static_cast<size_t>(cursor++)];
      }
    }
    for (long i = 0; i < out.dim; ++i) {
      for (long j = 0; j < out.dim; ++j) {
        c(i, j) = nums[static_cast<size_t>(cursor++)];
      }
    }
    out.drift[static_cast<size_t>(mode)] = std::move(a);
    out.covariance[static_cast<size_t>(mode)] = std::move(c);
  }
  return out;
}

} // namespace

void requireTrotterSprings(const std::string &springs, const char *use) {
  const std::string s = lower(springs);
  if (s == "eco" || s == "economised") {
    throw std::invalid_argument(
        std::string(use) +
        " requires Trotter springs; economised springs are refused");
  }
}

VectorXd trotterEigenvalues(long nBeads) {
  if (nBeads < 1) {
    throw std::invalid_argument("bead count must be positive");
  }
  VectorXd eva(nBeads);
  for (long k = 0; k < nBeads; ++k) {
    eva[k] = 2.0 * std::sin(kPi * static_cast<double>(k) /
                            static_cast<double>(nBeads));
  }
  return eva;
}

VectorXd ecoEigenvalues(long nBeads, double xmax) {
  if (nBeads < 1) {
    throw std::invalid_argument("bead count must be positive");
  }
  VectorXd eva = VectorXd::Zero(nBeads);
  if (nBeads == 1) {
    return eva;
  }
  if (!(xmax > 0.0)) {
    throw std::invalid_argument(
        "economised springs need a positive maximum frequency");
  }
  const VectorXd y = ecoFit(nBeads, xmax);
  for (long k = 1; k < nBeads; ++k) {
    const long pair = std::min(k, nBeads - k) - 1;
    eva[k] = y[pair] / static_cast<double>(nBeads);
  }
  return eva;
}

MatrixXd normalModeMatrix(long nBeads) {
  if (nBeads < 1) {
    throw std::invalid_argument("bead count must be positive");
  }
  MatrixXd b = MatrixXd::Zero(nBeads, nBeads);
  const double n = static_cast<double>(nBeads);
  for (long j = 0; j < nBeads; ++j) {
    b(0, j) = 1.0;
    for (long i = 1; i <= nBeads / 2; ++i) {
      b(i, j) = std::sqrt(2.0) * std::cos(2.0 * kPi * static_cast<double>(j) *
                                          static_cast<double>(i) / n);
    }
    for (long i = nBeads / 2 + 1; i < nBeads; ++i) {
      b(i, j) = std::sqrt(2.0) * std::sin(2.0 * kPi * static_cast<double>(j) *
                                          static_cast<double>(i) / n);
    }
  }
  if (nBeads % 2 == 0) {
    const long mid = nBeads / 2;
    for (long j = 0; j < nBeads; ++j) {
      b(mid, j) = (j % 2 == 0) ? 1.0 : -1.0;
    }
  }
  b /= std::sqrt(n);
  return b;
}

RingPolymer::RingPolymer(long nAtoms, std::vector<double> masses,
                         std::vector<int> atomicNumbers, std::vector<char> free,
                         Options opt)
    : opt_(std::move(opt)),
      nAtoms_(nAtoms),
      nDof_(3 * nAtoms),
      nBeads_(opt_.beads),
      atomicNumbers_(std::move(atomicNumbers)),
      free_(std::move(free)),
      rng_(opt_.seed == 0 ? 1 : opt_.seed) {
  if (nAtoms_ < 1 || nBeads_ < 1) {
    throw std::invalid_argument("path integral needs atoms and beads");
  }
  if (static_cast<long>(masses.size()) != nAtoms_ ||
      static_cast<long>(atomicNumbers_.size()) != nAtoms_ ||
      static_cast<long>(free_.size()) != nDof_) {
    throw std::invalid_argument("path integral mass, number or mask size");
  }
  if (!(opt_.temperature > 0.0) || !(opt_.kB > 0.0) || !(opt_.hbar > 0.0) ||
      !(opt_.dt > 0.0) || !(opt_.pileTau > 0.0) || !(opt_.pileScale > 0.0)) {
    throw std::invalid_argument(
        "path integral temperature, timestep and damping must be positive");
  }
  if (opt_.springs == Springs::Eco && opt_.thermostat == Thermostat::Piglet) {
    throw std::invalid_argument(
        "economised springs cannot be combined with a normal-mode GLE");
  }
  mass_.assign(static_cast<size_t>(nDof_), 0.0);
  nFree_ = 0;
  for (long i = 0; i < nAtoms_; ++i) {
    if (!(masses[static_cast<size_t>(i)] > 0.0)) {
      throw std::invalid_argument("path integral atom has no mass");
    }
    for (int axis = 0; axis < 3; ++axis) {
      const long a = 3 * i + axis;
      mass_[static_cast<size_t>(a)] = masses[static_cast<size_t>(i)];
      if (free_[static_cast<size_t>(a)]) {
        freeIndex_.push_back(a);
        ++nFree_;
      }
    }
  }
  if (nFree_ < 1) {
    throw std::invalid_argument("path integral has no free coordinate");
  }
  modes_ = normalModeMatrix(nBeads_);
  const double omegan =
      static_cast<double>(nBeads_) * opt_.kB * opt_.temperature / opt_.hbar;
  if (opt_.springs == Springs::Eco) {
    const double xmax =
        opt_.ecoOmegaMax * opt_.hbar / (opt_.kB * opt_.temperature);
    omegaK_ = omegan * ecoEigenvalues(nBeads_, xmax);
  } else {
    omegaK_ = omegan * trotterEigenvalues(nBeads_);
  }
  q_.assign(static_cast<size_t>(nBeads_), VectorXd::Zero(nDof_));
  p_.assign(static_cast<size_t>(nBeads_), VectorXd::Zero(nDof_));
  f_.assign(static_cast<size_t>(nBeads_), VectorXd::Zero(nDof_));
  qnm_.assign(static_cast<size_t>(nBeads_), VectorXd::Zero(nDof_));
  pnm_.assign(static_cast<size_t>(nBeads_), VectorXd::Zero(nDof_));
  if (opt_.thermostat == Thermostat::Piglet) {
    initGle();
  }
  thermalMomenta();
}

void RingPolymer::initGle() {
  if (nBeads_ == 1) {
    return;
  }
  if (opt_.gleFile.empty()) {
    throw std::invalid_argument("normal-mode GLE needs a matrix file");
  }
  const GleFile file = readGle(opt_.gleFile);
  long first = 0;
  if (file.nModes == nBeads_) {
    first = 1;
  } else if (file.nModes != nBeads_ - 1) {
    throw std::invalid_argument(
        "GLE matrix count must be the bead count or one less");
  }
  const double h = 0.5 * opt_.dt;
  gle_.resize(static_cast<size_t>(nBeads_ - 1));
  for (long k = 1; k < nBeads_; ++k) {
    const long src = first + (k - 1);
    ModeGle mode;
    mode.drift = file.drift[static_cast<size_t>(src)];
    mode.covariance = file.covariance[static_cast<size_t>(src)];
    mode.propagate = (-mode.drift * h).exp();
    const MatrixXd cov =
        opt_.kB * (mode.covariance - mode.propagate * mode.covariance *
                                         mode.propagate.transpose());
    mode.noise = factorCovariance(cov);
    mode.extended = MatrixXd::Zero(file.dim, nFree_);
    gle_[static_cast<size_t>(k - 1)] = std::move(mode);
  }
}

void RingPolymer::setAllBeads(const double *q) {
  if (q == nullptr) {
    throw std::invalid_argument("path integral positions are missing");
  }
  for (long bead = 0; bead < nBeads_; ++bead) {
    for (long a = 0; a < nDof_; ++a) {
      q_[static_cast<size_t>(bead)][a] = q[a];
    }
  }
  if (constrain_) {
    projectPosition();
  }
  haveForces_ = false;
}

void RingPolymer::setBeads(const std::vector<VectorXd> &beads) {
  if (static_cast<long>(beads.size()) != nBeads_) {
    throw std::invalid_argument("path integral bead count mismatch");
  }
  for (long bead = 0; bead < nBeads_; ++bead) {
    if (beads[static_cast<size_t>(bead)].size() != nDof_) {
      throw std::invalid_argument("path integral bead has the wrong length");
    }
    q_[static_cast<size_t>(bead)] = beads[static_cast<size_t>(bead)];
  }
  if (constrain_) {
    projectPosition();
  }
  haveForces_ = false;
}

void RingPolymer::setMomenta(const std::vector<VectorXd> &momenta) {
  if (static_cast<long>(momenta.size()) != nBeads_) {
    throw std::invalid_argument("path integral bead count mismatch");
  }
  for (long bead = 0; bead < nBeads_; ++bead) {
    if (momenta[static_cast<size_t>(bead)].size() != nDof_) {
      throw std::invalid_argument(
          "path integral momentum has the wrong length");
    }
    for (long a = 0; a < nDof_; ++a) {
      p_[static_cast<size_t>(bead)][a] =
          free_[static_cast<size_t>(a)] ? momenta[static_cast<size_t>(bead)][a]
                                        : 0.0;
    }
  }
}

void RingPolymer::setHyperplane(const VectorXd &normal,
                                const VectorXd &origin) {
  if (normal.size() != nDof_ || origin.size() != nDof_) {
    throw std::invalid_argument("hyperplane vectors have the wrong length");
  }
  planeNormal_ = VectorXd::Zero(nDof_);
  planeOrigin_ = origin;
  double norm = 0.0;
  for (long a : freeIndex_) {
    planeNormal_[a] = normal[a];
    norm += normal[a] * normal[a];
  }
  if (!(norm > 0.0)) {
    throw std::invalid_argument("hyperplane normal has no free component");
  }
  planeNormal_ /= std::sqrt(norm);
  constrain_ = true;
  projectPosition();
  projectMomentum();
}

double RingPolymer::gauss() {
  rng_ = rng_ * 6364136223846793005ULL + 1ULL;
  const double u1 =
      std::max((rng_ >> 11) * (1.0 / 9007199254740992.0), 1.0e-16);
  rng_ = rng_ * 6364136223846793005ULL + 1ULL;
  const double u2 = (rng_ >> 11) * (1.0 / 9007199254740992.0);
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * u2);
}

void RingPolymer::thermalMomenta() {
  for (long bead = 0; bead < nBeads_; ++bead) {
    p_[static_cast<size_t>(bead)].setZero();
  }
  toNormal(p_, pnm_);
  const double tSim = static_cast<double>(nBeads_) * opt_.kB * opt_.temperature;
  for (long k = 0; k < nBeads_; ++k) {
    const bool gleMode = opt_.thermostat == Thermostat::Piglet && k > 0;
    if (gleMode) {
      continue;
    }
    for (long a : freeIndex_) {
      pnm_[static_cast<size_t>(k)][a] =
          gauss() * std::sqrt(mass_[static_cast<size_t>(a)] * tSim);
    }
  }
  if (opt_.thermostat == Thermostat::Piglet) {
    for (long k = 1; k < nBeads_; ++k) {
      ModeGle &mode = gle_[static_cast<size_t>(k - 1)];
      const MatrixXd factor = factorCovariance(opt_.kB * mode.covariance);
      MatrixXd noise(mode.extended.rows(), mode.extended.cols());
      for (long r = 0; r < noise.rows(); ++r) {
        for (long c = 0; c < noise.cols(); ++c) {
          noise(r, c) = gauss();
        }
      }
      mode.extended = factor * noise;
      for (long col = 0; col < nFree_; ++col) {
        const long a = freeIndex_[static_cast<size_t>(col)];
        pnm_[static_cast<size_t>(k)][a] =
            mode.extended(0, col) * std::sqrt(mass_[static_cast<size_t>(a)]);
      }
    }
  }
  fromNormal(pnm_, p_);
  haveForces_ = false;
}

void RingPolymer::toNormal(const std::vector<VectorXd> &src,
                           std::vector<VectorXd> &dst) const {
  MatrixXd packed(nDof_, nBeads_);
  for (long j = 0; j < nBeads_; ++j) {
    packed.col(j) = src[static_cast<size_t>(j)];
  }
  const MatrixXd out = packed * modes_.transpose();
  for (long k = 0; k < nBeads_; ++k) {
    dst[static_cast<size_t>(k)] = out.col(k);
  }
}

void RingPolymer::fromNormal(const std::vector<VectorXd> &src,
                             std::vector<VectorXd> &dst) const {
  MatrixXd packed(nDof_, nBeads_);
  for (long k = 0; k < nBeads_; ++k) {
    packed.col(k) = src[static_cast<size_t>(k)];
  }
  const MatrixXd out = packed * modes_;
  for (long j = 0; j < nBeads_; ++j) {
    dst[static_cast<size_t>(j)] = out.col(j);
  }
}

void RingPolymer::forces(Potential &pot, const double *box) {
  double zeroBox[9] = {};
  const double *cell = box != nullptr ? box : zeroBox;
  std::vector<const double *> pos(static_cast<size_t>(nBeads_));
  std::vector<const int *> nrs(static_cast<size_t>(nBeads_));
  std::vector<double *> frc(static_cast<size_t>(nBeads_));
  std::vector<const double *> boxes(static_cast<size_t>(nBeads_), cell);
  for (long bead = 0; bead < nBeads_; ++bead) {
    pos[static_cast<size_t>(bead)] = q_[static_cast<size_t>(bead)].data();
    nrs[static_cast<size_t>(bead)] = atomicNumbers_.data();
    frc[static_cast<size_t>(bead)] = f_[static_cast<size_t>(bead)].data();
  }
  std::vector<double> energies(static_cast<size_t>(nBeads_), 0.0);
  std::vector<double> variances(static_cast<size_t>(nBeads_), 0.0);
  pot.forceBatch(nBeads_, nAtoms_, pos.data(), nrs.data(), frc.data(),
                 energies.data(), variances.data(), boxes.data());
  ++batches_;
  for (long bead = 0; bead < nBeads_; ++bead) {
    for (long a = 0; a < nDof_; ++a) {
      if (!free_[static_cast<size_t>(a)]) {
        f_[static_cast<size_t>(bead)][a] = 0.0;
      }
    }
  }
  haveForces_ = true;
}

VectorXd RingPolymer::centroid() const {
  VectorXd c = VectorXd::Zero(nDof_);
  for (long bead = 0; bead < nBeads_; ++bead) {
    c += q_[static_cast<size_t>(bead)];
  }
  c /= static_cast<double>(nBeads_);
  return c;
}

VectorXd RingPolymer::centroidVelocity() const {
  std::vector<VectorXd> pnm(static_cast<size_t>(nBeads_),
                            VectorXd::Zero(nDof_));
  toNormal(p_, pnm);
  VectorXd v = VectorXd::Zero(nDof_);
  const double scale = std::sqrt(static_cast<double>(nBeads_));
  for (long a = 0; a < nDof_; ++a) {
    const double m = mass_[static_cast<size_t>(a)];
    if (m > 0.0) {
      v[a] = pnm[0][a] / (scale * m);
    }
  }
  return v;
}

double RingPolymer::kineticCv() const {
  double k = 0.5 * static_cast<double>(nFree_) * opt_.kB * opt_.temperature;
  if (!haveForces_) {
    return k;
  }
  const VectorXd c = centroid();
  double virial = 0.0;
  for (long bead = 0; bead < nBeads_; ++bead) {
    for (long a : freeIndex_) {
      virial += (q_[static_cast<size_t>(bead)][a] - c[a]) *
                f_[static_cast<size_t>(bead)][a];
    }
  }
  k += -0.5 / static_cast<double>(nBeads_) * virial;
  return k;
}

double RingPolymer::meanForce() const {
  if (!constrain_ || recorded_ < 1) {
    return 0.0;
  }
  return forceSum_ / static_cast<double>(recorded_);
}

void RingPolymer::resetAverages() {
  recorded_ = 0;
  kineticSum_ = 0.0;
  forceSum_ = 0.0;
}

void RingPolymer::thermostat(double h) {
  toNormal(p_, pnm_);
  const double tSim = static_cast<double>(nBeads_) * opt_.kB * opt_.temperature;
  auto langevin = [&](long k, double tau) {
    const double damp = std::exp(-h / tau);
    const double noise = std::sqrt(tSim * (1.0 - damp * damp));
    for (long a : freeIndex_) {
      const double sm = std::sqrt(mass_[static_cast<size_t>(a)]);
      double pms = pnm_[static_cast<size_t>(k)][a] / sm;
      pms = damp * pms + noise * gauss();
      pnm_[static_cast<size_t>(k)][a] = pms * sm;
    }
  };
  langevin(0, opt_.pileTau);
  for (long k = 1; k < nBeads_; ++k) {
    if (opt_.thermostat == Thermostat::Piglet) {
      ModeGle &mode = gle_[static_cast<size_t>(k - 1)];
      for (long col = 0; col < nFree_; ++col) {
        const long a = freeIndex_[static_cast<size_t>(col)];
        mode.extended(0, col) = pnm_[static_cast<size_t>(k)][a] /
                                std::sqrt(mass_[static_cast<size_t>(a)]);
      }
      MatrixXd noise(mode.extended.rows(), mode.extended.cols());
      for (long r = 0; r < noise.rows(); ++r) {
        for (long c = 0; c < noise.cols(); ++c) {
          noise(r, c) = gauss();
        }
      }
      mode.extended = mode.propagate * mode.extended + mode.noise * noise;
      for (long col = 0; col < nFree_; ++col) {
        const long a = freeIndex_[static_cast<size_t>(col)];
        pnm_[static_cast<size_t>(k)][a] =
            mode.extended(0, col) * std::sqrt(mass_[static_cast<size_t>(a)]);
      }
    } else {
      const double tau = 1.0 / (2.0 * opt_.pileScale * omegaK_[k]);
      langevin(k, tau);
    }
  }
  fromNormal(pnm_, p_);
}

void RingPolymer::kick(double h, bool dropParallel) {
  VectorXd removal = VectorXd::Zero(nDof_);
  if (dropParallel && constrain_) {
    VectorXd fc = VectorXd::Zero(nDof_);
    for (long bead = 0; bead < nBeads_; ++bead) {
      fc += f_[static_cast<size_t>(bead)];
    }
    fc /= static_cast<double>(nBeads_);
    double along = 0.0;
    for (long a : freeIndex_) {
      along += planeNormal_[a] * fc[a];
    }
    removal = along * planeNormal_;
  }
  for (long bead = 0; bead < nBeads_; ++bead) {
    for (long a : freeIndex_) {
      p_[static_cast<size_t>(bead)][a] +=
          (f_[static_cast<size_t>(bead)][a] - removal[a]) * h;
    }
  }
}

void RingPolymer::propagate(double h) {
  toNormal(q_, qnm_);
  toNormal(p_, pnm_);
  for (long a : freeIndex_) {
    const double m = mass_[static_cast<size_t>(a)];
    qnm_[0][a] += pnm_[0][a] / m * h;
    for (long k = 1; k < nBeads_; ++k) {
      const double omega = omegaK_[k];
      const double c = std::cos(omega * h);
      const double s = std::sin(omega * h);
      const double pk = pnm_[static_cast<size_t>(k)][a];
      const double qk = qnm_[static_cast<size_t>(k)][a];
      pnm_[static_cast<size_t>(k)][a] = c * pk - m * omega * s * qk;
      qnm_[static_cast<size_t>(k)][a] = c * qk + s * pk / (omega * m);
    }
  }
  fromNormal(qnm_, q_);
  fromNormal(pnm_, p_);
}

// The centroid mode is the bead sum over sqrt(P), so a change d of that
// mode moves every bead by d / sqrt(P) and leaves the other modes alone.
void RingPolymer::projectPosition() {
  if (!constrain_) {
    return;
  }
  const VectorXd c = centroid();
  double sigma = 0.0;
  for (long a : freeIndex_) {
    sigma += planeNormal_[a] * (c[a] - planeOrigin_[a]);
  }
  for (long bead = 0; bead < nBeads_; ++bead) {
    for (long a : freeIndex_) {
      q_[static_cast<size_t>(bead)][a] -= sigma * planeNormal_[a];
    }
  }
}

void RingPolymer::projectMomentum() {
  if (!constrain_) {
    return;
  }
  VectorXd sum = VectorXd::Zero(nDof_);
  for (long bead = 0; bead < nBeads_; ++bead) {
    sum += p_[static_cast<size_t>(bead)];
  }
  // p_0 = sum / sqrt(P); remove lambda n from p_0.
  const double scale = std::sqrt(static_cast<double>(nBeads_));
  double num = 0.0;
  double den = 0.0;
  for (long a : freeIndex_) {
    const double m = mass_[static_cast<size_t>(a)];
    num += planeNormal_[a] * (sum[a] / scale) / m;
    den += planeNormal_[a] * planeNormal_[a] / m;
  }
  if (den > 0.0) {
    const double shift = num / den / scale;
    for (long bead = 0; bead < nBeads_; ++bead) {
      for (long a : freeIndex_) {
        p_[static_cast<size_t>(bead)][a] -= shift * planeNormal_[a];
      }
    }
  }
}

void RingPolymer::step(Potential &pot, const double *box, bool record) {
  const double half = 0.5 * opt_.dt;
  thermostat(half);
  projectMomentum();
  forces(pot, box);
  kick(half, true);
  projectMomentum();
  propagate(half);
  propagate(half);
  projectPosition();
  forces(pot, box);
  if (record) {
    kineticSum_ += kineticCv();
    if (constrain_) {
      VectorXd fc = VectorXd::Zero(nDof_);
      for (long bead = 0; bead < nBeads_; ++bead) {
        fc += f_[static_cast<size_t>(bead)];
      }
      fc /= static_cast<double>(nBeads_);
      double along = 0.0;
      for (long a : freeIndex_) {
        along += planeNormal_[a] * fc[a];
      }
      forceSum_ += along;
    }
    ++recorded_;
  }
  kick(half, true);
  projectMomentum();
  thermostat(half);
  projectMomentum();
}

void RingPolymer::nveStep(Potential &pot, const double *box) {
  if (constrain_) {
    throw std::logic_error("path integral NVE step with a hyperplane set");
  }
  const double half = 0.5 * opt_.dt;
  if (!haveForces_) {
    forces(pot, box);
  }
  kick(half, false);
  propagate(opt_.dt);
  forces(pot, box);
  kick(half, false);
}

Sample RingPolymer::sample(Potential &pot, const double *box,
                           long equilibration, long production) {
  if (equilibration < 0 || production < 1) {
    throw std::invalid_argument("path integral sample length");
  }
  for (long step = 0; step < equilibration; ++step) {
    this->step(pot, box, false);
  }
  resetAverages();
  const long batchesBefore = batches_;
  for (long step = 0; step < production; ++step) {
    this->step(pot, box, true);
  }
  Sample out;
  out.kineticCv = kineticSum_ / static_cast<double>(recorded_);
  out.meanForce = meanForce();
  out.batches = batches_ - batchesBefore;
  return out;
}

} // namespace eonc::pathintegral

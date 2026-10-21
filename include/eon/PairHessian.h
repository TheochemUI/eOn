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
#pragma once

#include "BaseStructures.h"
#include "Eigen.h"
#include "ObjectiveFunction.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace eonc {
namespace pairhess {

inline void addPairBlock(Eigen::MatrixXd &P, int i, int j,
                         const Eigen::Vector3d &rij, double k_par,
                         double k_perp) {
  const double r = rij.norm();
  if (r < 1.0e-14 || (!std::isfinite(k_par) && !std::isfinite(k_perp))) {
    return;
  }
  if (k_par == 0.0 && k_perp == 0.0) {
    return;
  }
  const Eigen::Vector3d u = rij / r;
  const Eigen::Matrix3d H = k_perp * Eigen::Matrix3d::Identity() +
                            (k_par - k_perp) * (u * u.transpose());
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      const int ia = 3 * i + a;
      const int ib = 3 * i + b;
      const int ja = 3 * j + a;
      const int jb = 3 * j + b;
      P(ia, ib) += H(a, b);
      P(ja, jb) += H(a, b);
      P(ia, jb) -= H(a, b);
      P(ja, ib) -= H(a, b);
    }
  }
}

inline void ljDerivs(double r, double eps, double sigma, double &vp,
                     double &vpp) {
  const double inv = sigma / r;
  const double x6 = inv * inv * inv * inv * inv * inv;
  const double x12 = x6 * x6;
  vp = 4.0 * eps * (-12.0 * x12 + 6.0 * x6) / r;
  vpp = 4.0 * eps * (156.0 * x12 - 42.0 * x6) / (r * r);
}

inline void morseDerivs(double r, double D, double a, double re, double &vp,
                        double &vpp) {
  const double e1 = std::exp(-a * (r - re));
  const double e2 = e1 * e1;
  vp = 2.0 * D * a * (e1 - e2);
  vpp = 2.0 * D * a * a * (2.0 * e2 - e1);
}

inline bool isModelHess(const std::string &kind) {
  return kind == "lindh_full" || kind == "swart" || kind == "fischer" ||
         kind == "schlegel";
}

/// Rank-1 update \f$k\,(\nabla q)(\nabla q)^\top\f$ for a scalar internal.
inline void addScalarInternal(Eigen::MatrixXd &P, const int *atoms,
                              const Eigen::Vector3d *g, int n_atoms, double k) {
  if (!(k > 0.0) || !std::isfinite(k)) {
    return;
  }
  for (int a = 0; a < n_atoms; ++a) {
    for (int b = 0; b < n_atoms; ++b) {
      for (int x = 0; x < 3; ++x) {
        for (int y = 0; y < 3; ++y) {
          P(3 * atoms[a] + x, 3 * atoms[b] + y) += k * g[a][x] * g[b][y];
        }
      }
    }
  }
}

/// Wilson \(B\) rows for the valence angle \(i{-}j{-}k\) (\(j\) central).
inline bool angleGrads(const Eigen::Vector3d &ri, const Eigen::Vector3d &rj,
                       const Eigen::Vector3d &rk, Eigen::Vector3d g[3]) {
  Eigen::Vector3d u = ri - rj;
  Eigen::Vector3d v = rk - rj;
  const double pu = u.norm();
  const double pv = v.norm();
  if (pu < 1.0e-14 || pv < 1.0e-14) {
    return false;
  }
  u /= pu;
  v /= pv;
  double c = u.dot(v);
  c = std::max(-1.0 + 1.0e-14, std::min(1.0 - 1.0e-14, c));
  const double s = std::sqrt(std::max(0.0, 1.0 - c * c));
  if (s < 1.0e-8) {
    return false;
  }
  g[0] = (c * u - v) / (pu * s);
  g[2] = (c * v - u) / (pv * s);
  g[1] = -g[0] - g[2];
  return true;
}

/// Wilson \(B\) rows for the torsion \(i{-}j{-}k{-}l\).
inline bool torsionGrads(const Eigen::Vector3d &ri, const Eigen::Vector3d &rj,
                         const Eigen::Vector3d &rk, const Eigen::Vector3d &rl,
                         Eigen::Vector3d g[4]) {
  const Eigen::Vector3d vij = ri - rj;
  const Eigen::Vector3d vkj = rk - rj;
  const Eigen::Vector3d vkl = rk - rl;
  Eigen::Vector3d n1 = vij.cross(vkj);
  Eigen::Vector3d n2 = vkl.cross(vkj);
  const double n1n = n1.norm();
  const double n2n = n2.norm();
  const double b2 = vkj.norm();
  if (n1n < 1.0e-14 || n2n < 1.0e-14 || b2 < 1.0e-14) {
    return false;
  }
  g[0] = -(b2 / (n1n * n1n)) * n1;
  g[3] = (b2 / (n2n * n2n)) * n2;
  const double c1 = vij.dot(vkj) / (b2 * b2);
  const double c2 = vkl.dot(vkj) / (b2 * b2);
  g[1] = -g[0] + c1 * g[0] - c2 * g[3];
  g[2] = -g[3] - c1 * g[0] + c2 * g[3];
  return true;
}

inline double lindhRho(double r, double r_nn, double alpha) {
  return std::exp(alpha * (r_nn * r_nn - r * r));
}

inline double swartScreen(double r, double r_nn) {
  return 1.0 / (1.0 + std::exp(12.0 * (r / r_nn - 1.2)));
}

inline void addModelHess(Eigen::MatrixXd &P, const Eigen::VectorXd &pos,
                         const std::string &kind, double A, double mu,
                         double r_nn, double rcut, ObjectiveFunction &objf) {
  const int nat = static_cast<int>(pos.size()) / 3;
  const double alpha = A > 0.0 ? A : 1.0;
  const double scale = mu > 0.0 ? mu : 1.0;
  const double r_bond = std::min(rcut, 1.35 * r_nn);
  const double r_cov = 0.5 * r_nn;
  auto mic = [&](int i, int j) {
    Eigen::Vector3d dr = pos.segment<3>(3 * i) - pos.segment<3>(3 * j);
    objf.minimumImage(dr);
    return dr;
  };
  std::vector<std::vector<int>> nb(static_cast<size_t>(nat));
  const double r_pair = (kind == "swart") ? rcut : r_bond;
  for (int i = 0; i < nat; ++i) {
    for (int j = i + 1; j < nat; ++j) {
      const Eigen::Vector3d rij = mic(i, j);
      const double r = rij.norm();
      if (r >= r_pair || r < 1.0e-14) {
        continue;
      }
      double k = 0.0;
      if (kind == "schlegel") {
        const double den = std::max(r - 0.45 * r_nn, 0.15);
        k = scale * 1.734 / (den * den * den);
      } else if (kind == "fischer") {
        k = scale * 0.3601 * std::exp(-1.944 * (r - r_cov));
      } else {
        k = scale * lindhRho(r, r_nn, alpha);
        if (kind == "swart") {
          k *= swartScreen(r, r_nn);
        }
      }
      addPairBlock(P, i, j, rij, k, 0.0);
      if (r < r_bond) {
        nb[static_cast<size_t>(i)].push_back(j);
        nb[static_cast<size_t>(j)].push_back(i);
      }
    }
  }
  for (int j = 0; j < nat; ++j) {
    const auto &nbr = nb[static_cast<size_t>(j)];
    for (size_t a = 0; a < nbr.size(); ++a) {
      for (size_t b = a + 1; b < nbr.size(); ++b) {
        const int i = nbr[a];
        const int k = nbr[b];
        const Eigen::Vector3d ri = pos.segment<3>(3 * i);
        const Eigen::Vector3d rj = pos.segment<3>(3 * j);
        const Eigen::Vector3d rk = pos.segment<3>(3 * k);
        Eigen::Vector3d g[3];
        Eigen::Vector3d rji = ri - rj;
        Eigen::Vector3d rjk = rk - rj;
        objf.minimumImage(rji);
        objf.minimumImage(rjk);
        if (!angleGrads(rj + rji, rj, rj + rjk, g)) {
          continue;
        }
        const double rij = rji.norm();
        const double rkj = rjk.norm();
        double kang = 0.16 * scale;
        if (kind == "fischer") {
          kang =
              scale *
              (0.089 + 0.11 * std::pow(r_cov * r_cov, -0.42) *
                           std::exp(-0.44 * ((rij - r_cov) + (rkj - r_cov))));
        } else if (kind == "lindh_full" || kind == "swart") {
          kang = 0.15 * scale * lindhRho(rij, r_nn, alpha) *
                 lindhRho(rkj, r_nn, alpha);
          if (kind == "swart") {
            kang *= swartScreen(rij, r_nn) * swartScreen(rkj, r_nn);
          }
        }
        const int atoms[3] = {i, j, k};
        addScalarInternal(P, atoms, g, 3, kang);
      }
    }
  }
  for (int j = 0; j < nat; ++j) {
    const auto &nj = nb[static_cast<size_t>(j)];
    for (int k : nj) {
      if (k <= j) {
        continue;
      }
      const auto &nk = nb[static_cast<size_t>(k)];
      for (int i : nj) {
        if (i == k) {
          continue;
        }
        for (int l : nk) {
          if (l == j || l == i) {
            continue;
          }
          Eigen::Vector3d rijv = pos.segment<3>(3 * i) - pos.segment<3>(3 * j);
          Eigen::Vector3d rkjv = pos.segment<3>(3 * k) - pos.segment<3>(3 * j);
          Eigen::Vector3d rklv = pos.segment<3>(3 * k) - pos.segment<3>(3 * l);
          objf.minimumImage(rijv);
          objf.minimumImage(rkjv);
          objf.minimumImage(rklv);
          const Eigen::Vector3d rj = pos.segment<3>(3 * j);
          Eigen::Vector3d g[4];
          if (!torsionGrads(rj + rijv, rj, rj + rkjv, rj + rkjv - rklv, g)) {
            continue;
          }
          const double rij = rijv.norm();
          const double rjk = rkjv.norm();
          const double rkl = rklv.norm();
          double ktor = 0.01 * scale;
          if (kind == "fischer") {
            ktor = 0.0015 * scale;
          } else if (kind == "lindh_full" || kind == "swart") {
            ktor = 0.005 * scale * lindhRho(rij, r_nn, alpha) *
                   lindhRho(rjk, r_nn, alpha) * lindhRho(rkl, r_nn, alpha);
            if (kind == "swart") {
              ktor *= swartScreen(rij, r_nn) * swartScreen(rjk, r_nn) *
                      swartScreen(rkl, r_nn);
            }
          }
          const int atoms[4] = {i, j, k, l};
          addScalarInternal(P, atoms, g, 4, ktor);
        }
      }
    }
  }
}

/// Analytic pair or model Hessian. `kind` is pair / pair_abs / pair_full /
/// lindh / lindh_full / exp / c1 / fischer / schlegel / swart.
inline Eigen::MatrixXd build(const Eigen::VectorXd &pos,
                             const std::string &kind, PotType pot, double A,
                             double mu, double rcut_in,
                             ObjectiveFunction &objf) {
  const int n = static_cast<int>(pos.size());
  const int nat = n / 3;
  std::vector<double> nn(static_cast<size_t>(nat),
                         std::numeric_limits<double>::infinity());
  auto mic = [&](int i, int j) {
    Eigen::Vector3d dr = pos.segment<3>(3 * i) - pos.segment<3>(3 * j);
    objf.minimumImage(dr);
    return dr;
  };
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < nat; ++j) {
      if (i == j) {
        continue;
      }
      const double r = mic(i, j).norm();
      if (r > 1.0e-12 && r < nn[static_cast<size_t>(i)]) {
        nn[static_cast<size_t>(i)] = r;
      }
    }
  }
  double r_nn = 0.0;
  for (double v : nn) {
    if (std::isfinite(v)) {
      r_nn = std::max(r_nn, v);
    }
  }
  if (r_nn <= 0.0) {
    r_nn = 1.0;
  }
  const double rcut = rcut_in > 0.0 ? rcut_in : 2.5 * r_nn;
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n, n);
  if (isModelHess(kind)) {
    addModelHess(P, pos, kind, A, mu, r_nn, rcut, objf);
    const double shift = 1.0e-6 * std::max(1.0, P.diagonal().mean());
    P.diagonal().array() += shift;
    return P;
  }
  for (int i = 0; i < nat; ++i) {
    for (int j = i + 1; j < nat; ++j) {
      const Eigen::Vector3d rij = mic(i, j);
      const double r = rij.norm();
      if (r >= rcut || r < 1.0e-14) {
        continue;
      }
      if (kind == "pair" || kind == "pair_abs" || kind == "pair_full") {
        double vp = 0.0;
        double vpp = 0.0;
        if (pot == PotType::MORSE_PT) {
          morseDerivs(r, 0.7102, 1.6047, 2.897, vp, vpp);
        } else {
          ljDerivs(r, 1.0, 1.0, vp, vpp);
        }
        double k_par = std::max(vpp, 0.0);
        double k_perp = std::max(vp / r, 0.0);
        if (kind == "pair_abs") {
          k_par = std::abs(vpp);
          k_perp = std::abs(vp / r);
        } else if (kind == "pair_full") {
          k_par = vpp;
          k_perp = vp / r;
        }
        addPairBlock(P, i, j, rij, k_par, k_perp);
      } else if (kind == "lindh") {
        const double alpha = A > 0.0 ? A : 1.0;
        const double k = mu * std::exp(alpha * (r_nn * r_nn - r * r));
        addPairBlock(P, i, j, rij, k, 0.0);
      } else if (kind == "c1") {
        const double x = r / rcut;
        const double w = mu * (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
        addPairBlock(P, i, j, rij, w, w);
      } else {
        const double w = mu * std::exp(-A * (r / r_nn - 1.0));
        addPairBlock(P, i, j, rij, w, w);
      }
    }
  }
  const double shift = 1.0e-6 * std::max(1.0, P.diagonal().mean());
  P.diagonal().array() += shift;
  return P;
}

/// Pair potential whose analytic force matches an objective force.
/// Lennard-Jones uses epsilon = sigma = 1. Morse uses D0 0.7102, alpha
/// 1.6047, r0 2.897. Both are the curvatures `build` already evaluates.
enum class KnownPair { None, LennardJones, Morse };

inline PotType knownPot(KnownPair known) {
  return known == KnownPair::Morse ? PotType::MORSE_PT : PotType::LJ;
}

inline Eigen::VectorXd withoutNetForce(Eigen::VectorXd force) {
  const int nat = static_cast<int>(force.size() / 3);
  if (nat <= 1) {
    return force;
  }
  Eigen::Vector3d mean = Eigen::Vector3d::Zero();
  for (int i = 0; i < nat; ++i) {
    mean += force.segment<3>(3L * i);
  }
  mean /= static_cast<double>(nat);
  for (int i = 0; i < nat; ++i) {
    force.segment<3>(3L * i) -= mean;
  }
  return force;
}

/// Analytic pair force, f = -dV/dx, with the mean force removed.
inline Eigen::VectorXd analyticPairForce(const Eigen::VectorXd &pos,
                                         ObjectiveFunction &objf, bool morse,
                                         double cutoff) {
  const int nat = static_cast<int>(pos.size() / 3);
  Eigen::VectorXd force = Eigen::VectorXd::Zero(pos.size());
  for (int i = 0; i < nat; ++i) {
    for (int j = i + 1; j < nat; ++j) {
      Eigen::Vector3d dr = pos.segment<3>(3L * i) - pos.segment<3>(3L * j);
      objf.minimumImage(dr);
      const double r = dr.norm();
      if (!(r > 1.0e-14 && r < cutoff)) {
        continue;
      }
      double vp = 0.0;
      double vpp = 0.0;
      if (morse) {
        morseDerivs(r, 0.7102, 1.6047, 2.897, vp, vpp);
      } else {
        ljDerivs(r, 1.0, 1.0, vp, vpp);
      }
      if (!std::isfinite(vp)) {
        continue;
      }
      const Eigen::Vector3d u = dr / r;
      force.segment<3>(3L * i) -= vp * u;
      force.segment<3>(3L * j) += vp * u;
    }
  }
  return withoutNetForce(std::move(force));
}

struct PairMatch {
  KnownPair kind{KnownPair::None};
  /// Cutoff of the candidate whose force matched. Zero when nothing matched.
  double cutoff{0.0};
};

/// Lennard-Jones cutoffs 15 and 40, and Morse cutoffs 7 and 9.5.
/// A match needs a force above 1 and a relative residual below 1e-5,
/// so a near-zero force is not claimed as every pair model at once.
inline PairMatch matchKnownPair(const Eigen::VectorXd &pos,
                                const Eigen::VectorXd &force,
                                ObjectiveFunction &objf) {
  PairMatch best;
  if (pos.size() < 6 || pos.size() % 3 != 0 || force.size() != pos.size()) {
    return best;
  }
  const Eigen::VectorXd observed = withoutNetForce(force);
  const double scale = observed.norm();
  if (!(scale > 1.0) || !observed.allFinite()) {
    return best;
  }
  struct Cand {
    KnownPair kind;
    bool morse;
    double cutoff;
  };
  const Cand cands[] = {
      {KnownPair::LennardJones, false, 15.0},
      {KnownPair::LennardJones, false, 40.0},
      {KnownPair::Morse, true, 7.0},
      {KnownPair::Morse, true, 9.5},
  };
  double bestRel = 1.0e-5;
  for (const Cand &cand : cands) {
    const Eigen::VectorXd model =
        analyticPairForce(pos, objf, cand.morse, cand.cutoff);
    if (!model.allFinite()) {
      continue;
    }
    const double rel = (model - observed).norm() / scale;
    if (rel < bestRel) {
      bestRel = rel;
      best.kind = cand.kind;
      best.cutoff = cand.cutoff;
    }
  }
  return best;
}

/// P^{-1} f for the positive pair Hessian. Returns `force` if the
/// factorisation is not finite.
inline Eigen::VectorXd pairNewtonForce(const Eigen::VectorXd &pos,
                                       const Eigen::VectorXd &force,
                                       KnownPair known,
                                       ObjectiveFunction &objf,
                                       double cutoff) {
  if (known == KnownPair::None || pos.size() < 6 || pos.size() % 3 != 0) {
    return force;
  }
  const Eigen::MatrixXd P =
      build(pos, "pair", knownPot(known), 3.0, 1.0, cutoff, objf);
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(P);
  const Eigen::VectorXd step = lu.solve(force);
  if (!step.allFinite()) {
    return force;
  }
  return step;
}

/// Unit vector of the lowest signed pair curvature, when it lies below
/// `-threshold`.
inline std::optional<Eigen::VectorXd>
pairNegativeMode(const Eigen::VectorXd &pos, KnownPair known,
                 ObjectiveFunction &objf, double threshold, double cutoff) {
  if (known == KnownPair::None || !(threshold > 0.0) || pos.size() < 6 ||
      pos.size() % 3 != 0) {
    return std::nullopt;
  }
  const Eigen::MatrixXd P =
      build(pos, "pair_full", knownPot(known), 3.0, 1.0, cutoff, objf);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(P);
  if (solver.info() != Eigen::Success || solver.eigenvalues().size() == 0) {
    return std::nullopt;
  }
  if (solver.eigenvalues()(0) >= -threshold) {
    return std::nullopt;
  }
  return solver.eigenvectors().col(0);
}

} // namespace pairhess
} // namespace eonc

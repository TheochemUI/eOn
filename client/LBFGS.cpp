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
// Two-loop L-BFGS (Nocedal 1980; Liu-Nocedal 1989; Nocedal-Wright 7.4).
// Optional Zhang-Xu modified secant, Li-Fukushima cautious update,
// Powell / Al-Baali damping, Packwood-Kermode pair preconditioner,
// Al-Baali extra-updates, and Grippo nonmonotone accept.

#include "eon/LBFGS.h"
#include "eon/GeometryAnalysis.h"
#include "eon/PairHessian.h"
#include "eon/SafeMath.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {
void maybeProjectRigid(Eigen::VectorXd &vec, const Eigen::VectorXd &pos,
                       bool enabled) {
  if (!enabled || pos.size() < 6 || pos.size() % 3 != 0) {
    return;
  }
  const long nat = pos.size() / 3;
  const AtomMatrix coords = Eigen::Map<const AtomMatrix>(pos.data(), nat, 3);
  eonc::geometry::projectOutRotTrans(vec, coords);
}

std::vector<int> twoLoopIndex(int nPairs, int extra) {
  std::vector<int> idx;
  idx.reserve(static_cast<size_t>(nPairs + std::max(extra, 0)));
  for (int i = 0; i < nPairs; ++i) {
    idx.push_back(i);
  }
  for (int k = 0; k < extra && nPairs > 0; ++k) {
    idx.push_back(nPairs - 1);
  }
  return idx;
}

void addPairBlock(Eigen::MatrixXd &P, int i, int j, const Eigen::Vector3d &rij,
                  double k_par, double k_perp) {
  const double r = rij.norm();
  if (r < 1.0e-14 || (!std::isfinite(k_par) && !std::isfinite(k_perp))) {
    return;
  }
  if (k_par == 0.0 && k_perp == 0.0) {
    return;
  }
  const Eigen::Vector3d u = rij / r;
  const Eigen::Matrix3d H =
      k_perp * Eigen::Matrix3d::Identity() +
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

void ljDerivs(double r, double eps, double sigma, double &vp, double &vpp) {
  const double inv = sigma / r;
  const double x6 = inv * inv * inv * inv * inv * inv;
  const double x12 = x6 * x6;
  vp = 4.0 * eps * (-12.0 * x12 + 6.0 * x6) / r;
  vpp = 4.0 * eps * (156.0 * x12 - 42.0 * x6) / (r * r);
}

void morseDerivs(double r, double D, double a, double re, double &vp,
                 double &vpp) {
  const double e1 = std::exp(-a * (r - re));
  const double e2 = e1 * e1;
  vp = 2.0 * D * a * (e1 - e2);
  vpp = 2.0 * D * a * a * (2.0 * e2 - e1);
}
} // namespace

Eigen::Vector3d LBFGS::micRij(const Eigen::VectorXd &pos, int i, int j) const {
  Eigen::Vector3d dr = pos.segment<3>(3 * i) - pos.segment<3>(3 * j);
  m_objf->minimumImage(dr);
  return dr;
}

bool LBFGS::usesPrecon() const {
  const std::string &p = m_optConfig.opts.lbfgs.precon;
  return p == "exp" || p == "c1" || p == "lindh" || p == "pair" ||
         p == "pair_abs" || p == "pair_full";
}

Eigen::MatrixXd LBFGS::buildPrecon(const Eigen::VectorXd &pos) const {
  const int n = static_cast<int>(pos.size());
  const int nat = n / 3;
  const auto &cfg = m_optConfig.opts.lbfgs;
  std::vector<double> nn(static_cast<size_t>(nat),
                         std::numeric_limits<double>::infinity());
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < nat; ++j) {
      if (i == j) {
        continue;
      }
      const double r = micRij(pos, i, j).norm();
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
  const double rcut = cfg.precon_rcut > 0.0 ? cfg.precon_rcut : 2.5 * r_nn;
  const double A = cfg.precon_A;
  const double mu = cfg.precon_mu;
  const std::string &kind = cfg.precon;
  const PotType pot = m_optConfig.potential;

  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < nat; ++i) {
    for (int j = i + 1; j < nat; ++j) {
      const Eigen::Vector3d rij = micRij(pos, i, j);
      const double r = rij.norm();
      if (r >= rcut || r < 1.0e-14) {
        continue;
      }
      if (kind == "pair") {
        double vp = 0.0;
        double vpp = 0.0;
        if (pot == PotType::MORSE_PT) {
          // eOn Morse_Pt / LAMMPS pair_coeff 0.7102 1.6047 2.897
          morseDerivs(r, 0.7102, 1.6047, 2.897, vp, vpp);
        } else {
          // rgpot LJ / LJCluster: 4 u0 ((psi/r)^12 - (psi/r)^6), u0=psi=1
          ljDerivs(r, 1.0, 1.0, vp, vpp);
        }
        // Mones 2018 Sci. Rep.: keep the positive-definite pair pieces.
        // pair_abs uses |V''| and |V'/r| so attractive LJ neighbours
        // still contribute (the well sits past the LJ inflection).
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
        // Lindh, Bernhardsson, Karlstrom, Malmqvist, Chem. Phys. Lett. 1995.
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

Eigen::VectorXd LBFGS::applyH0(const Eigen::VectorXd &q, double H0,
                               const Eigen::VectorXd &pos) const {
  const auto &cfg = m_optConfig.opts.lbfgs;
  if (usesPrecon() && pos.size() >= 6 && pos.size() % 3 == 0) {
    const Eigen::MatrixXd P = buildPrecon(pos);
    Eigen::LDLT<Eigen::MatrixXd> ldlt(P);
    if (ldlt.info() == Eigen::Success) {
      return ldlt.solve(q);
    }
    QUILL_LOG_DEBUG(m_log, "[LBFGS] Packwood P failed to factor, using H0 I");
  }
  return H0 * q;
}

Eigen::VectorXd LBFGS::hessianStep(double a_maxMove,
                                   const Eigen::VectorXd &a_f) {
  Eigen::VectorXd r = m_objf->getPositions();
  Eigen::VectorXd f = a_f;
  if (!usesPrecon()) {
    return eonc::helpers::maxAtomMotionAppliedV(
        m_optConfig.opts.lbfgs.inverse_curvature * f, a_maxMove);
  }
  maybeProjectRigid(f, r, m_optConfig.opts.lbfgs.project_rigid);
  const Eigen::VectorXd g = -f;
  Eigen::MatrixXd H = buildPrecon(r);
  Eigen::VectorXd d;
  if (m_optConfig.opts.lbfgs.step == "rfo") {
    // Banerjee, Adams, Simons, Shepard, J. Phys. Chem. 1985.
    // Baker, J. Comput. Chem. 1986. Lowest mode of the augmented Hessian.
    const int n = static_cast<int>(H.rows());
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    A.topLeftCorner(n, n) = H;
    A.col(n).head(n) = g;
    A.row(n).head(n) = g.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    if (es.info() != Eigen::Success) {
      return eonc::helpers::maxAtomMotionAppliedV(f, a_maxMove);
    }
    const Eigen::VectorXd v = es.eigenvectors().col(0);
    if (std::abs(v(n)) < 1.0e-14) {
      return eonc::helpers::maxAtomMotionAppliedV(f, a_maxMove);
    }
    d = v.head(n) / v(n);
  } else {
    // Regularized Newton: H + mu I with mu = max(0, eps - lambda_min).
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    if (es.info() != Eigen::Success) {
      return eonc::helpers::maxAtomMotionAppliedV(f, a_maxMove);
    }
    Eigen::VectorXd ev = es.eigenvalues();
    const double lmin = ev.minCoeff();
    const double mu = (lmin < 1.0e-8) ? (1.0e-8 - lmin) : 0.0;
    ev.array() += mu;
    d = -es.eigenvectors() *
        (es.eigenvectors().transpose() * g).cwiseQuotient(ev);
  }
  maybeProjectRigid(d, r, m_optConfig.opts.lbfgs.project_rigid);
  return eonc::helpers::maxAtomMotionAppliedV(d, a_maxMove);
}

Eigen::VectorXd LBFGS::getStep(double a_maxMove, const Eigen::VectorXd &a_f) {
  const std::string &step = m_optConfig.opts.lbfgs.step;
  if (step == "newton" || step == "rfo") {
    return hessianStep(a_maxMove, a_f);
  }
  double H0 = m_optConfig.opts.lbfgs.inverse_curvature;
  Eigen::VectorXd r = m_objf->getPositions();
  Eigen::VectorXd f = a_f;
  maybeProjectRigid(f, r, m_optConfig.opts.lbfgs.project_rigid);

  if (m_iteration > 0) {
    Eigen::VectorXd dr = m_objf->difference(r, m_rPrev);
    const Eigen::VectorXd y = m_fPrev - f;
    const double sy = dr.dot(y);
    const double yy = y.squaredNorm();
    const double ss = dr.squaredNorm();
    const double C = eonc::safemath::safe_div(yy, sy, -1.0);
    if (C < 0) {
      QUILL_LOG_DEBUG(
          m_log, "[LBFGS] Negative curvature: {:.4f} eV/A^2 take max move step",
          C);
      if (m_optConfig.opts.lbfgs.curvature == "reset") {
        reset();
      }
      return eonc::helpers::maxAtomMotionAppliedV(1000 * f, a_maxMove);
    }

    if (m_optConfig.opts.lbfgs.auto_scale && sy > 0.0 && yy > 0.0) {
      const double g_syyy = sy / yy;
      const double g_sssy = ss / sy;
      const std::string &h0 = m_optConfig.opts.lbfgs.h0;
      if (h0 == "ss_sy") {
        H0 = g_sssy;
      } else if (h0 == "adaptive") {
        H0 = std::min(g_syyy, g_sssy);
      } else {
        H0 = g_syyy;
      }
      if (m_optConfig.opts.lbfgs.max_inverse_curvature > 0.0) {
        H0 = std::min(H0, m_optConfig.opts.lbfgs.max_inverse_curvature);
      }
      QUILL_LOG_DEBUG(m_log, "[LBFGS] H0={:.4e} ({} scale)", H0, h0);
    }
  }

  const auto idx =
      twoLoopIndex(static_cast<int>(m_s.size()),
                   static_cast<int>(m_optConfig.opts.lbfgs.extra_updates));
  const int loopmax = static_cast<int>(idx.size());
  std::vector<double> a(static_cast<size_t>(loopmax));

  Eigen::VectorXd q = -f;

  for (int k = loopmax - 1; k >= 0; k--) {
    const int i = idx[static_cast<size_t>(k)];
    a[static_cast<size_t>(k)] = m_rho[static_cast<size_t>(i)] * m_s[static_cast<size_t>(i)].dot(q);
    q -= a[static_cast<size_t>(k)] * m_y[static_cast<size_t>(i)];
  }

  Eigen::VectorXd z = applyH0(q, H0, r);

  for (int k = 0; k < loopmax; k++) {
    const int i = idx[static_cast<size_t>(k)];
    const double b =
        m_rho[static_cast<size_t>(i)] * m_y[static_cast<size_t>(i)].dot(z);
    z += m_s[static_cast<size_t>(i)] * (a[static_cast<size_t>(k)] - b);
  }

  Eigen::VectorXd d = -z;

  double distance = eonc::helpers::maxAtomMotionV(d);
  if (distance >= a_maxMove && m_optConfig.opts.lbfgs.distance_reset) {
    QUILL_LOG_DEBUG(m_log,
                    "[LBFGS] reset memory, proposed step too large: {:.4f}",
                    distance);
    reset();
    return eonc::helpers::maxAtomMotionAppliedV(H0 * f, a_maxMove);
  }

  double vd = eonc::safemath::safe_normalized(d).dot(
      eonc::safemath::safe_normalized(f));
  if (vd > 1.0)
    vd = 1.0;
  if (vd < -1.0)
    vd = -1.0;
  double angle = eonc::safemath::safe_acos(vd) * (180.0 / eonc::helpers::pi);
  if (angle > 90.0 && m_optConfig.opts.lbfgs.angle_reset) {
    QUILL_LOG_DEBUG(m_log,
                    "[LBFGS] reset memory, angle between LBFGS angle and "
                    "force too large: {:.4f}",
                    angle);
    reset();
    return eonc::helpers::maxAtomMotionAppliedV(H0 * f, a_maxMove);
  }

  maybeProjectRigid(d, r, m_optConfig.opts.lbfgs.project_rigid);
  return eonc::helpers::maxAtomMotionAppliedV(d, a_maxMove);
}

void LBFGS::reset() {
  m_s.clear();
  m_y.clear();
  m_rho.clear();
  m_eHist.clear();
}

int LBFGS::update(const Eigen::VectorXd &a_r1, const Eigen::VectorXd &a_r0,
                  const Eigen::VectorXd &a_f1, const Eigen::VectorXd &a_f0,
                  double a_e1, double a_e0) {
  Eigen::VectorXd s0 = m_objf->difference(a_r1, a_r0);

  // y0 is the change in the gradient, not the force
  Eigen::VectorXd y0 = a_f0 - a_f1;
  double sy = s0.dot(y0);
  const auto &cfg = m_optConfig.opts.lbfgs;
  const std::string &curv = cfg.curvature;
  const double H0 = std::max(cfg.inverse_curvature, 1.0e-16);
  const double ss = s0.squaredNorm();

  if (cfg.secant == "zhangxu" && ss > LBFGS_EPS) {
    // Zhang, Deng, Chen, JOTA 102, 147 (1999); Zhang and Xu, JOTA 2001.
    // t = 6(f_k - f_{k+1}) + 3(g_k + g_{k+1})·s, ŷ = y + (t - y·s)/||s||^2 s.
    // Forces f = -g, so (g_k + g_{k+1})·s = -(f0 + f1)·s.
    const double t = 6.0 * (a_e0 - a_e1) - 3.0 * (a_f0 + a_f1).dot(s0);
    const double theta = t - sy;
    y0 += (theta / ss) * s0;
    sy = s0.dot(y0);
    QUILL_LOG_DEBUG(m_log, "[LBFGS] Zhang-Xu θ={:.4e} s·ŷ={:.4e}", theta, sy);
  }

  double sBs = ss / H0;
  if (usesPrecon() && a_r1.size() >= 6 && a_r1.size() % 3 == 0) {
    const Eigen::MatrixXd P = buildPrecon(a_r1);
    sBs = s0.dot(P * s0);
  }

  const double gnorm = a_f1.norm();
  const double sn = std::sqrt(ss);
  const double yn = y0.norm();

  if (curv == "cautious") {
    // Li and Fukushima, J. Comput. Appl. Math. 129, 15 (2001).
    const double thresh =
        cfg.cautious_eps * ss * std::pow(std::max(gnorm, 1.0e-30), cfg.cautious_alpha);
    if (sy < thresh) {
      QUILL_LOG_DEBUG(m_log,
                      "[LBFGS] Li-Fukushima skip, s·y={:.4e} < {:.4e}", sy,
                      thresh);
      return 0;
    }
  } else if (std::abs(sy) < LBFGS_EPS || (curv != "reset" && sy < 0.2 * sBs)) {
    if (curv == "skip") {
      QUILL_LOG_DEBUG(m_log, "[LBFGS] skip pair, s·y={:.4e}", sy);
      return 0;
    }
    if (curv == "damped" && sBs > sy) {
      const double theta = std::clamp(0.8 * sBs / (sBs - sy), 0.0, 1.0);
      const Eigen::VectorXd B0s =
          (usesPrecon() && a_r1.size() >= 6 && a_r1.size() % 3 == 0)
              ? Eigen::VectorXd(buildPrecon(a_r1) * s0)
              : Eigen::VectorXd(s0 / H0);
      y0 = theta * y0 + (1.0 - theta) * B0s;
      sy = s0.dot(y0);
      QUILL_LOG_DEBUG(m_log, "[LBFGS] Powell damp θ={:.3f} s·ŷ={:.4e}", theta,
                      sy);
    } else if (std::abs(sy) < LBFGS_EPS) {
      QUILL_LOG_WARNING(m_log,
                        "[LBFGS] s0.y0 too small ({:.4e}), resetting memory",
                        s0.dot(y0));
      reset();
      return 0;
    }
  }

  // Relative curvature skip (xtsci / compact-Hessian safeguard).
  if (!std::isfinite(sy) || sy <= 1.0e-8 * sn * yn) {
    QUILL_LOG_DEBUG(m_log, "[LBFGS] overlap skip, s·y={:.4e}", sy);
    return 0;
  }

  m_rho.push_back(eonc::safemath::safe_recip(sy, 0.0));
  m_s.push_back(std::move(s0));
  m_y.push_back(std::move(y0));

  if (static_cast<int>(m_s.size()) > m_memory) {
    m_s.pop_front();
    m_y.pop_front();
    m_rho.pop_front();
  }
  return 0;
}

int LBFGS::step(double a_maxMove) {
  int status = 0;
  Eigen::VectorXd r = m_objf->getPositions();
  Eigen::VectorXd f = -m_objf->getGradient();
  const double e0 = m_objf->getEnergy();

  if (m_iteration > 0) {
    status = update(r, m_rPrev, f, m_fPrev, e0, m_ePrev);
  }
  if (status < 0)
    return -1;

  Eigen::VectorXd dr = getStep(a_maxMove, f);
  constexpr double max_erise = 1.0e-8;
  double ref = e0;
  if (m_optConfig.opts.lbfgs.accept == "nonmonotone" && !m_eHist.empty()) {
    ref = *std::max_element(m_eHist.begin(), m_eHist.end());
  }
  double alpha = 1.0;
  bool accepted = false;
  double e_acc = e0;
  for (int dec = 0; dec < 10; ++dec) {
    m_objf->setPositions(r + alpha * dr);
    e_acc = m_objf->getEnergy();
    if (e_acc - ref <= max_erise) {
      accepted = true;
      break;
    }
    alpha *= 0.5;
  }
  if (!accepted) {
    m_objf->setPositions(r);
    reset();
    m_objf->setPositions(
        r + eonc::helpers::maxAtomMotionAppliedV(0.1 * f, a_maxMove));
    e_acc = m_objf->getEnergy();
  }

  m_eHist.push_back(e_acc);
  if (m_eHist.size() > 5) {
    m_eHist.pop_front();
  }

  m_rPrev = r;
  m_fPrev = f;
  m_ePrev = e0;

  m_iteration++;

  return m_objf->isConverged() ? 1 : 0;
}

int LBFGS::run(size_t a_maxSteps, double a_maxMove) {
  int status;
  while (!m_objf->isConverged() && m_iteration < a_maxSteps) {
    status = step(a_maxMove);
    if (status < 0)
      return -1;
  }
  return m_objf->isConverged() ? 1 : 0;
}

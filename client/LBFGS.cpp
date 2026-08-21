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
// Based on the LBFGS minimizer written in ASE.

#include "eon/LBFGS.h"
#include "eon/GeometryAnalysis.h"
#include "eon/SafeMath.h"

#include <algorithm>
#include <cmath>

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
} // namespace

Eigen::VectorXd LBFGS::getStep(double a_maxMove, const Eigen::VectorXd &a_f) {
  double H0 = m_optConfig.opts.lbfgs.inverse_curvature;
  Eigen::VectorXd r = m_objf->getPositions();
  Eigen::VectorXd f = a_f;
  maybeProjectRigid(f, r, m_optConfig.opts.lbfgs.project_rigid);

  if (m_iteration > 0) {
    Eigen::VectorXd dr = m_objf->difference(r, m_rPrev);
    double C = eonc::safemath::safe_div((m_fPrev - f).dot(m_fPrev - f),
                                        dr.dot(m_fPrev - f), -1.0);
    if (C < 0) {
      QUILL_LOG_DEBUG(
          m_log, "[LBFGS] Negative curvature: {:.4f} eV/A^2 take max move step",
          C);
      // Li-Fukushima / Powell: keep older pairs. reset is the ASE wipe.
      if (m_optConfig.opts.lbfgs.curvature == "reset") {
        reset();
      }
      return eonc::helpers::maxAtomMotionAppliedV(1000 * f, a_maxMove);
    }

    if (m_optConfig.opts.lbfgs.auto_scale) {
      H0 = eonc::safemath::safe_recip(C, -1.0);
      QUILL_LOG_DEBUG(m_log, "[LBFGS] Curvature: {:.4e} eV/A^2", C);
    }
  }

  int loopmax = m_s.size();
  std::vector<double> a(loopmax);

  Eigen::VectorXd q = -f;

  for (int i = loopmax - 1; i >= 0; i--) {
    a[i] = m_rho[i] * m_s[i].dot(q);
    q -= a[i] * m_y[i];
  }

  Eigen::VectorXd z = H0 * q;

  for (int i = 0; i < loopmax; i++) {
    double b = m_rho[i] * m_y[i].dot(z);
    z += m_s[i] * (a[i] - b);
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
}

int LBFGS::update(const Eigen::VectorXd &a_r1, const Eigen::VectorXd &a_r0,
                  const Eigen::VectorXd &a_f1, const Eigen::VectorXd &a_f0) {
  Eigen::VectorXd s0 = m_objf->difference(a_r1, a_r0);

  // y0 is the change in the gradient, not the force
  Eigen::VectorXd y0 = a_f0 - a_f1;
  double sy = s0.dot(y0);
  const std::string &curv = m_optConfig.opts.lbfgs.curvature;
  const double H0 = std::max(m_optConfig.opts.lbfgs.inverse_curvature, 1.0e-16);
  // B0 = I/H0, so s^T B0 s = ||s||^2 / H0 (Powell / Nocedal-Wright §18.3).
  const double sBs = s0.squaredNorm() / H0;

  if (std::abs(sy) < LBFGS_EPS || (curv != "reset" && sy < 0.2 * sBs)) {
    if (curv == "skip") {
      QUILL_LOG_DEBUG(m_log, "[LBFGS] Li-Fukushima skip, s·y={:.4e}", sy);
      return 0;
    }
    if (curv == "damped" && sBs > sy) {
      const double theta = std::clamp(0.8 * sBs / (sBs - sy), 0.0, 1.0);
      y0 = theta * y0 + (1.0 - theta) * (s0 / H0);
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

  if (std::abs(sy) < LBFGS_EPS) {
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
    status = update(r, m_rPrev, f, m_fPrev);
  }
  if (status < 0)
    return -1;

  Eigen::VectorXd dr = getStep(a_maxMove, f);
  // OPTIM MYLBFGS: no Wolfe search; shrink the step if the energy rises.
  constexpr double max_erise = 1.0e-8;
  double alpha = 1.0;
  bool accepted = false;
  for (int dec = 0; dec < 10; ++dec) {
    m_objf->setPositions(r + alpha * dr);
    if (m_objf->getEnergy() - e0 <= max_erise) {
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
  }

  m_rPrev = r;
  m_fPrev = f;

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

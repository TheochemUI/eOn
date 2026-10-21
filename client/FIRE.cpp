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
#include "eon/FIRE.h"
#include "eon/ExactReduction.h"
#include "eon/GeometryAnalysis.h"
#include "eon/HelperFunctions.h"

#include <algorithm>
#include <format>
#include <stdexcept>

namespace eonc {

void FIRE::noteModel(const Eigen::VectorXd &force) {
  if (m_modelChecked) {
    return;
  }
  m_modelChecked = true;
  const pairhess::PairMatch match =
      pairhess::matchKnownPair(m_objf->getPositions(), force, *m_objf);
  m_known = match.kind;
  m_cutoff = match.cutoff;
  if (m_known != pairhess::KnownPair::None) {
    QUILL_LOG_DEBUG(m_log, "[FIRE] pair metric from the objective force");
  }
}

Eigen::VectorXd FIRE::metric(const Eigen::VectorXd &force) const {
  if (m_known == pairhess::KnownPair::None) {
    return force;
  }
  return pairhess::pairNewtonForce(m_objf->getPositions(), force, m_known,
                                   *m_objf, m_cutoff);
}

bool FIRE::leaveSaddle(double a_maxMove) {
  bool moved = false;
  const double threshold = m_optConfig.opts.converged_force;
  while (m_known != pairhess::KnownPair::None && m_escapes < 4 &&
         m_objf->isConverged()) {
    const Eigen::VectorXd pos = m_objf->getPositions();
    const std::optional<Eigen::VectorXd> mode =
        pairhess::pairNegativeMode(pos, m_known, *m_objf, threshold, m_cutoff);
    if (!mode) {
      break;
    }
    Eigen::VectorXd step =
        eonc::geometry::maxAtomMotionAppliedV(*mode, a_maxMove);
    const Eigen::VectorXd force = -m_objf->getGradient();
    if (eonc::exact::dot(step, force) < 0.0) {
      step = -step;
    }
    m_objf->setPositions(pos + step);
    m_vel.setZero();
    m_alpha = m_alpha_start;
    m_N = 0;
    ++m_escapes;
    moved = true;
    QUILL_LOG_DEBUG(m_log, "[FIRE] leave negative pair curvature");
  }
  return moved;
}

int FIRE::step(double a_maxMove) {
  double P = 0;
  // Check convergence.
  if (m_objf->isConverged()) {
    return leaveSaddle(a_maxMove) ? 0 : 1;
  }

  // Velocity Verlet. When the force is a known pair force, integrate
  // P^{-1} f so the step has the length of a pair Newton step.
  Eigen::VectorXd f = -m_objf->getGradient();
  noteModel(f);
  Eigen::VectorXd x = m_objf->getPositions();
  const Eigen::VectorXd drive = metric(f);

  m_vel += drive * m_dt;
  Eigen::VectorXd dx = m_vel * m_dt;

  dx = eonc::geometry::maxAtomMotionAppliedV(dx, a_maxMove);
  m_objf->setPositions(x + dx);

  f = -m_objf->getGradient();
  if (m_objf->isConverged()) {
    return leaveSaddle(a_maxMove) ? 0 : 1;
  }
  const Eigen::VectorXd mixed = metric(f);
  const double fnorm = eonc::exact::norm(mixed);
  Eigen::VectorXd f_unit = Eigen::VectorXd::Zero(mixed.size());
  if (fnorm > 0.0) {
    f_unit = mixed / fnorm;
  }

  // FIRE
  P = eonc::exact::dot(mixed, m_vel);
  m_vel = (1 - m_alpha) * m_vel + m_alpha * f_unit * eonc::exact::norm(m_vel);
  QUILL_LOG_DEBUG(m_log,
                  "P: {:.4f}, v: {:.4f}, m_dt: {:.4f}, m_alpha: {:.4f}, N: {}",
                  P, m_vel.norm(), m_dt, m_alpha, m_N);
  if (P >= 0) {
    m_N++;
    if (m_N > m_N_min) {
      m_dt = std::min(m_dt * m_f_inc, m_dt_max);
      m_alpha = m_alpha * m_f_a;
    }
  } else {
    m_dt = m_dt * m_f_dec;
    m_vel.setZero();
    m_alpha = m_alpha_start;
    m_N = 0;
  }

  // add a sanity check on m_dt
  if (m_dt < 1e-6) {
    QUILL_LOG_CRITICAL(m_log, "[FIRE] [critical] m_dt is too small: {:.4f}",
                       m_dt);
    throw std::runtime_error(
        std::format("[FIRE] m_dt is too small: {:.4f}", m_dt));
  }

  m_iteration++;
  return m_objf->isConverged() ? 1 : 0;
}

int FIRE::run(size_t a_maxIterations, double a_maxMove) {
  while (!m_objf->isConverged() && m_iteration < a_maxIterations) {
    step(a_maxMove);
  }
  return m_objf->isConverged() ? 1 : 0;
}

} // namespace eonc

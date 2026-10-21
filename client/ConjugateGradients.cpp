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
#include "eon/ConjugateGradients.h"
#include "eon/ExactReduction.h"
#include "eon/SafeMath.h"

#include <cmath>

namespace eonc {

Eigen::VectorXd ConjugateGradients::getStep() {
  double a = std::abs(eonc::exact::dot(m_force, m_forceOld));
  double b = eonc::exact::squaredNorm(m_forceOld);
  double gamma = 0.0;
  if (a < 0.5 * b) {
    // Polak-Ribiere way to determine how much to mix in of old direction
    gamma = eonc::safemath::safe_div(
        eonc::exact::dot(m_force, m_force - m_forceOld), b, 0.0);
  } else {
    gamma = 0;
  }
  m_direction = m_force + gamma * m_directionOld;
  m_directionNorm = m_direction;
  m_directionNorm = eonc::exact::normalized(m_directionNorm);
  m_directionOld = m_direction;
  m_forceOld = m_force;

  // Only if value for max nr of iteration before reset
  if (m_optConfig.opts.cg.max_iter_before_reset > 0 &&
      m_optConfig.opts.cg.max_iter_before_reset <= m_cg_i) {
    m_cg_i = 0;
    m_forceOld.setZero();
    m_directionOld.setZero();
  }
  m_cg_i += 1;

  return m_direction;
}

int ConjugateGradients::step(double a_maxMove) {
  bool converged;
  if (m_optConfig.opts.cg.line_search) {
    converged = line_search(a_maxMove);
  } else {
    converged = single_step(a_maxMove);
  }
  if (converged)
    return 1;
  return 0;
}

int ConjugateGradients::line_search(double a_maxMove) {
  Eigen::VectorXd pos;
  Eigen::VectorXd posStep;
  Eigen::VectorXd forceBeforeStep;
  double stepSize;
  double projectedForce;
  double projectedForceBeforeStep;
  double curvature;

  forceBeforeStep = -m_objf->getGradient();
  m_force = forceBeforeStep;
  getStep();
  pos = m_objf->getPositions();
  projectedForceBeforeStep = eonc::exact::dot(m_force, m_directionNorm);

  // move system an infinitesimal step to determine the optimal step size along
  // the search line
  posStep = pos + m_directionNorm * m_optConfig.finiteDifference;
  m_objf->setPositions(posStep);
  m_force = -m_objf->getGradient(true);
  projectedForce = eonc::exact::dot(m_force, m_directionNorm);
  stepSize = m_optConfig.finiteDifference;

  int line_i = 0;
  do {
    // Determine curvature from last step (Secant method)
    curvature = std::abs(eonc::safemath::safe_div(
        projectedForceBeforeStep - projectedForce, stepSize, 0.0));
    // A negative max move is a length (bowl breakout), not a backward step.
    const double maxMove = std::fabs(a_maxMove);
    stepSize = eonc::safemath::safe_div(projectedForce, curvature, maxMove);

    if (maxMove < std::fabs(stepSize)) {
      // first part get the sign of stepSize
      stepSize = ((stepSize > 0) - (stepSize < 0)) * maxMove;
    }

    forceBeforeStep = m_force;
    projectedForceBeforeStep = projectedForce;

    pos += stepSize * m_directionNorm;
    m_objf->setPositions(pos);
    m_force = -m_objf->getGradient();
    projectedForce = eonc::exact::dot(m_force, m_directionNorm);

    line_i += 1;

    // Line search considered converged based in the ratio between the projected
    // force and the norm of the true force
  } while (m_optConfig.opts.cg.line_converged <
               std::abs(projectedForce) /
                   (std::sqrt(eonc::exact::dot(m_force, m_force) +
                              m_optConfig.opts.cg.line_converged)) &&
           line_i < m_optConfig.opts.cg.line_search_max_iter);
  if (m_objf->isConverged())
    return 1;
  return 0;
}

int ConjugateGradients::single_step(double a_maxMove) {
  Eigen::VectorXd pos;
  Eigen::VectorXd posStep;
  Eigen::VectorXd forceAfterStep;

  m_force = -m_objf->getGradient();
  pos = m_objf->getPositions();
  const Eigen::VectorXd posBefore = pos;
  getStep();

  if (!m_objf->reuseDisplacementCurvature()) {
    // Minimization and dimer translation. One finite-difference probe
    // along the search line, then a Newton length, every iteration.
    posStep = pos + m_directionNorm * m_optConfig.finiteDifference;
    m_objf->setPositions(posStep);
    forceAfterStep = -m_objf->getGradient(true);

    double projectedForce1 = eonc::exact::dot(m_force, m_directionNorm);
    double projectedForce2 = eonc::exact::dot(forceAfterStep, m_directionNorm);
    double curvature = eonc::safemath::safe_div(
        projectedForce1 - projectedForce2, m_optConfig.finiteDifference, 0.0);

    double stepSize = a_maxMove;
    if (curvature > 0.0) {
      stepSize =
          eonc::safemath::safe_div(projectedForce1, curvature, a_maxMove);
    }
    if (m_optConfig.bowlBreakout && a_maxMove < 0.0) {
      stepSize = -a_maxMove;
      a_maxMove = -a_maxMove;
    }
    if (!m_optConfig.opts.cg.no_overshooting) {
      if (m_optConfig.bowlBreakout) {
        pos += eonc::geometry::maxMotionAppliedV(stepSize * m_directionNorm,
                                                 a_maxMove);
      } else {
        pos += eonc::geometry::maxAtomMotionAppliedV(stepSize * m_directionNorm,
                                                     a_maxMove);
      }
      m_objf->setPositions(pos);
    } else {
      double passedMinimum = -1.;
      double forceChange = 0.;
      while (passedMinimum < 0.0 &&
             0.1 * std::abs(projectedForce1) < std::abs(projectedForce2)) {
        posStep = pos + eonc::geometry::maxAtomMotionAppliedV(
                            stepSize * m_directionNorm, a_maxMove);
        m_objf->setPositions(posStep);
        forceAfterStep = -m_objf->getGradient(true);
        projectedForce2 = eonc::exact::dot(forceAfterStep, m_directionNorm);
        passedMinimum = projectedForce1 * projectedForce2;
        if (passedMinimum < 0.0 &&
            0.1 * std::abs(projectedForce1) < std::abs(projectedForce2)) {
          forceChange = (projectedForce1 - projectedForce2);
          stepSize =
              eonc::safemath::safe_div(projectedForce1, forceChange, 0.0) *
              stepSize;
          QUILL_LOG_DEBUG(m_log, "Force changed {}, step size adjusted to {}",
                          forceChange, stepSize);
        }
      }
    }
    if (m_optConfig.opts.cg.knock_out_max_move) {
      if (stepSize >= a_maxMove) {
        m_directionOld.setZero();
        m_forceOld.setZero();
        QUILL_LOG_DEBUG(m_log, "Resetting the old search direction");
      }
    }
    return m_objf->isConverged() ? 1 : 0;
  }

  const double maxMove = std::fabs(a_maxMove);
  double projectedForce1 = eonc::exact::dot(m_force, m_directionNorm);
  if (projectedForce1 <= 0.0 && eonc::exact::squaredNorm(m_force) > 0.0) {
    // A restarted mix can point uphill. Step along the force.
    m_directionNorm = eonc::exact::normalized(m_force);
    projectedForce1 = eonc::exact::dot(m_force, m_directionNorm);
  }
  double stepSize = maxMove;

  if (m_optConfig.opts.cg.no_overshooting) {
    // The overshoot guard still probes. The stock band path does not.
    posStep = pos + m_directionNorm * m_optConfig.finiteDifference;
    m_objf->setPositions(posStep);
    forceAfterStep = -m_objf->getGradient(true);
    double projectedForce2 = eonc::exact::dot(forceAfterStep, m_directionNorm);
    double curvature = eonc::safemath::safe_div(
        projectedForce1 - projectedForce2, m_optConfig.finiteDifference, 0.0);
    if (curvature > 0.0) {
      stepSize =
          eonc::safemath::safe_div(projectedForce1, curvature, a_maxMove);
    }
    if (m_optConfig.bowlBreakout && a_maxMove < 0.0) {
      stepSize = -a_maxMove;
      a_maxMove = -a_maxMove;
    }
    double passedMinimum = -1.;
    double forceChange = 0.;
    while (passedMinimum < 0.0 &&
           0.1 * std::abs(projectedForce1) < std::abs(projectedForce2)) {
      posStep = pos + eonc::geometry::maxAtomMotionAppliedV(
                          stepSize * m_directionNorm, a_maxMove);
      m_objf->setPositions(posStep);
      forceAfterStep = -m_objf->getGradient(true);
      projectedForce2 = eonc::exact::dot(forceAfterStep, m_directionNorm);

      passedMinimum = projectedForce1 * projectedForce2;
      if (passedMinimum < 0.0 &&
          0.1 * std::abs(projectedForce1) < std::abs(projectedForce2)) {
        forceChange = (projectedForce1 - projectedForce2);
        stepSize = eonc::safemath::safe_div(projectedForce1, forceChange, 0.0) *
                   stepSize;
        QUILL_LOG_DEBUG(m_log, "Force changed {}, step size adjusted to {}",
                        forceChange, stepSize);
      }
    }
  } else {
    // First step: one finite-difference probe, then a Newton length.
    // Later steps: secant curvature of the previous real displacement.
    // A restart only drops the old search direction. The trust region
    // shrinks when that step raised the force, so a bad length is not
    // repeated at the max-move cap.
    double trust = maxMove;
    const double fd = m_optConfig.finiteDifference;
    if (!m_haveSecant && fd > 0.0) {
      posStep = pos + m_directionNorm * fd;
      m_objf->setPositions(posStep);
      forceAfterStep = -m_objf->getGradient(true);
      const double projectedForce2 =
          eonc::exact::dot(forceAfterStep, m_directionNorm);
      const double curvature = (projectedForce1 - projectedForce2) / fd;
      m_objf->setPositions(pos);
      if (curvature > 0.0 && projectedForce1 > 0.0) {
        stepSize = projectedForce1 / curvature;
      } else {
        stepSize = 0.25 * maxMove;
      }
    } else if (m_haveSecant) {
      const Eigen::VectorXd disp = m_objf->difference(pos, m_secantPos);
      const double ds = eonc::exact::dot(disp, m_secantDir);
      if (std::abs(ds) > 1e-14) {
        const double projNow = eonc::exact::dot(m_force, m_secantDir);
        const double curvature = (m_secantProj - projNow) / ds;
        if (curvature > 0.0 && projectedForce1 > 0.0) {
          stepSize = projectedForce1 / curvature;
        } else {
          stepSize = 0.25 * maxMove;
        }
      }
      if (m_secantForce.size() == m_force.size() && m_force.size() >= 6 &&
          m_force.size() % 3 == 0) {
        const long nAtoms = m_force.size() / 3;
        if (static_cast<long>(m_atomCurv.size()) != nAtoms) {
          m_atomCurv.assign(static_cast<size_t>(nAtoms), 0.0);
        }
        for (long a = 0; a < nAtoms; ++a) {
          const auto dsAtom = disp.segment(3 * a, 3);
          const double dsn = dsAtom.squaredNorm();
          if (dsn < 1e-12) {
            continue;
          }
          const auto df =
              m_secantForce.segment(3 * a, 3) - m_force.segment(3 * a, 3);
          const double c = df.dot(dsAtom) / dsn;
          if (c > 0.05) {
            m_atomCurv[static_cast<size_t>(a)] = c;
          }
        }
      }
      if (m_lastMotion > 0.0) {
        const double fNorm = std::sqrt(eonc::exact::squaredNorm(m_force));
        // Halve the trust when the force rose. A slower grow stalls the
        // stiff pairs at the iteration cap, so a falling force may double.
        const double grow = (fNorm > m_lastForceNorm * 1.05) ? 0.5 : 2.0;
        trust = std::min(maxMove, grow * m_lastMotion);
      }
    }

    Eigen::VectorXd delta = m_directionNorm;
    if (m_optConfig.bowlBreakout && a_maxMove < 0.0) {
      stepSize = maxMove;
      a_maxMove = maxMove;
      trust = maxMove;
      delta *= stepSize;
    } else {
      delta *= stepSize;
    }

    if (m_optConfig.bowlBreakout) {
      pos += eonc::geometry::maxMotionAppliedV(delta, a_maxMove);
    } else if (!m_atomCurv.empty() &&
               static_cast<long>(m_atomCurv.size()) * 3 == delta.size()) {
      // Per-atom Newton length, then a cap. A scalar step is stuck on
      // the stiffest atom and crawls along the soft ones.
      const long nAtoms = static_cast<long>(m_atomCurv.size());
      for (long a = 0; a < nAtoms; ++a) {
        const double c = m_atomCurv[static_cast<size_t>(a)];
        auto seg = delta.segment(3 * a, 3);
        const double motion = seg.norm();
        if (!(c > 0.0) || motion < 1e-14) {
          continue;
        }
        const double newton = m_force.segment(3 * a, 3).norm() / c;
        const double target = std::min(trust, std::max(newton, 0.0));
        const double scale = std::min(4.0, target / motion);
        if (scale > 0.0 && std::abs(scale - 1.0) > 1e-6) {
          seg *= scale;
        }
        const double capped = seg.norm();
        if (capped > trust) {
          seg *= trust / capped;
        }
      }
      pos += delta;
    } else {
      pos += eonc::geometry::maxAtomMotionAppliedV(delta, trust);
    }
    m_objf->setPositions(pos);
    m_lastForceNorm = std::sqrt(eonc::exact::squaredNorm(m_force));
    m_lastMotion =
        eonc::geometry::maxAtomMotionV(m_objf->difference(pos, posBefore));
    m_secantPos = posBefore;
    m_secantForce = m_force;
    m_secantDir = m_directionNorm;
    m_secantProj = projectedForce1;
    m_haveSecant = true;
  }
  if (m_optConfig.opts.cg.knock_out_max_move) {
    if (stepSize >= a_maxMove) {
      // knockout old search direction
      m_directionOld.setZero();
      m_forceOld.setZero();
      m_haveSecant = false;
      m_atomCurv.clear();
      QUILL_LOG_DEBUG(m_log, "Resetting the old search direction");
    }
  }

  return m_objf->isConverged() ? 1 : 0;
}

int ConjugateGradients::run(size_t a_maxIterations, double a_maxMove) {
  size_t iterations = 0;
  while (!m_objf->isConverged() && iterations < a_maxIterations) {
    step(a_maxMove);
    iterations++;
  }
  return m_objf->isConverged() ? 1 : 0;
}

} // namespace eonc

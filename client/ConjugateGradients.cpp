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
#include "ConjugateGradients.h"
#include "HelperFunctions.h"
#include "ObjectiveFunction.h"
namespace eonc {
VectorType ConjugateGradients::getStep() {
  double a = 0, b = 0, gamma = 0;
  a = std::fabs(m_force.dot(m_forceOld));
  b = m_forceOld.squaredNorm();
  if (a < 0.5 * b) {
    // Polak-Ribiere way to determine how much to mix in of old direction
    gamma = m_force.dot(m_force - m_forceOld) / b;
  } else {
    gamma = 0;
  }
  m_direction = m_force + gamma * m_directionOld;
  m_directionNorm = m_direction;
  m_directionNorm.normalize();
  m_directionOld = m_direction;
  m_forceOld = m_force;

  // Only if value for max nr of iteration before reset
  if ((m_p.max_iter_before_reset > 0) and (m_p.max_iter_before_reset <= m_cg_i))
  // or gamma == 0))
  {
    m_cg_i = 0;
    m_forceOld = m_objf.getPositions() * 0.0; // setZero
    m_directionOld = m_objf.getPositions() * 0.0;
    //        std::cout<<"reset\n";
  }
  m_cg_i += 1;

  return m_direction;
}

bool ConjugateGradients::stepImpl(double a_maxMove) {
  bool converged;
  if (m_p.line_search) {
    converged = line_search(a_maxMove);
  } else {
    converged = single_step(a_maxMove);
  }
  if (converged)
    return 1;
  return 0;
}

int ConjugateGradients::line_search(double a_maxMove) {
  VectorType pos;
  VectorType posStep;
  VectorType forceBeforeStep;
  double stepSize;
  double projectedForce;
  double projectedForceBeforeStep;
  double curvature;

  forceBeforeStep = -m_objf.getGradient();
  m_force = forceBeforeStep;
  getStep();
  pos = m_objf.getPositions();
  projectedForceBeforeStep = m_force.dot(m_directionNorm);

  // move system an infinitesimal step to determine the optimal step size along
  // the search line
  posStep = pos + m_directionNorm * m_p.finite_diff;
  const_cast<ObjectiveFunction &>(m_objf).setPositions(posStep);
  m_force = -m_objf.getGradient(true);
  projectedForce = m_force.dot(m_directionNorm);
  stepSize = m_p.finite_diff;

  int line_i = 0;
  do {
    // Determine curvature from last step (Secant method)
    curvature = fabs((projectedForceBeforeStep - projectedForce) / stepSize);
    stepSize = projectedForce / curvature;
    // stepSize = projectedForceBeforeStep / curvature;

    if (a_maxMove < fabs(stepSize)) {
      // first part get the sign of stepSize
      stepSize = ((stepSize > 0) - (stepSize < 0)) * a_maxMove;
    }

    forceBeforeStep = m_force;
    projectedForceBeforeStep = projectedForce;

    pos += stepSize * m_directionNorm;
    const_cast<ObjectiveFunction &>(m_objf).setPositions(pos);
    m_force = -m_objf.getGradient();
    projectedForce = m_force.dot(m_directionNorm);

    line_i += 1;

    // Line search considered converged based in the ratio between the projected
    // force and the norm of the true force
  } while (m_p.line_convergence <
               fabs(projectedForce) /
                   (sqrt(m_force.dot(m_force) + m_p.line_convergence)) and
           (line_i < m_p.max_line_serach_iter));
  //    return objf.isConverged();
  if (m_objf.isConverged())
    return 1;
  return 0;
}

int ConjugateGradients::single_step(double a_maxMove) {
  VectorType pos;
  VectorType posStep;
  VectorType forceAfterStep;

  m_force = -m_objf.getGradient();
  pos = m_objf.getPositions();
  getStep();

  // move system an infinitesimal step to determine the optimal step size along
  // the search line
  posStep = pos + m_directionNorm * m_p.finite_diff;
  const_cast<ObjectiveFunction &>(m_objf).setPositions(posStep);
  forceAfterStep = -m_objf.getGradient(true);

  // Determine curvature
  double projectedForce1 = m_force.dot(m_directionNorm);
  double projectedForce2 = forceAfterStep.dot(m_directionNorm);
  double curvature = (projectedForce1 - projectedForce2) / m_p.finite_diff;

  double stepSize = a_maxMove;

  if (curvature > 0.0) {
    stepSize = projectedForce1 / curvature;
  }

  if (m_p.saddle_bowl_breakout and a_maxMove < 0.0) {
    stepSize = -a_maxMove;
    a_maxMove = -a_maxMove;
  }

  if (!m_p.no_overshooting) {
    if (m_p.saddle_bowl_breakout) {
      // max displacement is based on system not single atom
      pos += helper_functions::maxMotionAppliedV(stepSize * m_directionNorm,
                                                 a_maxMove);
    } else {
      pos += helper_functions::maxAtomMotionAppliedV(stepSize * m_directionNorm,
                                                     a_maxMove);
    }
    const_cast<ObjectiveFunction &>(m_objf).setPositions(pos);
  } else {
    // negative if product of the projected forces before and after the step are
    // in opposite directions
    double passedMinimum = -1.;
    double forceChange = 0.;
    while (passedMinimum < 0. and
           (0.1 * fabs(projectedForce1) < fabs(projectedForce2))) {
      posStep = pos + helper_functions::maxAtomMotionAppliedV(
                          stepSize * m_directionNorm, a_maxMove);
      const_cast<ObjectiveFunction &>(m_objf).setPositions(posStep);
      forceAfterStep = -m_objf.getGradient(true);
      projectedForce2 = forceAfterStep.dot(m_directionNorm);

      passedMinimum = projectedForce1 * projectedForce2;
      if (passedMinimum < 0. and
          (0.1 * fabs(projectedForce1) < fabs(projectedForce2))) {
        forceChange = (projectedForce1 - projectedForce2);
        stepSize = (projectedForce1 / forceChange) * stepSize;
        SPDLOG_LOGGER_DEBUG(m_log, "Force changed {}, step size adjusted to {}",
                            forceChange, stepSize);
      }
    }
  }
  if (m_p.knock_out_max_move) {
    if (stepSize >= a_maxMove) {
      // knockout old search direction
      m_directionOld = m_objf.getPositions() * 0.0;
      m_forceOld = m_objf.getPositions() * 0.0;
      SPDLOG_LOGGER_DEBUG(m_log, "Resetting the old search direction");
    }
  }

  return m_objf.isConverged() ? 1 : 0;
}

bool ConjugateGradients::runOptImpl(size_t maxIterations, ScalarType maxMove) {
  int iterations = 0;
  while (!m_objf.isConverged() && iterations <= maxIterations) {
    stepImpl(maxMove);
    iterations++;
  }
  //    return objf->isConverged();
  if (m_objf.isConverged())
    return 1;
  return 0;
}

} // namespace eonc

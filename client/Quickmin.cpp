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
#include "eon/Quickmin.h"
#include "eon/HelperFunctions.h"

#include <cmath>


namespace eonc {

int Quickmin::step(double a_maxMove) {
  Eigen::VectorXd force = -m_objf->getGradient();
  const double fn = force.norm();
  if (!(fn > 0.0) || !std::isfinite(fn)) {
    ++m_iteration;
    return 1;
  }
  if (m_optConfig.opts.quickmin.steepest_descent) {
    m_vel.setZero();
  } else {
    if (m_vel.dot(force) < 0) {
      m_vel.setZero();
    } else {
      const Eigen::VectorXd f_unit = force / fn;
      m_vel = m_vel.dot(f_unit) * f_unit;
    }
  }

  m_vel += force * m_dt;
  Eigen::VectorXd dr = eonc::helpers::maxAtomMotionAppliedV(
      m_vel * m_dt,
      a_maxMove); // used to be m_optConfig.opts.max_time_step
  QUILL_LOG_INFO(m_log, "{} M_Vel.norm() is {}", m_iteration, m_vel.norm());
  m_objf->setPositions(m_objf->getPositions() + dr);
  m_iteration++;
  return m_objf->isConverged() ? 1 : 0;
}

int Quickmin::run(size_t a_maxSteps, double a_maxMove) {
  while (!m_objf->isConverged() && m_iteration < a_maxSteps) {
    step(a_maxMove);
  }
  return m_objf->isConverged() ? 1 : 0;
}

} // namespace eonc

#include "Quickmin.h"
#include "HelperFunctions.h"

int Quickmin::step(double a_maxMove) {
  Eigen::VectorXd force = -m_objf->getGradient();
  if (m_params->optQMSteepestDecent) {
    m_vel.setZero();
  } else {
    if (m_vel.dot(force) < 0) {
      m_vel.setZero();
    } else {
      Eigen::VectorXd f_unit = force.normalized();
      m_vel = m_vel.dot(f_unit) * f_unit;
    }
  }

  m_vel += force * m_dt;
  Eigen::VectorXd dr = helper_functions::maxAtomMotionAppliedV(
      m_vel * m_dt, a_maxMove); // used to be m_params->optMaxTimeStep
  SPDLOG_LOGGER_INFO(m_log, "{} M_Vel.norm() is {}", m_iteration,
                     m_vel.norm());
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

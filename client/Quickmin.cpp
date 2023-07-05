#include "Quickmin.h"
#include "HelperFunctions.h"

int Quickmin::step(double m_maxMove) {
  Eigen::VectorXd force = -m_objf->getGradient();
  if (m_params->optQMSteepestDecent) {
    m_vel.setZero();
  } else {
    if (m_vel.dot(force) < 0) {
      m_vel.setZero();
    } else {
      Eigen::VectorXd f_unit = force / force.norm();
      m_vel = m_vel.dot(f_unit) * f_unit;
    }
  }

  m_vel += force * m_dt;
  Eigen::VectorXd dr = helper_functions::maxAtomMotionAppliedV(
      m_vel * m_dt, m_maxMove); // used to be m_params->optMaxTimeStep
  SPDLOG_LOGGER_INFO(m_log, "{} M_Vel is {}", m_iteration,
                     fmt::streamed(m_vel));
  m_objf->setPositions(m_objf->getPositions() + dr);
  m_iteration++;
  return m_objf->isConverged();
}

int Quickmin::run(size_t m_maxSteps, double m_maxMove) {
  while (!m_objf->isConverged() && m_iteration < m_maxSteps) {
    step(m_maxMove);
  }
  return m_objf->isConverged();
}

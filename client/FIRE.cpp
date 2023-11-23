#include "FIRE.h"
#include "HelperFunctions.h"

int FIRE::step(double a_maxMove) {
  double P = 0;
  // Check convergence.
  if (m_objf->isConverged()) {
    return 1;
  }

  // Velocity Verlet
  Eigen::VectorXd f = -m_objf->getGradient();
  Eigen::VectorXd x = m_objf->getPositions();

  m_vel += f * m_dt;
  Eigen::VectorXd dx = m_vel * m_dt;

  dx = helper_functions::maxAtomMotionAppliedV(dx, m_max_move);
  m_objf->setPositions(x + dx);

  f = -m_objf->getGradient();
  Eigen::VectorXd f_unit = f / f.norm();

  // FIRE
  P = f.dot(m_vel);
  m_vel = (1 - m_alpha) * m_vel + m_alpha * f_unit * m_vel.norm();
  SPDLOG_LOGGER_DEBUG(
      m_log, "P: {:.4f}, v: {:.4f}, m_dt: {:.4f}, m_alpha: {:.4f}, N: {}", P,
      m_vel.norm(), m_dt, m_alpha, m_N);
  if (P >= 0) {
    m_N++;
    if (m_N > m_N_min) {
      m_dt = min(m_dt * m_f_inc, m_dt_max);
      m_alpha = m_alpha * m_f_a;
    }
  } else {
    m_dt = m_dt * m_f_dec;
    m_vel = m_vel * 0.0;
    m_alpha = m_alpha_start;
    m_N = 0;
  }

  // add a sanity check on m_dt
  if (m_dt < 1e-6) {
    SPDLOG_LOGGER_CRITICAL(m_log, "[FIRE] [critical] m_dt is too small: {:.4f}",
                           m_dt);
    std::exit(1);
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

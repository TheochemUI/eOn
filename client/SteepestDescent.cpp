
// Based on the SteepestDescent minimizer written in ASE.

#include "SteepestDescent.h"

int SteepestDescent::step(double a_maxMove) {
  Eigen::VectorXd r = m_objf->getPositions();
  Eigen::VectorXd f = -m_objf->getGradient();

  Eigen::VectorXd dr;
  double alpha = m_params->optSDAlpha;
  if (m_params->optSDTwoPoint == true && iteration > 0) {
    Eigen::VectorXd dx = r - m_rPrev;
    Eigen::VectorXd dg = -f + m_fPrev;
    alpha = dx.dot(dx) / dx.dot(dg);
    if (alpha < 0) {
      alpha = m_params->optSDAlpha;
    }
    SPDLOG_LOGGER_DEBUG(m_log, "[SD] alpha: {:.4e}", alpha);
  }

  dr = alpha * f;
  dr = helper_functions::maxAtomMotionAppliedV(dr, a_maxMove);

  m_objf->setPositions(r + dr);

  m_rPrev = r;
  m_fPrev = f;

  iteration++;

  return m_objf->isConverged() ? 1 : 0;
}

int SteepestDescent::run(size_t a_maxIteration, double a_maxMove) {
  while (!m_objf->isConverged() && iteration < a_maxIteration) {
    step(a_maxMove);
  }
  return m_objf->isConverged() ? 1 : 0;
}

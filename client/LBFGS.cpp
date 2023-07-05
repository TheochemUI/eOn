// Based on the LBFGS minimizer written in ASE.

#include "LBFGS.h"

Eigen::VectorXd LBFGS::getStep(double a_maxMove, Eigen::VectorXd a_f) {
  double H0 = m_params->optLBFGSInverseCurvature;
  Eigen::VectorXd r = m_objf->getPositions();

  if (m_iteration > 0) {
    Eigen::VectorXd dr = m_objf->difference(r, m_rPrev);
    // double C = dr.dot(fPrev-f)/dr.dot(dr);
    double C = (m_fPrev - a_f).dot(m_fPrev - a_f) / dr.dot(m_fPrev - a_f);
    if (C < 0) {
      SPDLOG_LOGGER_DEBUG(
          m_log, "[LBFGS] Negative curvature: {:.4f} eV/A^2 take max move step",
          C);
      reset();
      return helper_functions::maxAtomMotionAppliedV(1000 * a_f, a_maxMove);
    }

    if (m_params->optLBFGSAutoScale) {
      H0 = 1. / C;
      SPDLOG_LOGGER_DEBUG(m_log, "[LBFGS] Curvature: {:.4e} eV/A^2", C);
    }
  }

  if (m_iteration == 0 && m_params->optLBFGSAutoScale) {
    m_objf->setPositions(r + m_params->finiteDifference * a_f.normalized());
    Eigen::VectorXd dg = m_objf->getGradient(true) + a_f;
    double C = dg.dot(a_f.normalized()) / m_params->finiteDifference;
    H0 = 1.0 / C;
    m_objf->setPositions(r);
    if (H0 < 0) {
      SPDLOG_LOGGER_WARN(m_log,
                         "[LBFGS] Negative curvature calculated via FD: {:.4e} "
                         "eV/A^2, take max move step",
                         C);
      reset();
      return helper_functions::maxAtomMotionAppliedV(1000 * a_f, a_maxMove);
    } else {
      SPDLOG_LOGGER_DEBUG(
          m_log, "[LBFGS] Curvature calculated via FD: {:.4e} eV/A^2", C);
    }
  }

  int loopmax = m_s.size();
  double a[loopmax];

  Eigen::VectorXd q = -a_f;

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

  double distance = helper_functions::maxAtomMotionV(d);
  if (distance >= a_maxMove && m_params->optLBFGSDistanceReset) {
    SPDLOG_LOGGER_DEBUG(m_log,
                        "[LBFGS] reset memory, proposed step too large: {:.4f}",
                        distance);
    reset();
    return helper_functions::maxAtomMotionAppliedV(H0 * a_f, a_maxMove);
  }

  double vd = d.normalized().dot(a_f.normalized());
  if (vd > 1.0)
    vd = 1.0;
  if (vd < -1.0)
    vd = -1.0;
  double angle = acos(vd) * (180.0 / M_PI);
  if (angle > 90.0 && m_params->optLBFGSAngleReset) {
    SPDLOG_LOGGER_DEBUG(m_log,
                        "[LBFGS] reset memory, angle between LBFGS angle and "
                        "force too large: {:.4f}",
                        angle);
    reset();
    return helper_functions::maxAtomMotionAppliedV(H0 * a_f, a_maxMove);
  }

  return helper_functions::maxAtomMotionAppliedV(d, a_maxMove);
}

void LBFGS::reset(void) {
  m_s.clear();
  m_y.clear();
  m_rho.clear();
}

int LBFGS::update(Eigen::VectorXd a_r1, Eigen::VectorXd a_r0,
                  Eigen::VectorXd a_f1, Eigen::VectorXd a_f0) {
  Eigen::VectorXd s0 = m_objf->difference(a_r1, a_r0);

  // y0 is the change in the gradient, not the force
  Eigen::VectorXd y0 = a_f0 - a_f1;

  // GH: added to prevent crashing
  if (abs(s0.dot(y0)) < LBFGS_EPS) {
    SPDLOG_LOGGER_ERROR(m_log, "[LBFGS] error, s0.y0 is too small: {:.4f}",
                        s0.dot(y0));
    return -1;
  }

  m_s.push_back(s0);
  m_y.push_back(y0);
  m_rho.push_back(1.0 / (s0.dot(y0)));

  if ((int)m_s.size() > m_memory) {
    m_s.erase(m_s.begin());
    m_y.erase(m_y.begin());
    m_rho.erase(m_rho.begin());
  }
  return 0;
}

int LBFGS::step(double a_maxMove) {
  int status = 0;
  Eigen::VectorXd r = m_objf->getPositions();
  Eigen::VectorXd f = -m_objf->getGradient();

  if (m_iteration > 0) {
    status = update(r, m_rPrev, f, m_fPrev);
  }
  if (status < 0)
    return -1;

  Eigen::VectorXd dr = getStep(a_maxMove, f);

  m_objf->setPositions(r + dr);

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

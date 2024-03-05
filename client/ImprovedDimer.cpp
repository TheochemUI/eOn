// An implementation of Johannes Kästner and Paul Sherwood's improved dimer.
// An attempt to keep to the variable names in their 2008 paper has been made.

#include "ImprovedDimer.h"
#include "HelperFunctions.h"
#include "LowestEigenmode.h"

using namespace helper_functions;

const char ImprovedDimer::OPT_SD[] = "sd";
const char ImprovedDimer::OPT_CG[] = "cg";
const char ImprovedDimer::OPT_LBFGS[] = "lbfgs";

ImprovedDimer::ImprovedDimer(std::shared_ptr<Matter> matter,
                             std::shared_ptr<Parameters> params,
                             std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  x0 = std::make_shared<Matter>(pot, params);
  x1 = std::make_shared<Matter>(pot, params);
  *x0 = *matter;
  *x1 = *matter;
  tau.resize(3 * matter->numberOfAtoms());
  tau.setZero();
  totalForceCalls = 0;

  if (params->dimerOptMethod == OPT_CG) {
    init_cg = true;
  }
  log = spdlog::get("combi");
}

void ImprovedDimer::compute(std::shared_ptr<Matter> matter,
                            AtomMatrix initialDirectionAtomMatrix) {

  VectorXd initialDirection = VectorXd::Map(initialDirectionAtomMatrix.data(),
                                            3 * matter->numberOfAtoms());
  tau = initialDirection.array() * matter->getFreeV().array();
  tau = initialDirection;
  tau.normalize();
  *x0 = *matter;
  *x1 = *matter;
  VectorXd x0_r = x0->getPositionsV();

  x1->setPositionsV(x0_r + params->finiteDifference * tau);

  if (params->dimerOptMethod == OPT_LBFGS) {
    s.clear();
    y.clear();
    rho.clear();
    init_lbfgs = true;
  }

  // other vectors
  VectorXd x1_rp;
  VectorXd x1_r;
  VectorXd tau_prime;
  VectorXd tau_Old;
  VectorXd g1_prime;

  double delta = params->finiteDifference;
  double phi_tol = M_PI * (params->dimerConvergedAngle / 180.0);
  double phi_prime = 0.0;
  double phi_min = 0.0;

  statsRotations = 0;

  auto x1p = std::make_shared<Matter>(pot, params);

  // Melander, Laasonen, Jonsson, JCTC, 11(3), 1055–1062, 2015
  // http://doi.org/10.1021/ct501155k
  if (params->dimerRemoveRotation) {
    rotationRemove(MatrixXd::Map(x0_r.data(), x0->numberOfAtoms(), 3), x1);
    x1_r = x1->getPositionsV();
    tau = x1_r - x0_r;
    tau.normalize();
    x1_r = x0_r + tau * delta;
  }

  // Calculate the gradients on x0 and x1, g0 and g1, respectively.
  VectorXd g0 = -x0->getForcesV();
  VectorXd g1 = -x1->getForcesV();

  positions.clear();
  gradients.clear();

  positions.push_back(x0->getPositionsV());
  positions.push_back(x1->getPositionsV());
  gradients.push_back(g0);
  gradients.push_back(g1);

  do // while we have not reached phi_tol or maximum rotations.
  {

    // Calculate the rotational force, F_R.
    F_R = -2.0 * (g1 - g0) + 2.0 * ((g1 - g0).dot(tau)) * tau;

    statsTorque = F_R.norm() / (params->finiteDifference * 2.0);

    // Determine the step direction, theta
    if (params->dimerOptMethod == OPT_SD) // steepest descent
    {
      theta = F_R / F_R.norm();
    } else if (params->dimerOptMethod == OPT_CG) // conjugate gradients
    {
      if (init_cg) {
        init_cg = false;
        gamma = 0.0;
      } else {
        a = fabs(F_R.dot(F_R_Old));
        b = F_R_Old.squaredNorm();
        if (a < 0.5 * b) {
          gamma = F_R.dot(F_R - F_R_Old) / b;
        } else {
          gamma = 0.0;
        }
      }

      if (gamma == 0.0) {
        theta = F_R;
      } else {
        // theta = F_R + thetaOld * F_R_Old.norm() * gamma;
        // xph: use thetaOld before normalized, not F_R_Old
        theta = F_R + thetaOld * gamma;
      }

      theta = theta - theta.dot(tau) * tau;
      // xph: use thetaOld before normalized
      thetaOld = theta;
      theta.normalize();

      F_R_Old = F_R;
    } else if (params->dimerOptMethod == OPT_LBFGS) // quasi-newton
    {
      if (init_lbfgs == false) {
        // xph: s0 should the difference between tau and tau_Old, which are
        // normalized vectors. VectorXd s0 = x1->getPositionsV() - rPrev;
        VectorXd s0 = tau - tau_Old;
        s.push_back(s0);
        // xph: rescale the force; or rescale s0 = s0*delta
        VectorXd y0 = (F_R_Old - F_R) / delta;
        y.push_back(y0);
        rho.push_back(1.0 / (s0.dot(y0)));
      } else {
        init_lbfgs = false;
      }

      // xph: 60 is better than 10. The default H0 in ASE is 70.
      double H0 = 1. / 60.;

      int loopmax = s.size();
      double a[loopmax];

      VectorXd q = -F_R;

      for (int i = loopmax - 1; i >= 0; i--) {
        a[i] = rho[i] * s[i].dot(q);
        q -= a[i] * y[i];
      }

      VectorXd z = H0 * q;

      for (int i = 0; i < loopmax; i++) {
        double b = rho[i] * y[i].dot(z);
        z += s[i] * (a[i] - b);
      }

      double vd = -z.normalized().dot(F_R.normalized());
      if (vd > 1.0)
        vd = 1.0;
      if (vd < -1.0)
        vd = -1.0;
      double angle = acos(vd) * (180.0 / M_PI);

      // xph: larger angle is allowed to avoid frequent restart (was 70.0)
      if (angle > 87.0) {
        s.clear();
        y.clear();
        rho.clear();
        z = -F_R;
      }

      theta = -z.normalized();
      // xph:  rethogonalize theta to tau
      theta = theta - theta.dot(tau) * tau;
      theta.normalize();

      thetaOld = theta;
      F_R_Old = F_R;
      tau_Old = tau;
    }

    // Calculate the curvature along tau, C_tau.
    C_tau = (g1 - g0).dot(tau) / delta;

    // Calculate a rough estimate (phi_prime) of the optimum rotation angle.
    double d_C_tau_d_phi = 2.0 * (g1 - g0).dot(theta) / delta;
    phi_prime = -0.5 * atan(d_C_tau_d_phi / (2.0 * abs(C_tau)));
    statsAngle = phi_prime * (180.0 / M_PI);

    if (abs(phi_prime) > phi_tol) {
      double b1 = 0.5 * d_C_tau_d_phi;

      // Calculate g1_prime.
      x0_r = x0->getPositionsV();
      // xph: renormalize the new tangent after rotating phi_prime
      // x1_rp = x0_r + (tau * cos(phi_prime) + theta * sin(phi_prime)) * delta;
      tau_prime = tau * cos(phi_prime) + theta * sin(phi_prime);
      tau_prime = tau_prime.normalized();
      x1_rp = x0_r + tau_prime * delta;

      *x1p = *x1;
      x1p->setPositionsV(x1_rp);
      g1_prime = -x1p->getForcesV();

      // Update position and curvature histories
      positions.push_back(x1_rp);
      gradients.push_back(g1_prime);

      // Calculate C_tau_prime.
      // tau_prime = (x1_rp - x0_r) / (x1_rp - x0_r).norm(); //xph
      double C_tau_prime = (g1_prime - g0).dot(tau_prime) / delta;

      // Calculate phi_min.
      double a1 = (C_tau - C_tau_prime + b1 * sin(2.0 * phi_prime)) /
                  (1.0 - cos(2.0 * phi_prime));
      double a0 = 2.0 * (C_tau - a1);
      phi_min = 0.5 * atan(b1 / a1);

      // Determine the curvature for phi_min.
      double C_tau_min =
          0.5 * a0 + a1 * cos(2.0 * phi_min) + b1 * sin(2.0 * phi_min);

      // If the curvature is being maximized, push it over pi/2.
      if (C_tau_min > C_tau) {
        phi_min += M_PI * 0.5;
        C_tau_min =
            0.5 * a0 + a1 * cos(2.0 * phi_min) + b1 * sin(2.0 * phi_min);
      }

      // xph: for accurate LBFGS
      if (phi_min > M_PI * 0.5) {
        phi_min -= M_PI;
      }
      statsAngle = phi_min * (180.0 / M_PI);

      // Update x1, tau, and C_tau.
      // xph: normalize first
      tau = tau * cos(phi_min) + theta * sin(phi_min);
      tau = tau.normalized();
      x1_r = x0_r + tau * delta;

      // Melander, Laasonen, Jonsson, JCTC, 11(3), 1055–1062, 2015
      // http://doi.org/10.1021/ct501155k
      if (params->dimerRemoveRotation) {
        x1->setPositionsV(x1_r);
        rotationRemove(MatrixXd::Map(x0_r.data(), x0->numberOfAtoms(), 3), x1);
        x1_r = x1->getPositionsV();
        tau = x1_r - x0_r;
        tau.normalize();
        x1_r = x0_r + tau * delta;
      }

      // x1_r = x0_r + (tau * cos(phi_min) + theta * sin(phi_min)) * delta;
      x1->setPositionsV(x1_r);
      // tau = (x1_r - x0_r) / (x1_r - x0_r).norm();

      C_tau = C_tau_min;

      // Calculate the new g1.
      g1 = g1 * (sin(phi_prime - phi_min) / sin(phi_prime)) +
           g1_prime * (sin(phi_min) / sin(phi_prime)) +
           g0 * (1.0 - cos(phi_min) - sin(phi_min) * tan(phi_prime * 0.5));

      statsTorque = F_R.norm() / (2.0 * params->finiteDifference);
      statsRotations += 1;
      SPDLOG_LOGGER_INFO(
          log,
          "[IDimerRot]  -----   ---------   ----------   ------------------   "
          "{:9.4f}   {:7.3f}   {:6.2f}   {:4}",
          C_tau, statsTorque, statsAngle, statsRotations);
    } else {
      SPDLOG_LOGGER_INFO(
          log,
          "[IDimerRot]  -----   ---------   ----------   ------------------   "
          "{:9.4f}   {:7.3f}   ------   ----",
          C_tau, F_R.norm() / delta);
    }

  } while (abs(phi_prime) > abs(phi_tol) and abs(phi_min) > abs(phi_tol) and
           statsRotations < params->dimerRotationsMax);
}

double ImprovedDimer::getEigenvalue() { return C_tau; }

AtomMatrix ImprovedDimer::getEigenvector() {
  return MatrixXd::Map(tau.data(), x0->numberOfAtoms(), 3);
}

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
#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "eon/ConFileIO.h"
#include "eon/EonLogger.h"
#include "eon/GleThermostat.h"
#include "eon/HelperFunctions.h"
#include "eon/OHTSTJob.h"

namespace eonc {

// Optimal-hyperplane TST. Johánnesson and Jónsson, "Optimization of
// hyperplanar transition states", J. Chem. Phys. 115, 9644 (2001),
// doi:10.1063/1.1415499. The colored-noise thermostat option follows
// Ceriotti, Bussi and Parrinello, Phys. Rev. Lett. 102, 020601 (2009),
// doi:10.1103/PhysRevLett.102.020601.
//
// This lifts the harmonic-TST rate of adaptive KMC (Henkelman and
// Jónsson, J. Chem. Phys. 115, 9657 (2001), doi:10.1063/1.1415500;
// Pedersen and Jónsson, Math. Comput. Simul. 80, 1487 (2010),
// doi:10.1016/j.matcom.2009.02.010) to a free-energy barrier, which
// matters where a barrier is not large against kT: harmonic TST is
// quoted as reliable only for dE > ~5 kT.
//
// The dividing
// surface is the hyperplane n.(x - Gamma(s)) = 0 in the 3N_free
// configuration space, with Gamma(s) = R + s*u on the straight
// reactant->product guideline u = (P - R)/|P - R|. Sign conventions
// follow Appendix A: the plane translates along the driving force
// -<F.n> and its normal is driven along +<(n.F) R/(alpha |R|^2)>
// (verified to tilt the normal onto the ridge normal for a tilted-
// ridge model), so both damped Verlet integrations climb the free
// energy and settle where the mean forces vanish.

double OHTSTJob::uniformDraw() {
  // Deterministic LCG: the job must be reproducible for a given
  // random_seed across MPI farm workers.
  m_seedState = (1664525L * m_seedState + 1013904223L) & 0x7fffffffL;
  return static_cast<double>(m_seedState) / 2147483648.0;
}

double OHTSTJob::gaussDraw() {
  if (m_gaussHave) {
    m_gaussHave = false;
    return m_gaussSpare;
  }
  double u1 = std::max(uniformDraw(), 1e-12);
  double u2 = uniformDraw();
  double r = std::sqrt(-2.0 * std::log(u1));
  m_gaussSpare = r * std::sin(2.0 * helpers::pi * u2);
  m_gaussHave = true;
  return r * std::cos(2.0 * helpers::pi * u2);
}

void OHTSTJob::drawThermalVelocities(VectorXd &vel, const VectorXd *normal) {
  for (long k = 0; k < vel.size(); ++k) {
    vel[k] = std::sqrt(m_kbt / m_masses3N[k]) * gaussDraw();
  }
  if (normal != nullptr) {
    vel -= (*normal) * normal->dot(vel);
  }
}

bool OHTSTJob::symmetryReflect(const VectorXd &xR, VectorXd &x, VectorXd &v,
                               const VectorXd &xOld, const VectorXd *normal) {
  if (m_symDirs.size() < 2) {
    return false;
  }
  // Eq 15: distance from the configuration to each half-line
  // l_i = { R + t p_i, t >= 0 }.
  const VectorXd rel = x - xR;
  const double rel2 = rel.squaredNorm();
  double dPrimary = 0.0;
  size_t closest = 0;
  double dMin = 0.0;
  for (size_t i = 0; i < m_symDirs.size(); ++i) {
    const double proj = rel.dot(m_symDirs[i]);
    const double d2 = std::max(0.0, rel2 - proj * proj);
    const double d = std::sqrt(d2);
    if (i == 0) {
      dPrimary = d;
      dMin = d;
    } else if (d < dMin) {
      dMin = d;
      closest = i;
    }
  }
  if (closest == 0 || dPrimary <= dMin) {
    return false;
  }
  // Eqs 16-17: step back and reflect the velocity about the mirror
  // that maps p_1 onto the offending p_i; the component of the mirror
  // normal along the hyperplane normal is removed so the reflected
  // velocity stays within the plane.
  VectorXd q = m_symDirs[0] - m_symDirs[closest];
  if (normal != nullptr) {
    q -= (*normal) * normal->dot(q);
  }
  const double qn = q.norm();
  if (qn < 1e-12) {
    return false;
  }
  q /= qn;
  x = xOld;
  v -= 2.0 * v.dot(q) * q;
  return true;
}

OHTSTJob::PlaneAverages OHTSTJob::samplePlane(Matter &matter,
                                              const VectorXd &gamma,
                                              const VectorXd &normal) {
  const long equilSteps = params.oh_tst_options.equil_steps;
  const long sampleSteps = params.oh_tst_options.sample_steps;
  const double alphaRot = params.oh_tst_options.alpha_rot;

  // Constrain the current geometry exactly onto the plane.
  VectorXd x = matter.getPositionsFreeV();
  x -= normal * normal.dot(x - gamma);
  matter.setPositionsFreeV(x);

  VectorXd v(x.size());
  drawThermalVelocities(v, &normal);

  // Colored-noise option: an exact OU half-step before and after each
  // Verlet step (the auxiliary momenta start fresh per sampling
  // block); velocities re-project onto the plane after every kick.
  std::unique_ptr<GleThermostat> gle;
  if (params.oh_tst_options.thermostat == "gle") {
    const MatrixXd a =
        GleThermostat::loadDriftMatrix(params.oh_tst_options.gle_a_file);
    gle = std::make_unique<GleThermostat>(a, m_kbt, 0.5 * m_dt, x.size());
    if (!gle->valid()) {
      EONC_LOG_CRITICAL("OH-TST gle thermostat unusable (gle_a_file = {})",
                        params.oh_tst_options.gle_a_file);
      throw std::runtime_error("oh_tst: gle thermostat unusable");
    }
  }
  const auto gauss = [this]() { return gaussDraw(); };

  VectorXd f = matter.getForcesFreeV();
  VectorXd fPlane = f - normal * normal.dot(f);

  PlaneAverages avg;
  avg.rotNorm = VectorXd::Zero(x.size());
  avg.rotRaw = VectorXd::Zero(x.size());
  avg.pos = VectorXd::Zero(x.size());
  long nAccum = 0;

  VectorXd xPrev = x;
  for (long step = 0; step < equilSteps + sampleSteps; ++step) {
    // Velocity Verlet on the plane: forces and velocities projected,
    // positions corrected back onto the constraint (RATTLE for a
    // linear constraint is an exact projection).
    xPrev = x;
    if (gle) {
      gle->apply(v, m_masses3N, gauss);
      v -= normal * normal.dot(v);
    }
    VectorXd a = fPlane.cwiseQuotient(m_masses3N);
    x += m_dt * v + 0.5 * m_dt * m_dt * a;
    x -= normal * normal.dot(x - gamma);
    matter.setPositionsFreeV(x);
    f = matter.getForcesFreeV();
    fPlane = f - normal * normal.dot(f);
    VectorXd aNew = fPlane.cwiseQuotient(m_masses3N);
    v += 0.5 * m_dt * (a + aNew);
    v -= normal * normal.dot(v);
    // Eqs 14-17: keep the sampling in the primary product subregion.
    if (symmetryReflect(m_symXR, x, v, xPrev, &normal)) {
      x -= normal * normal.dot(x - gamma);
      matter.setPositionsFreeV(x);
      f = matter.getForcesFreeV();
      fPlane = f - normal * normal.dot(f);
    }
    if (gle) {
      gle->apply(v, m_masses3N, gauss);
      v -= normal * normal.dot(v);
    } else if (uniformDraw() < m_andersenProb) {
      // Andersen collisions keep the constrained ensemble canonical.
      drawThermalVelocities(v, &normal);
    }
    if (step >= equilSteps) {
      const double fn = normal.dot(f);
      const VectorXd arm = x - gamma;
      const double arm2 = arm.squaredNorm();
      avg.fn += fn;
      if (arm2 > 1e-16) {
        avg.rotNorm.noalias() += (fn / (alphaRot * arm2)) * arm;
      }
      avg.rotRaw.noalias() += fn * arm;
      avg.pos.noalias() += x;
      ++nAccum;
    }
  }
  if (nAccum > 0) {
    const double inv = 1.0 / static_cast<double>(nAccum);
    avg.fn *= inv;
    avg.rotNorm *= inv;
    avg.rotRaw *= inv;
    avg.pos *= inv;
  }
  return avg;
}

double OHTSTJob::reactantQRatio(Matter &matter, const VectorXd &gammaR,
                                const VectorXd &normal) {
  // Eq 22: Q^ZR/Q^R = (dt / t_tot) * sum over plane crossings of
  // 1 / |(r_{i+1} - r_i).n|. The trajectory is unconstrained,
  // thermostatted, and stays in the reactant basin by construction
  // (it starts there and the barrier is >> kT).
  const long steps = params.oh_tst_options.reactant_md_steps;
  const long equilSteps = params.oh_tst_options.equil_steps;
  VectorXd x = matter.getPositionsFreeV();
  VectorXd v(x.size());
  drawThermalVelocities(v, nullptr);
  std::unique_ptr<GleThermostat> gle;
  if (params.oh_tst_options.thermostat == "gle") {
    const MatrixXd a =
        GleThermostat::loadDriftMatrix(params.oh_tst_options.gle_a_file);
    gle = std::make_unique<GleThermostat>(a, m_kbt, 0.5 * m_dt, x.size());
    if (!gle->valid()) {
      throw std::runtime_error("oh_tst: gle thermostat unusable");
    }
  }
  const auto gauss = [this]() { return gaussDraw(); };
  VectorXd f = matter.getForcesFreeV();
  double side = normal.dot(x - gammaR);
  double crossingSum = 0.0;
  long counted = 0;
  for (long step = 0; step < equilSteps + steps; ++step) {
    const VectorXd xOld = x;
    if (gle) {
      gle->apply(v, m_masses3N, gauss);
    }
    VectorXd a = f.cwiseQuotient(m_masses3N);
    x += m_dt * v + 0.5 * m_dt * m_dt * a;
    matter.setPositionsFreeV(x);
    f = matter.getForcesFreeV();
    VectorXd aNew = f.cwiseQuotient(m_masses3N);
    v += 0.5 * m_dt * (a + aNew);
    if (symmetryReflect(m_symXR, x, v, xOld, nullptr)) {
      matter.setPositionsFreeV(x);
      f = matter.getForcesFreeV();
    }
    if (gle) {
      gle->apply(v, m_masses3N, gauss);
    } else if (uniformDraw() < m_andersenProb) {
      drawThermalVelocities(v, nullptr);
    }
    const double sideNew = normal.dot(x - gammaR);
    if (step >= equilSteps) {
      if (side * sideNew < 0.0) {
        const double proj = std::fabs(normal.dot(x - xOld));
        if (proj > 1e-14) {
          crossingSum += 1.0 / proj;
        }
      }
      ++counted;
    }
    side = sideNew;
  }
  return counted > 0 ? crossingSum / static_cast<double>(counted) : 0.0;
}

std::vector<std::string> OHTSTJob::run(void) {
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  if (!eonc::io::io_ok(
          reactant->con2matter(params.oh_tst_options.reactant_filename))) {
    EONC_LOG_CRITICAL("OH-TST failed to load {}",
                      params.oh_tst_options.reactant_filename);
    throw std::runtime_error("oh_tst: failed to load reactant");
  }
  if (!eonc::io::io_ok(
          product->con2matter(params.oh_tst_options.product_filename))) {
    EONC_LOG_CRITICAL("OH-TST failed to load {}",
                      params.oh_tst_options.product_filename);
    throw std::runtime_error("oh_tst: failed to load product");
  }

  const double temperature = params.main_options.temperature;
  EONC_LOG_INFO("[oh_tst] thermostat = {}{}", params.oh_tst_options.thermostat,
                params.oh_tst_options.thermostat == "gle"
                    ? std::string(" (drift: ") +
                          params.oh_tst_options.gle_a_file + ")"
                    : std::string());
  m_kbt = params.constants.kB * temperature;
  m_dt = params.oh_tst_options.time_step / params.constants.timeUnit;
  m_seedState = (params.main_options.randomSeed > 0)
                    ? params.main_options.randomSeed
                    : 12345;
  // Per-step collision probability from the Andersen collision period.
  const double tcol =
      params.thermostat_options.andersen_tcol_input / params.constants.timeUnit;
  m_andersenProb = (tcol > 0.0) ? std::min(1.0, m_dt / tcol) : 0.1;

  // Free-DOF mass vector (amu per coordinate).
  const long nAtoms = reactant->numberOfAtoms();
  auto masses = reactant->getMasses();
  std::vector<double> m3;
  m3.reserve(3 * nAtoms);
  for (long i = 0; i < nAtoms; ++i) {
    if (!reactant->getFixed(i)) {
      for (int j = 0; j < 3; ++j)
        m3.push_back(masses[i]);
    }
  }
  m_masses3N = VectorXd::Map(m3.data(), static_cast<long>(m3.size()));

  // Straight guideline in 3N space, minimum-image on the endpoint
  // difference so the plane never strides a cell boundary.
  const VectorXd xR = reactant->getPositionsFreeV();
  VectorXd diff = product->getPositionsFreeV() - xR;
  {
    AtomMatrix d(AtomMatrix::Map(diff.data(), diff.size() / 3, 3));
    // Rigid-translation alignment under PBC: AV-embedded state frames
    // ride the defect (drift = one hop vector, a/2<011> observed), so
    // atoms near cell boundaries min-image inconsistently until the
    // drift is gone. Iterate min-image -> de-drift to the fixed point;
    // the converged residual is the localized reaction coordinate.
    Eigen::RowVector3d total_drift = Eigen::RowVector3d::Zero();
    for (int pass = 0; pass < 6; ++pass) {
      d = reactant->pbc(d);
      const Eigen::RowVector3d drift =
          d.colwise().sum() / static_cast<double>(d.rows());
      d.rowwise() -= drift;
      total_drift += drift;
      if (drift.norm() < 1e-6) {
        break;
      }
    }
    EONC_LOG_INFO("[oh_tst] rigid drift removed: ({:.4f}, {:.4f}, "
                  "{:.4f}) A per atom",
                  total_drift[0], total_drift[1], total_drift[2]);
    diff = VectorXd::Map(d.data(), diff.size());
  }
  const double guideLen = diff.norm();
  EONC_LOG_INFO("[oh_tst] guideline length |P - R| = {:.4f} A over {} free "
                "DOF",
                guideLen, xR.size());
  // A sub-Angstrom guideline is an on-site shuffle (dumbbell rotation,
  // flicker partner): the plane progression has no room to climb and
  // the run would "converge" at s = 0 measuring nothing.
  if (guideLen < 0.5) {
    EONC_LOG_CRITICAL("[oh_tst] guideline too short ({:.4f} A) -- pick a "
                      "translation-class endpoint pair",
                      guideLen);
    throw std::runtime_error("oh_tst: degenerate guideline");
  }
  const VectorXd u = diff / guideLen;

  // Eqs 13-14: unit vectors to every symmetry-equivalent product;
  // index 0 is the primary product the guideline points to.
  m_symXR = xR;
  m_symDirs.clear();
  m_symDirs.push_back(u);
  if (!params.oh_tst_options.symmetry_products.empty()) {
    std::string rest = params.oh_tst_options.symmetry_products;
    while (!rest.empty()) {
      const auto comma = rest.find(',');
      std::string fname = rest.substr(0, comma);
      rest = (comma == std::string::npos) ? "" : rest.substr(comma + 1);
      if (fname.empty()) {
        continue;
      }
      Matter other(pot, params);
      if (!eonc::io::io_ok(other.con2matter(fname))) {
        EONC_LOG_CRITICAL("OH-TST failed to load symmetry product {}", fname);
        throw std::runtime_error("oh_tst: failed to load symmetry product");
      }
      VectorXd d = other.getPositionsFreeV() - xR;
      AtomMatrix dm(AtomMatrix::Map(d.data(), d.size() / 3, 3));
      dm = reactant->pbc(dm);
      d = VectorXd::Map(dm.data(), d.size());
      const double dn = d.norm();
      if (dn > 1e-8) {
        m_symDirs.push_back(d / dn);
      }
    }
    EONC_LOG_INFO("[oh_tst] symmetry restriction active over {} product "
                  "directions",
                  m_symDirs.size());
  }

  // Plane state: progression coordinate s, normal n, their conjugate
  // velocities, and the previous-iteration driving forces for the
  // two-force velocity Verlet updates (Eqs 6-7 and 9-10).
  double s = params.oh_tst_options.s_init * guideLen;
  double vS = 0.0;
  VectorXd n = u;
  VectorXd omega = VectorXd::Zero(n.size());
  const double mS = params.oh_tst_options.plane_mass;
  const double dtPlane = params.oh_tst_options.plane_time_step;
  const double dsMax = params.oh_tst_options.ds_max;
  const double dThetaMax = params.oh_tst_options.dtheta_max;
  const double fTol = params.oh_tst_options.force_tol;

  Matter walker(*reactant);

  // Reversible-work accumulators and the previous plane's averages
  // for the trapezoid rules of Eqs 18-19.
  double aTrans = 0.0, aRot = 0.0, aBest = 0.0, sBest = s;
  VectorXd nBest = n;
  bool havePrev = false;
  double fnPrev = 0.0;
  VectorXd rotRawPrev, posPrev, nPrev;
  double gSPrev = 0.0;
  VectorXd gRotPrev;
  int sideSign = 0; // sign of <F.n> at plane 1 (grooming, Sec IIE)

  // RAII: the divergence guard and samplePlane both throw out of the plane
  // loop below, so the handle has to close itself.
  std::ofstream prog("oh_tst_progression.dat");
  if (prog) {
    prog << "# plane  s/L  <F.n> (eV/A)  dA_trans (eV)  dA_rot (eV)  "
            "A (eV)  n.u\n";
  } else {
    EONC_LOG_ERROR("[oh_tst] cannot open oh_tst_progression.dat");
  }

  const bool scanMode = params.oh_tst_options.pmf_scan;
  const long nScan = std::max(2L, params.oh_tst_options.scan_planes);
  const double dsScan = scanMode ? (guideLen / (double)(nScan - 1)) : 0.0;
  const long nPlanes = scanMode ? nScan : params.oh_tst_options.max_planes;
  long plane = 0;
  bool converged = false;
  // Sec IIC guideline refinement: after the translational force first
  // changes sign the guideline is re-anchored every step through the
  // previous plane's average position along the current normal, so a
  // small rotational force no longer demands a huge s move (the
  // failure mode that made small-inertia mechanism discovery diverge).
  bool guidelineMoving = false;
  VectorXd gOrigin = xR;
  VectorXd gDir = u;
  long rotOnlySteps = 0;
  for (; plane < nPlanes; ++plane) {
    const VectorXd gamma = gOrigin + s * gDir;
    PlaneAverages avg = samplePlane(walker, gamma, n);

    // Driving forces: translation climbs against <F.n> (Eq 5 with the
    // reversed-force convention); the normal is driven along
    // +<(n.F) R/(alpha |R|^2)> (Appendix A, Eq A1), projected onto
    // the tangent space of the unit sphere.
    const double gS = -avg.fn / mS;
    VectorXd gRot = avg.rotNorm - n * n.dot(avg.rotNorm);

    if (plane == 0) {
      sideSign = (avg.fn < 0.0) ? -1 : 1;
    } else if (!guidelineMoving && ((avg.fn < 0.0) ? -1 : 1) != sideSign) {
      guidelineMoving = true;
      EONC_LOG_DEBUG("[oh_tst] plane {}: translational force changed "
                     "sign; guideline now follows <r> along the normal",
                     plane);
    }

    // Reversible work (Eqs 18-19), trapezoid between consecutive
    // planes: the translational path is the piecewise-linear track of
    // the average configuration <r>, and the rotational work pairs
    // the unnormalized <(n.F) R> with the actual change of the
    // normal. Grooming (Sec IIE): only planes on the reactant side of
    // the ridge (same sign of <F.n> as plane 1) contribute.
    if (havePrev) {
      // Scan mode integrates the whole reactant->product path; the
      // adaptive run grooms to the reactant side of the ridge only.
      const bool sameSide = scanMode || ((fnPrev < 0.0 ? -1 : 1) == sideSign &&
                                         (avg.fn < 0.0 ? -1 : 1) == sideSign);
      if (sameSide) {
        const VectorXd fParMean = 0.5 * (fnPrev * nPrev + avg.fn * n);
        aTrans += -fParMean.dot(avg.pos - posPrev);
        const VectorXd rotMean = 0.5 * (rotRawPrev + avg.rotRaw);
        aRot += rotMean.dot(n - nPrev);
      }
    }

    const double aTotal = aTrans + aRot;
    if (aTotal > aBest) {
      aBest = aTotal;
      sBest = s;
      nBest = n;
    }
    if (prog) {
      prog << std::format("{:6} {:10.6f} {:14.6e} {:12.6f} {:12.6f} {:12.6f} "
                          "{:10.6f}\n",
                          plane, s / guideLen, avg.fn, aTrans, aRot, aTotal,
                          n.dot(u));
      prog.flush();
    }
    EONC_LOG_DEBUG("[oh_tst] plane {} s/L {:.4f} <F.n> {:.4e} A {:.4f} eV",
                   plane, s / guideLen, avg.fn, aTotal);

    // Divergence guard: reversible work beyond any physical barrier
    // means the endpoints were not minimized (static relaxation
    // forces leak into <F.n>) or the plane is chasing a drifting
    // ensemble; 300 planes of that is pure waste.
    if (aTotal > params.oh_tst_options.max_delta_a) {
      EONC_LOG_CRITICAL(
          "[oh_tst] accumulated work {:.2f} eV exceeds max_delta_a "
          "{:.2f} eV at plane {} -- endpoints likely unminimized",
          aTotal, params.oh_tst_options.max_delta_a, plane);
      throw std::runtime_error("oh_tst: diverging reversible work");
    }
    // A stationary plane only counts as the variational maximum
    // after the progression has actually climbed: the reactant basin
    // bottom is also force-free, and a short guideline puts the first
    // plane inside it. Two thermal energies of accumulated reversible
    // work is the cheapest evidence of a ridge.
    if (!scanMode && plane > 2 && aTotal > 2.0 * m_kbt &&
        std::fabs(avg.fn) < fTol &&
        gRot.norm() * params.oh_tst_options.alpha_rot < fTol) {
      converged = true;
      // The converged plane is the optimal one even if sampling noise
      // put an earlier plane marginally higher.
      aBest = aTotal;
      sBest = s;
      nBest = n;
      if (prog) {
        prog << std::format("# converged at plane {}\n", plane);
      }
      break;
    }

    if (scanMode) {
      // Fixed-increment PMF scan: normal stays along the guideline,
      // s advances uniformly, and the walker is projected onto the
      // next constraint plane for samplePlane to re-equilibrate. The
      // A(s) profile is the PMF; aBest already tracks its maximum.
      fnPrev = avg.fn;
      nPrev = n;
      rotRawPrev = avg.rotRaw;
      posPrev = avg.pos;
      gSPrev = gS;
      gRotPrev = gRot;
      havePrev = true;
      s = (double)(plane + 1) * dsScan;
      const VectorXd gammaNext = xR + s * u;
      VectorXd xStart = walker.getPositionsFreeV();
      xStart -= u * (u.dot(xStart - gammaNext));
      walker.setPositionsFreeV(xStart);
      continue;
    }
    // Damped two-force Verlet on s (Eqs 6-7): the velocity is zeroed
    // when it opposes the driving force, so the plane settles at the
    // free-energy maximum instead of oscillating over it ("the
    // velocity of the plane is zeroed if the plane has gone past the
    // maximum free energy position").
    // Pure-rotation step (Sec IIC): right after a guideline change a
    // spiked rotational force is relaxed at fixed s before the plane
    // translates again.
    const bool rotationOnly =
        guidelineMoving &&
        gRot.norm() * params.oh_tst_options.alpha_rot > 5.0 * fTol &&
        rotOnlySteps < 10;
    if (rotationOnly) {
      ++rotOnlySteps;
    } else {
      rotOnlySteps = 0;
    }
    if (havePrev) {
      vS += 0.5 * dtPlane * (gS + gSPrev);
    } else {
      vS += dtPlane * gS;
    }
    if (vS * gS < 0.0)
      vS = 0.0;
    double ds = dtPlane * vS + 0.5 * dtPlane * dtPlane * gS;
    ds = std::clamp(ds, -dsMax, dsMax);
    if (rotationOnly) {
      ds = 0.0;
      vS = 0.0;
    }
    const double sNew =
        guidelineMoving ? s + ds : std::clamp(s + ds, 0.0, guideLen);

    // Damped rotation (Eqs 9-10 with the Appendix A driving): the
    // angular velocity keeps only its projection along the current
    // driving force while aligned with it, and is zeroed otherwise.
    if (havePrev && gRotPrev.size() == gRot.size()) {
      omega += 0.5 * dtPlane * (gRot + gRotPrev);
    } else {
      omega += dtPlane * gRot;
    }
    omega -= n * n.dot(omega);
    const double gNorm = gRot.norm();
    if (gNorm > 1e-14) {
      const VectorXd gHat = gRot / gNorm;
      const double along = omega.dot(gHat);
      if (along > 0.0) {
        omega = gHat * along;
      } else {
        omega.setZero();
      }
    } else {
      omega.setZero();
    }
    VectorXd dn = dtPlane * omega + 0.5 * dtPlane * dtPlane * gRot;
    const double dTheta = dn.norm();
    if (dTheta > dThetaMax) {
      dn *= dThetaMax / dTheta;
    }
    const VectorXd nOld = n;
    n = (n + dn).normalized();
    if (n.dot(gDir) < 0.0) {
      // Keep the normal pointing towards increasing s.
      n = -n;
    }

    // Eq 11-12 restart: place the walker in the new plane at the
    // rotated image of the average arm, rescaled so rotation does not
    // change the arm length.
    VectorXd arm = avg.pos - gamma;
    VectorXd armNew = arm - nOld * arm.dot(n - nOld);
    const double armLen = arm.norm();
    const double armNewLen = armNew.norm();
    if (armLen > 1e-12 && armNewLen > 1e-12) {
      armNew *= armLen / armNewLen;
    }
    VectorXd gammaNew;
    if (guidelineMoving) {
      // Re-anchor: the guideline passes through the previous plane's
      // average position along the CURRENT normal, and the plane
      // ADVANCES by this iteration's ds along the fresh line. (A
      // dead re-assignment here used to pin the plane at <r> forever:
      // the walker drifted downhill and the trans-work integral ran
      // away by ~20 eV per plane on the Cu doorway.)
      gOrigin = avg.pos;
      gDir = n;
      s = ds;
      gammaNew = gOrigin + s * gDir;
    } else {
      gammaNew = gOrigin + sNew * gDir;
    }
    VectorXd xStart = gammaNew + armNew;
    xStart -= n * n.dot(xStart - gammaNew);
    walker.setPositionsFreeV(xStart);

    fnPrev = avg.fn;
    nPrev = nOld;
    rotRawPrev = avg.rotRaw;
    posPrev = avg.pos;
    gSPrev = gS;
    gRotPrev = gRot;
    havePrev = true;
    if (!guidelineMoving) {
      s = sNew;
    }
  }
  bool progOk = prog.is_open();
  if (progOk) {
    prog.close();
    progOk = static_cast<bool>(prog);
    if (!progOk) {
      EONC_LOG_ERROR("[oh_tst] failed to write oh_tst_progression.dat");
    }
  }
  // A completed scan yields a valid PMF even with no interior
  // maximum; report it as converged so downstream tooling accepts
  // the profile (barrier = aBest, the max of A(s)).
  if (scanMode && plane >= nPlanes)
    converged = true;

  // Direction-dependent effective mass (Eq 24) and the one-sided
  // thermal flux factor sqrt(kBT / 2 pi mu) (Eq 23).
  const double mu = (m_masses3N.array() * nBest.array().square()).sum();
  const double vFlux = std::sqrt(m_kbt / (2.0 * helpers::pi * mu));

  // Q^ZR/Q^R at the FIRST plane of the progression (Z^R), whose
  // reversible work to the optimal plane is what aBest measures.
  Matter rWalker(*reactant);
  const VectorXd gammaR = xR + (params.oh_tst_options.s_init * guideLen) * u;
  const double qRatio = reactantQRatio(rWalker, gammaR, u);

  // Rate in internal units (1/internal-time), then SI.
  const double kInternal = vFlux * qRatio * std::exp(-aBest / m_kbt);
  const double kSI = kInternal / (params.constants.timeUnit * 1.0e-15);

  std::vector<std::string> returnFiles;
  std::ofstream out("results.dat");
  if (!out) {
    EONC_LOG_ERROR("[oh_tst] cannot open results.dat");
    throw std::runtime_error("oh_tst: cannot open results.dat");
  }
  out << "oh_tst job_type\n";
  out << std::format("{} converged\n", converged ? 1 : 0);
  out << std::format("{} planes_used\n", plane);
  out << std::format("{:.8f} free_energy_barrier_eV\n", aBest);
  out << std::format("{:.8f} delta_a_trans_eV\n", aTrans);
  out << std::format("{:.8f} delta_a_rot_eV\n", aRot);
  out << std::format("{:.8f} s_star_over_L\n", sBest / guideLen);
  out << std::format("{:.8f} guideline_length_A\n", guideLen);
  out << std::format("{:.8f} normal_overlap_with_guideline\n", nBest.dot(u));
  out << std::format("{:.8e} effective_mass_amu\n", mu);
  out << std::format("{:.8e} q_ratio_per_A\n", qRatio);
  out << std::format("{:.8e} rate_ohtst_per_s\n", kSI);
  out << std::format("{:.4f} temperature_K\n", temperature);
  out.close();
  if (!out) {
    EONC_LOG_ERROR("[oh_tst] failed to write results.dat");
    throw std::runtime_error("oh_tst: failed to write results.dat");
  }
  returnFiles.push_back("results.dat");
  if (progOk) {
    returnFiles.push_back("oh_tst_progression.dat");
  }
  EONC_LOG_INFO("[oh_tst] {} after {} planes: A = {:.4f} eV at s/L = {:.4f}, "
                "k = {:.4e} 1/s at {:.1f} K",
                converged ? "converged" : "max planes", plane, aBest,
                sBest / guideLen, kSI, temperature);
  return returnFiles;
}

} // namespace eonc

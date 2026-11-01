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
#include "eon/PrfoRefine.h"

#include "eon/EonLogger.h"

#include <stdexcept>
#include <string>
#include <vector>

#ifdef WITH_RGSADDLE
#include <rgsaddle.h>
#endif

namespace eonc {

double modeCurvature(const Matter &matter, const AtomMatrix &mode, double dr) {
  const AtomMatrix free = matter.getFree();
  const AtomMatrix t = mode.cwiseProduct(free);
  const double norm = t.norm();
  if (!(norm > 0.0) || !(dr > 0.0)) {
    throw std::invalid_argument("modeCurvature: zero mode or step");
  }
  const AtomMatrix unit = t / norm;
  Matter probe(matter);
  probe.setPositions(matter.getPositions() + dr * unit);
  return -(probe.getForces() - matter.getForces()).cwiseProduct(unit).sum() /
         dr;
}

// Application binary interface (ABI) minor 15 adds
// rgsaddle_sella_saddle_seed_mode.
#if defined(WITH_RGSADDLE) && RGSADDLE_ABI_MINOR >= 15
namespace {

// The session's coordinates are unwrapped, and Matter wraps what it is
// given, so the last evaluation is keyed on the session's own vector.
struct RefineSurface {
  Matter *matter;
  VectorXd mask;
  VectorXd start;
  VectorXd x;
  double energy;
  VectorXd gradient;
};

rgsaddle_status_t refineCallback(void *user, rgsaddle_surface_request_t *req) {
  auto *s = static_cast<RefineSurface *>(user);
  if (s == nullptr || req == nullptr || req->positions == nullptr ||
      req->energies == nullptr || req->gradients == nullptr ||
      req->n_images != 1) {
    return RGSADDLE_SURFACE_FAILED;
  }
  const auto dof = s->x.size();
  if (3 * req->n_atoms != dof) {
    return RGSADDLE_SHAPE;
  }
  try {
    const Eigen::Map<const VectorXd> x(req->positions, dof);
    if (x != s->x) {
      // A fixed coordinate of a free atom keeps its starting value.
      s->matter->setPositionsFreeV((s->mask.array() > 0.5).select(x, s->start));
      s->x = x;
      s->energy = s->matter->getPotentialEnergy();
      s->gradient = -s->matter->getForcesFreeV();
    }
    req->energies[0] = s->energy;
    Eigen::Map<VectorXd>(req->gradients, dof) = s->gradient;
    return RGSADDLE_OK;
  } catch (...) {
    return RGSADDLE_SURFACE_FAILED;
  }
}

void check(int rc, const char *what) {
  if (rc != RGSADDLE_OK) {
    throw std::runtime_error(
        std::string(what) + ": " +
        rgsaddle_status_name(static_cast<rgsaddle_status_t>(rc)));
  }
}

// The rows getPositionsFree() returns, in its order.
std::vector<long> freeAtoms(const Matter &matter) {
  const AtomMatrix free = matter.getFree();
  std::vector<long> rows;
  for (long i = 0; i < free.rows(); ++i) {
    if (free.row(i).sum() > 0.5) {
      rows.push_back(i);
    }
  }
  return rows;
}

} // namespace

PrfoRefineResult prfoRefineSaddle(Matter &matter, const AtomMatrix &mode,
                                  double curvature, const Parameters &params) {
  const auto &opt = params.saddle_search_options().prfo_refine;
  if (!(opt.force > 0.0) || opt.max_iterations < 1) {
    throw std::invalid_argument("prfo_refine_force must be positive and "
                                "prfo_refine_max_iterations at least 1");
  }
  rgsaddle_version_t stamp{};
  check(rgsaddle_abi_stamp(&stamp), "rgsaddle_abi_stamp");
  if (stamp.major != RGSADDLE_ABI_MAJOR || stamp.minor < 15) {
    throw std::runtime_error(
        "prfo_refine needs rgsaddle ABI 1.15 or newer at run time (loaded " +
        std::to_string(stamp.major) + "." + std::to_string(stamp.minor) + ")");
  }
  const std::vector<long> rows = freeAtoms(matter);
  const auto nFree = static_cast<long>(rows.size());
  const long callsBefore = matter.getForceCalls();
  const Matter start(matter);
  const double startForce = matter.maxForce();
  const AtomMatrix free = matter.getFree();
  AtomMatrix freeMode(nFree, 3);
  AtomMatrix freeMask(nFree, 3);
  VectorXd masses(nFree);
  for (long j = 0; j < nFree; ++j) {
    freeMask.row(j) = free.row(rows[j]);
    freeMode.row(j) = mode.row(rows[j]).cwiseProduct(free.row(rows[j]));
    masses[j] = matter.getMass(rows[j]);
  }

  RefineSurface surf{&matter,
                     VectorXd::Map(freeMask.data(), 3 * nFree),
                     matter.getPositionsFreeV(),
                     matter.getPositionsFreeV(),
                     matter.getPotentialEnergy(),
                     -matter.getForcesFreeV()};

  rgsaddle_sella_saddle_config_t config{};
  config.version = RGSADDLE_VERSION_INIT;
  config.force_tol = opt.force;
  config.force_gate = RGSADDLE_FORCE_MAX_ATOM;
  config.order = 1;
  // An empty chart: the default rigid quotient would project rotations
  // out of a periodic cell or a slab with fixed atoms.
  rgsaddle_constraints_config_t chartConfig{};
  chartConfig.version = RGSADDLE_VERSION_INIT;
  RgsaddleConstraints *chart = rgsaddle_constraints_create(&chartConfig, nFree);
  if (chart == nullptr) {
    throw std::runtime_error("rgsaddle_constraints_create failed");
  }
  RgsaddleSellaSaddle *session = rgsaddle_sella_saddle_create_on(
      &config, nFree, surf.start.data(), masses.data(), chart);
  rgsaddle_constraints_free(chart);
  if (session == nullptr) {
    throw std::runtime_error("rgsaddle_sella_saddle_create_on failed");
  }

  PrfoRefineResult result;
  try {
    // TS-BFGS (1), the transition-state Broyden-Fletcher-Goldfarb-Shanno
    // update. Plain BFGS (0) can stall at a saddle.
    check(rgsaddle_sella_saddle_set_hess_update(session, 1),
          "rgsaddle_sella_saddle_set_hess_update");
    check(rgsaddle_sella_saddle_seed_mode(session, freeMode.data(), curvature),
          "rgsaddle_sella_saddle_seed_mode");
    rgsaddle_report_t report{};
    while (!result.converged && result.steps < opt.max_iterations) {
      check(rgsaddle_sella_saddle_step(session, refineCallback, &surf, &report),
            "rgsaddle_sella_saddle_step");
      ++result.steps;
      result.converged = report.status == RGSADDLE_STATUS_CONVERGED;
    }
    VectorXd x(3 * nFree);
    check(rgsaddle_sella_saddle_position(session, x.data()),
          "rgsaddle_sella_saddle_position");
    if (x != surf.x) {
      matter.setPositionsFreeV((surf.mask.array() > 0.5).select(x, surf.start));
    }
  } catch (...) {
    rgsaddle_sella_saddle_free(session);
    throw;
  }
  rgsaddle_sella_saddle_free(session);

  result.maxForce = matter.maxForce();
  result.forceCalls = matter.getForceCalls() - callsBefore;
  if (result.maxForce > startForce) {
    matter.assignKeepingBias(start);
    result.restored = true;
    result.converged = false;
    result.maxForce = startForce;
  }
  result.energy = matter.getPotentialEnergy();
  EONC_LOG_INFO("[P-RFO] {} steps, {} force calls, max atom force {:.3e} -> "
                "{:.3e} eV/A{}",
                result.steps, result.forceCalls, startForce, result.maxForce,
                result.restored ? ", kept the starting point" : "");
  return result;
}

#else

PrfoRefineResult prfoRefineSaddle(Matter &, const AtomMatrix &, double,
                                  const Parameters &) {
  throw std::runtime_error("[Saddle Search] prfo_refine needs eOn built with "
                           "-Dwith_rgsaddle=true against rgsaddle ABI 1.15 or "
                           "newer");
}

#endif

} // namespace eonc

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
#include "eon/XtsciMinMode.h"

#include <rgsaddle.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace eonc {
namespace {

struct MinModeUser {
  // The host centre. It keeps its evaluation: a request at its own
  // coordinates reads the cache, and displaced probes go to `probe`.
  Matter *matter;
  Matter *probe;
  const AtomMatrix *fixed;
};

rgsaddle_status_t surfaceCallback(void *user, rgsaddle_surface_request_t *req) {
  auto *ctx = static_cast<MinModeUser *>(user);
  if (ctx == nullptr || ctx->matter == nullptr || ctx->probe == nullptr ||
      req == nullptr || req->positions == nullptr || req->energies == nullptr ||
      req->gradients == nullptr || req->n_images != 1) {
    return RGSADDLE_SURFACE_FAILED;
  }
  try {
    const auto nAtoms = ctx->matter->numberOfAtoms();
    if (req->n_atoms != nAtoms) {
      return RGSADDLE_SHAPE;
    }
    AtomMatrix pos = AtomMatrix::Map(req->positions, nAtoms, 3);
    if (ctx->fixed != nullptr) {
      for (int j = 0; j < nAtoms; ++j) {
        if (ctx->matter->getFixed(j)) {
          pos.row(j) = ctx->fixed->row(j);
        }
      }
    }
    Matter *at = ctx->matter;
    if (at->getPositions() != pos) {
      at = ctx->probe;
      if (at->getPositions() != pos) {
        at->setPositions(pos);
      }
    }
    req->energies[0] = at->getPotentialEnergy();
    const AtomMatrix &force = at->getForces();
    const auto dof = static_cast<Eigen::Index>(force.size());
    for (Eigen::Index k = 0; k < dof; ++k) {
      req->gradients[k] = -force.data()[k];
    }
    return RGSADDLE_OK;
  } catch (...) {
    return RGSADDLE_SURFACE_FAILED;
  }
}

void checkStatus(int rc, const char *what) {
  if (rc == RGSADDLE_OK) {
    return;
  }
  throw std::runtime_error(
      std::string(what) + ": " +
      rgsaddle_status_name(static_cast<rgsaddle_status_t>(rc)));
}

} // namespace

XtsciMinMode::XtsciMinMode(std::shared_ptr<Matter> matter,
                           const Parameters &params,
                           std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  if (matter) {
    m_eigenvector = AtomMatrix::Zero(matter->numberOfAtoms(), 3);
  }
}

void XtsciMinMode::compute(std::shared_ptr<Matter> matter,
                           AtomMatrix initialDirection) {
  if (!matter) {
    throw std::runtime_error("xtsci min-mode requires a geometry");
  }
  rgsaddle_version_t stamp{};
  checkStatus(rgsaddle_abi_stamp(&stamp), "rgsaddle_abi_stamp");
  if (stamp.major != RGSADDLE_ABI_MAJOR) {
    throw std::runtime_error("incompatible rgsaddle ABI");
  }
  const auto nAtoms = matter->numberOfAtoms();
  if (initialDirection.rows() != nAtoms || initialDirection.cols() != 3) {
    throw std::runtime_error("xtsci min-mode direction shape");
  }
  const AtomMatrix saved = matter->getPositions();
  m_fixed = saved;
  const auto &dim = params.dimer_options();
  rgsaddle_minmode_config_t config{};
  config.version = RGSADDLE_VERSION_INIT;
  config.flags = 0;
  // Dimer rotation. Lanczos is the same session with kind Lanczos when
  // the host dimer flag is off and the optimizer method is "lanczos".
  const bool lanczos = params.optimizer_options().xtsci.method == "lanczos";
  config.kind = lanczos ? RGSADDLE_MINMODE_LANCZOS : RGSADDLE_MINMODE_DIMER;
  config.method = params.optimizer_options().xtsci.method == "lbfgs"
                      ? RGSADDLE_METHOD_LBFGS
                      : RGSADDLE_METHOD_FIRE;
  config.dr = dim.rotation_angle > 0.0 ? dim.rotation_angle : 1.0e-3;
  config.rotation_tol = dim.torque_min;
  config.max_rotations = std::max<long>(1, dim.rotations_max);
  config.krylov_dim = 12;
  config.force_tol = 0.0;
  // Zero cap: rotation only. The climb uses the host optimizer.
  config.max_move = 0.0;

  RgsaddleMinMode *session = rgsaddle_minmode_create(
      &config, nAtoms, saved.data(), initialDirection.data());
  if (session == nullptr) {
    throw std::runtime_error("rgsaddle_minmode_create failed");
  }
  Matter probe(*matter);
  MinModeUser ctx{matter.get(), &probe, &m_fixed};
  const size_t callsBefore = matter->getPotentialCalls();
  rgsaddle_report_t report{};
  const int rc = rgsaddle_minmode_step(session, surfaceCallback, &ctx, &report);
  std::vector<double> mode(static_cast<size_t>(3 * nAtoms), 0.0);
  const int modeRc = rgsaddle_minmode_mode(session, mode.data());
  rgsaddle_minmode_free(session);
  checkStatus(rc, "rgsaddle_minmode_step");
  checkStatus(modeRc, "rgsaddle_minmode_mode");
  m_eigenvector = AtomMatrix::Map(mode.data(), nAtoms, 3);
  m_eigenvalue = report.curvature;
  statsCurvature = report.curvature;
  statsRotations = report.rotations;
  totalIterations += 1;
  totalForceCalls +=
      static_cast<long>(matter->getPotentialCalls() - callsBefore);
}

double XtsciMinMode::getEigenvalue() { return m_eigenvalue; }

AtomMatrix XtsciMinMode::getEigenvector() { return m_eigenvector; }

} // namespace eonc

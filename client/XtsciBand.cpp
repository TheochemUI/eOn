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
#include "eon/XtsciBand.h"

#include "eon/NudgedElasticBand.h"

#include <rgsaddle.h>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace eonc {

int xtsciBandSurface(void *user, void *request);

namespace {

int surfaceCallback(void *user, rgsaddle_surface_request_t *req) {
  auto *self = static_cast<XtsciBand *>(user);
  return self == nullptr ? -1 : xtsciBandSurface(self, req);
}

void checkStatus(int rc, const char *what) {
  if (rc == RGSADDLE_OK) {
    return;
  }
  throw std::runtime_error(std::string(what) + ": " +
                           rgsaddle_status_name(rc));
}

int32_t tangentKind(const neb_options_t &neb) {
  return neb.climbing_image.use_old_tangent ? RGSADDLE_TANGENT_SIMPLE
                                            : RGSADDLE_TANGENT_IMPROVED;
}

int32_t springKind(const neb_options_t &neb) {
  if (neb.spring.om.enabled) {
    return RGSADDLE_SPRING_ONSAGER_MACHLUP;
  }
  if (neb.spring.weighting.enabled) {
    return RGSADDLE_SPRING_WEIGHTED;
  }
  return RGSADDLE_SPRING_UNIFORM;
}

int32_t projectionKind(const neb_options_t &neb) {
  if (neb.spring.use_elastic_band) {
    return RGSADDLE_PROJECTION_PLAIN_EB;
  }
  if (neb.spring.doubly_nudged) {
    return RGSADDLE_PROJECTION_DNEB;
  }
  return RGSADDLE_PROJECTION_NEB;
}

int32_t methodKind(const Parameters::optimizer_options_t &opt) {
  const std::string &name = opt.xtsci.method;
  if (name == "lbfgs") {
    return RGSADDLE_METHOD_LBFGS;
  }
  return RGSADDLE_METHOD_FIRE;
}

void pinFixed(const AtomMatrix &ref, Matter &image, AtomMatrix *pos) {
  bool any = false;
  for (int j = 0; j < image.numberOfAtoms(); ++j) {
    if (!image.getFixed(j)) {
      continue;
    }
    pos->row(j) = ref.row(j);
    any = true;
  }
  if (any) {
    image.setPositions(*pos);
  }
}

} // namespace

int xtsciBandSurface(void *user, void *request) {
  auto *self = static_cast<XtsciBand *>(user);
  auto *req = static_cast<rgsaddle_surface_request_t *>(request);
  if (self == nullptr || req == nullptr || req->positions == nullptr ||
      req->energies == nullptr || req->gradients == nullptr) {
    return -1;
  }
  try {
    auto *neb = self->m_neb;
    const auto nImages = static_cast<int64_t>(neb->numImages + 2);
    const auto nAtoms = static_cast<int64_t>(neb->atoms);
    if (req->n_images != nImages || req->n_atoms != nAtoms) {
      return -1;
    }
    if (req->version.major != RGSADDLE_ABI_MAJOR) {
      return RGSADDLE_ABI_MISMATCH;
    }
    const auto dof = static_cast<Eigen::Index>(3 * nAtoms);
    for (int64_t i = 0; i < nImages; ++i) {
      AtomMatrix pos = AtomMatrix::Map(req->positions + i * dof, nAtoms, 3);
      neb->path[static_cast<size_t>(i)]->setPositions(pos);
      if (static_cast<size_t>(i) < self->m_fixed.size()) {
        pinFixed(self->m_fixed[static_cast<size_t>(i)], *neb->path[static_cast<size_t>(i)],
                 &pos);
      }
      req->energies[i] = neb->path[i]->getPotentialEnergy();
      const AtomMatrix &force = neb->path[i]->getForces();
      double *grad = req->gradients + i * dof;
      const double *src = force.data();
      for (Eigen::Index k = 0; k < dof; ++k) {
        grad[k] = -src[k];
      }
    }
    return 0;
  } catch (...) {
    return -1;
  }
}

XtsciBand::XtsciBand(NudgedElasticBand &neb, const Parameters &params)
    : m_neb{&neb},
      m_maxMove{params.optimizer_options().max_move} {
  rgsaddle_version_t stamp{};
  checkStatus(rgsaddle_abi_stamp(&stamp), "rgsaddle_abi_stamp");
  if (stamp.major != RGSADDLE_ABI_MAJOR) {
    throw std::runtime_error("incompatible rgsaddle ABI");
  }
  const auto &nebOpt = params.neb_options();
  const auto nImages = static_cast<int64_t>(neb.numImages + 2);
  const auto nAtoms = static_cast<int64_t>(neb.atoms);
  if (nImages < 3 || nAtoms < 1) {
    throw std::runtime_error("rgsaddle band needs at least 3 images");
  }
  if (springKind(nebOpt) == RGSADDLE_SPRING_WEIGHTED) {
    m_springKs.assign(static_cast<size_t>(nImages - 1),
                      nebOpt.spring.constant);
  }
  rgsaddle_band_config_t config{};
  config.version = RGSADDLE_VERSION_INIT;
  config.flags = 0;
  config.tangent = tangentKind(nebOpt);
  config.spring = springKind(nebOpt);
  config.projection = projectionKind(nebOpt);
  config.method = methodKind(params.optimizer_options());
  config.spring_k = nebOpt.spring.om.enabled
                        ? nebOpt.spring.constant * nebOpt.spring.om.k_scale
                        : nebOpt.spring.constant;
  config.spring_ks = m_springKs.empty() ? nullptr : m_springKs.data();
  if (nebOpt.climbing_image.enabled) {
    config.ci_trigger_factor = nebOpt.climbing_image.trigger_factor > 0.0
                                   ? nebOpt.climbing_image.trigger_factor
                                   : 1.0e-300;
    config.ci_trigger_force = nebOpt.climbing_image.trigger_force;
  } else {
    config.ci_trigger_factor = 0.0;
    config.ci_trigger_force = 0.0;
  }
  if (!neb.path.empty() && neb.path[0]->getPeriodic()) {
    const Matrix3d cell = neb.path[0]->getCell();
    m_cell.assign(cell.data(), cell.data() + 9);
    config.cell = m_cell.data();
  }
  config.force_tol = nebOpt.force_tolerance;
  config.max_move = m_maxMove > 0.0 ? m_maxMove : 0.2;
  config.memory = std::max<long>(1, params.optimizer_options().lbfgs.memory);

  std::vector<double> positions(
      static_cast<size_t>(nImages * 3 * nAtoms), 0.0);
  for (int64_t i = 0; i < nImages; ++i) {
    const AtomMatrix &pos = neb.path[i]->getPositions();
    std::copy(pos.data(), pos.data() + pos.size(),
              positions.begin() + static_cast<std::ptrdiff_t>(i * pos.size()));
  }
  m_band = rgsaddle_band_create(&config, nImages, nAtoms, positions.data());
  if (m_band == nullptr) {
    throw std::runtime_error("rgsaddle_band_create failed");
  }
}

XtsciBand::~XtsciBand() { rgsaddle_band_free(m_band); }

void XtsciBand::syncFromPath() {
  const auto nImages = static_cast<int64_t>(m_neb->numImages + 2);
  const auto nAtoms = static_cast<int64_t>(m_neb->atoms);
  std::vector<double> positions(
      static_cast<size_t>(nImages * 3 * nAtoms), 0.0);
  for (int64_t i = 0; i < nImages; ++i) {
    const AtomMatrix &pos = m_neb->path[i]->getPositions();
    std::copy(pos.data(), pos.data() + pos.size(),
              positions.begin() + static_cast<std::ptrdiff_t>(i * pos.size()));
  }
  checkStatus(rgsaddle_band_set_positions(m_band, positions.data()),
              "rgsaddle_band_set_positions");
}

void XtsciBand::reset() {
  checkStatus(rgsaddle_band_reset(m_band), "rgsaddle_band_reset");
}

void XtsciBand::step(double maxMove) {
  if (maxMove > 0.0) {
    m_maxMove = maxMove;
  }
  const auto nImages = static_cast<int64_t>(m_neb->numImages + 2);
  m_fixed.resize(static_cast<size_t>(nImages));
  for (int64_t i = 0; i < nImages; ++i) {
    m_fixed[static_cast<size_t>(i)] = m_neb->path[i]->getPositions();
  }
  syncFromPath();
  rgsaddle_report_t report{};
  checkStatus(rgsaddle_band_step(m_band, surfaceCallback, this, &report),
              "rgsaddle_band_step");
  std::vector<double> positions(m_fixed.size() * m_fixed[0].size(), 0.0);
  checkStatus(rgsaddle_band_positions(m_band, positions.data()),
              "rgsaddle_band_positions");
  const auto dof = m_fixed[0].size();
  for (size_t i = 0; i < m_fixed.size(); ++i) {
    AtomMatrix pos =
        AtomMatrix::Map(positions.data() + static_cast<std::ptrdiff_t>(i * dof),
                        m_neb->atoms, 3);
    // Endpoints are fixed by the session. Interior fixed atoms stay put.
    if (i == 0 || i + 1 == m_fixed.size()) {
      continue;
    }
    m_neb->path[i]->setPositions(pos);
    pinFixed(m_fixed[i], *m_neb->path[i], &pos);
  }
  m_neb->movedAfterForceCall = true;
}

} // namespace eonc

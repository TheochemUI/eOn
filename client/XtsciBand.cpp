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
#include <vector>

namespace eonc {

int xtsciBandSurface(void *user, void *request);

namespace {

rgsaddle_status_t surfaceCallback(void *user, rgsaddle_surface_request_t *req) {
  auto *self = static_cast<XtsciBand *>(user);
  return static_cast<rgsaddle_status_t>(
      self == nullptr ? RGSADDLE_SURFACE_FAILED : xtsciBandSurface(self, req));
}

void checkStatus(int rc, const char *what) {
  if (rc == RGSADDLE_OK) {
    return;
  }
  throw std::runtime_error(
      std::string(what) + ": " +
      rgsaddle_status_name(static_cast<rgsaddle_status_t>(rc)));
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

// Fixed atoms keep the coordinates they had before the step.
void pinRows(const AtomMatrix &ref, const Matter &image, AtomMatrix *pos) {
  for (int j = 0; j < image.numberOfAtoms(); ++j) {
    if (image.getFixed(j)) {
      pos->row(j) = ref.row(j);
    }
  }
}

// Moves an image only when the coordinates differ, so an image that
// already holds this geometry keeps its evaluation.
void moveTo(Matter &image, const AtomMatrix &pos) {
  if (image.getPositions() != pos) {
    image.setPositions(pos);
  }
}

} // namespace

int xtsciBandSurface(void *user, void *request) {
  auto *self = static_cast<XtsciBand *>(user);
  auto *req = static_cast<rgsaddle_surface_request_t *>(request);
  if (self == nullptr || req == nullptr || req->positions == nullptr ||
      req->energies == nullptr || req->gradients == nullptr) {
    return RGSADDLE_SURFACE_FAILED;
  }
  try {
    auto *neb = self->m_neb;
    const auto band = static_cast<int64_t>(neb->numImages + 2);
    const auto nAtoms = static_cast<int64_t>(neb->atoms);
    if (req->version.major != RGSADDLE_ABI_MAJOR) {
      return RGSADDLE_ABI_MISMATCH;
    }
    if (req->n_atoms != nAtoms) {
      return RGSADDLE_SHAPE;
    }
    // Rows carried by this request and the band image of row 0. The
    // session sends the whole band on its first evaluation, the interior
    // images 1 .. band - 2 afterwards, or one image with
    // RGSADDLE_REQ_ONE_IMAGE.
    int64_t rows = 0;
    int64_t first = 0;
    if ((req->flags & RGSADDLE_REQ_ONE_IMAGE) != 0) {
      if (req->n_images != band || req->image < 0 || req->image >= band) {
        return RGSADDLE_SHAPE;
      }
      rows = 1;
      first = req->image;
    } else if (req->n_images == band) {
      rows = band;
    } else if (req->n_images == band - 2) {
      rows = band - 2;
      first = 1;
    } else {
      return RGSADDLE_SHAPE;
    }
    const auto dof = static_cast<Eigen::Index>(3 * nAtoms);
    std::vector<Matter *> carried;
    carried.reserve(static_cast<size_t>(rows));
    for (int64_t r = 0; r < rows; ++r) {
      const auto i = static_cast<size_t>(first + r);
      AtomMatrix pos = AtomMatrix::Map(req->positions + r * dof, nAtoms, 3);
      if (i < self->m_fixed.size()) {
        pinRows(self->m_fixed[i], *neb->path[i], &pos);
      }
      // The fixed endpoints and a start point the band already evaluated
      // come back from the image's own cache.
      moveTo(*neb->path[i], pos);
      carried.push_back(neb->path[i].get());
    }
    // One batch for every carried image, so calculator groups or a batched
    // model take them together; a serial potential evaluates them in turn.
    eonc::evaluateTogether(*carried.front()->getPotential(), carried);
    for (int64_t r = 0; r < rows; ++r) {
      Matter &image = *carried[static_cast<size_t>(r)];
      req->energies[r] = image.getPotentialEnergy();
      const AtomMatrix &force = image.getForces();
      double *grad = req->gradients + r * dof;
      const double *src = force.data();
      for (Eigen::Index k = 0; k < dof; ++k) {
        grad[k] = -src[k];
      }
    }
    return RGSADDLE_OK;
  } catch (...) {
    return RGSADDLE_SURFACE_FAILED;
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
    m_springKs.assign(static_cast<size_t>(nImages - 1), nebOpt.spring.constant);
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

  std::vector<double> positions(static_cast<size_t>(nImages * 3 * nAtoms), 0.0);
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
  std::vector<double> positions(static_cast<size_t>(nImages * 3 * nAtoms), 0.0);
  for (int64_t i = 0; i < nImages; ++i) {
    const AtomMatrix &pos = m_neb->path[i]->getPositions();
    std::copy(pos.data(), pos.data() + pos.size(),
              positions.begin() + static_cast<std::ptrdiff_t>(i * pos.size()));
  }
  checkStatus(rgsaddle_band_set_positions(m_band, positions.data()),
              "rgsaddle_band_set_positions");
}

void XtsciBand::reset() {
#if RGSADDLE_ABI_MINOR >= 5
  // Same surface: drop the optimizer history and climbing state, keep the
  // endpoint energies and the last evaluation.
  checkStatus(rgsaddle_band_restart(m_band), "rgsaddle_band_restart");
#else
  checkStatus(rgsaddle_band_reset(m_band), "rgsaddle_band_reset");
#endif
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
  // The session holds the band of the last step. Resending it drops the
  // session's cached endpoint and start-point evaluations, so only a band
  // moved outside the session (a reparameterization) is resent.
#if RGSADDLE_ABI_MINOR >= 5
  // Unchanged rows cost nothing: the session keeps its cached values.
  syncFromPath();
#else
  if (!sessionMatchesPath()) {
    syncFromPath();
  }
#endif
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
    pinRows(m_fixed[i], *m_neb->path[i], &pos);
    // The accepted point is usually the band's last evaluation, which the
    // image still holds; the band update after the step then costs nothing.
    moveTo(*m_neb->path[i], pos);
  }
#if RGSADDLE_ABI_MINOR >= 5
  // An image whose last evaluation was elsewhere (the solver backed off to
  // a point it evaluated earlier) takes the session's evaluation of the
  // accepted band instead of a fresh force call.
  std::vector<double> energies(m_fixed.size());
  std::vector<double> gradients(m_fixed.size() * dof);
  if (rgsaddle_band_evaluation(m_band, energies.data(), gradients.data(),
                               nullptr) == RGSADDLE_OK) {
    for (size_t i = 1; i + 1 < m_fixed.size(); ++i) {
      Matter &image = *m_neb->path[i];
      if (!image.needsForceUpdate()) {
        continue;
      }
      const AtomMatrix grad = AtomMatrix::Map(
          gradients.data() + static_cast<std::ptrdiff_t>(i * dof), m_neb->atoms,
          3);
      image.setEvaluation(-grad, energies[i]);
    }
  }
#endif
  m_neb->movedAfterForceCall = true;
}

bool XtsciBand::sessionMatchesPath() const {
  const auto nImages = static_cast<int64_t>(m_neb->numImages + 2);
  const auto nAtoms = static_cast<int64_t>(m_neb->atoms);
  std::vector<double> positions(static_cast<size_t>(nImages * 3 * nAtoms));
  if (rgsaddle_band_positions(m_band, positions.data()) != RGSADDLE_OK) {
    return false;
  }
  const auto dof = static_cast<std::ptrdiff_t>(3 * nAtoms);
  for (int64_t i = 0; i < nImages; ++i) {
    const AtomMatrix &pos = m_neb->path[static_cast<size_t>(i)]->getPositions();
    if (!std::equal(pos.data(), pos.data() + dof,
                    positions.begin() + i * dof)) {
      return false;
    }
  }
  return true;
}

} // namespace eonc

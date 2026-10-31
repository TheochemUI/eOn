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

#include "ForEachImage.h"

#include <rgsaddle.h>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc {

int xtsciBandSurface(void *user, void *request);
int xtsciSolidBandSurface(void *user, void *request);

namespace {

#if RGSADDLE_ABI_MINOR >= 17
// Forward from the C callback to the friend that reads the band.
rgsaddle_status_t solidSurfaceCallback(void *user,
                                       rgsaddle_solid_surface_request_t *req) {
  auto *self = static_cast<XtsciBand *>(user);
  return static_cast<rgsaddle_status_t>(self == nullptr
                                            ? RGSADDLE_SURFACE_FAILED
                                            : xtsciSolidBandSurface(self, req));
}

// rgsaddle cells are row-major with rows the lattice vectors, as eOn's
// cell rows are; Eigen stores a Matrix3d column-major.
void cellToRowMajor(const Matrix3d &cell, double *out) {
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      out[3 * r + c] = cell(r, c);
    }
  }
}

Matrix3d cellFromRowMajor(const double *in) {
  Matrix3d cell;
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      cell(r, c) = in[3 * r + c];
    }
  }
  return cell;
}
#endif

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

#ifdef WITH_RGSADDLE_BAND_PRECON
bool bandHasFixedAtoms(const NudgedElasticBand &neb) {
  for (const auto &image : neb.path) {
    if (image->numberOfFreeAtoms() < image->numberOfAtoms()) {
      return true;
    }
  }
  return false;
}
#endif

// band_precon and band_stepper, refused here where rgsaddle 0.3.0 would
// refuse the metric. A library older than ABI minor 22 reads ODE12r as
// FIRE, so the loaded library is checked as well as the header.
void checkBandOptions(const NudgedElasticBand &neb, const neb_options_t &opt,
                      [[maybe_unused]] const rgsaddle_version_t &stamp) {
  const auto &band = opt.band;
  const bool exp = band.precon == "exp";
  const bool ode12r = band.stepper == "ode12r";
  if (!exp && band.precon != "none") {
    throw std::invalid_argument("band_precon is none or exp, not '" +
                                band.precon + "'");
  }
  if (!ode12r && band.stepper != "auto") {
    throw std::invalid_argument("band_stepper is auto or ode12r, not '" +
                                band.stepper + "'");
  }
  if (!exp && !ode12r) {
    return;
  }
#ifndef WITH_RGSADDLE_BAND_PRECON
  throw std::invalid_argument(
      "band_precon = exp and band_stepper = ode12r need eOn built against "
      "rgsaddle 0.3.0 or newer (built against ABI minor " +
      std::to_string(RGSADDLE_ABI_MINOR) + ")");
#else
  if (stamp.minor < 22) {
    throw std::runtime_error(
        "band_precon = exp and band_stepper = ode12r need rgsaddle ABI "
        "minor 22 or newer at run time (loaded minor " +
        std::to_string(stamp.minor) + ")");
  }
  if (!exp) {
    return;
  }
  if (neb.solidState()) {
    throw std::invalid_argument(
        "band_precon = exp is refused with solid_state; the metric has no "
        "pair geometry for the cell block");
  }
  if (neb.atoms < 2) {
    throw std::invalid_argument(
        "band_precon = exp needs at least two atoms per image");
  }
  if (springKind(opt) != RGSADDLE_SPRING_UNIFORM) {
    throw std::invalid_argument(
        "band_precon = exp needs a uniform spring; energy_weighted and "
        "onsager_machlup are refused");
  }
  if (bandHasFixedAtoms(neb)) {
    throw std::invalid_argument(
        "band_precon = exp is refused with fixed atoms: the rgsaddle band "
        "metric has no fixed-atom mask, so the step would move fixed atoms "
        "that eOn pins back. Free every atom, or keep band_precon = none");
  }
#endif
}

int32_t methodKind(const Parameters &params) {
#ifdef WITH_RGSADDLE_BAND_PRECON
  if (params.neb_options().band.stepper == "ode12r") {
    return RGSADDLE_METHOD_ODE12R;
  }
#endif
  const auto nebMethod = params.neb_options().opt_method;
  if (nebMethod == OptType::FIRE) {
    return RGSADDLE_METHOD_FIRE;
  }
  if (nebMethod == OptType::XTSCI &&
      params.optimizer_options().xtsci.method != "lbfgs") {
    return RGSADDLE_METHOD_FIRE;
  }
  return RGSADDLE_METHOD_LBFGS;
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
    } else if (req->n_images == 1) {
      // With RGSADDLE_BAND_NO_MIN_MODE set, only the estimate of the band
      // metric's energy scale sends one row: a point off the band near
      // image 1. A copy of image 1 takes it, so no band image moves. That
      // row comes batched, with image -1 and no flags.
      if (req->image >= 0 || req->flags != 0) {
        return RGSADDLE_SHAPE;
      }
      Matter probe(*neb->path[1]);
      moveTo(probe, AtomMatrix::Map(req->positions, nAtoms, 3));
      req->energies[0] = probe.getPotentialEnergy();
      const double *src = probe.getForces().data();
      for (int64_t k = 0; k < 3 * nAtoms; ++k) {
        req->gradients[k] = -src[k];
      }
      return RGSADDLE_OK;
    } else {
      return RGSADDLE_SHAPE;
    }
    const auto dof = static_cast<Eigen::Index>(3 * nAtoms);
    std::vector<Matter *> carried;
    std::vector<long> owners;
    carried.reserve(static_cast<size_t>(rows));
    owners.reserve(static_cast<size_t>(rows));
    for (int64_t r = 0; r < rows; ++r) {
      const auto i = static_cast<size_t>(first + r);
      carried.push_back(neb->path[i].get());
      owners.push_back(neb->imageOwner(i));
    }
    // Row r: the session's coordinates onto its image. The fixed endpoints
    // and a start point the band already evaluated come back from the
    // image's own cache.
    auto place = [&](int64_t r) {
      const auto i = static_cast<size_t>(first + r);
      AtomMatrix pos = AtomMatrix::Map(req->positions + r * dof, nAtoms, 3);
      if (i < self->m_fixed.size()) {
        pinRows(self->m_fixed[i], *neb->path[i], &pos);
      }
      moveTo(*neb->path[i], pos);
    };
    auto report = [&](int64_t r) {
      Matter &image = *carried[static_cast<size_t>(r)];
      req->energies[r] = image.getPotentialEnergy();
      const double *src = image.getForces().data();
      double *grad = req->gradients + r * dof;
      for (Eigen::Index k = 0; k < dof; ++k) {
        grad[k] = -src[k];
      }
    };
    Potential &pot = *carried.front()->getPotential();
    if (self->m_threads && rows > 1 && !pot.supportsBatchEvaluation() &&
        eonc::canComputeConcurrently(carried)) {
      // One pass on the image pool: each row is placed, evaluated and
      // read back on the thread that took it.
      eonc::forEachImage(static_cast<long>(rows), [&](long k) {
        const auto r = static_cast<int64_t>(k - 1);
        place(r);
        Matter &image = *carried[static_cast<size_t>(r)];
        if (image.needsForceUpdate()) {
          image.getForcesRaw();
        }
        report(r);
      });
    } else {
      for (int64_t r = 0; r < rows; ++r) {
        place(r);
      }
      // One batch for every carried image, so calculator groups or a
      // batched model take them together; otherwise in turn. Each image
      // carries its band identity, as in the classic update.
      eonc::evaluateTogether(pot, carried, false, owners);
      for (int64_t r = 0; r < rows; ++r) {
        report(r);
      }
    }
    return RGSADDLE_OK;
  } catch (...) {
    return RGSADDLE_SURFACE_FAILED;
  }
}

int xtsciSolidBandSurface([[maybe_unused]] void *user,
                          [[maybe_unused]] void *request) {
#if RGSADDLE_ABI_MINOR >= 17
  auto *self = static_cast<XtsciBand *>(user);
  auto *req = static_cast<rgsaddle_solid_surface_request_t *>(request);
  if (self == nullptr || req == nullptr || req->positions == nullptr ||
      req->cells == nullptr || req->energies == nullptr ||
      req->gradients == nullptr || req->stresses == nullptr) {
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
    for (int64_t r = 0; r < rows; ++r) {
      const auto i = static_cast<size_t>(first + r);
      Matter &image = *neb->path[i];
      const Matrix3d cell = cellFromRowMajor(req->cells + 9 * r);
      AtomMatrix pos = AtomMatrix::Map(req->positions + r * dof, nAtoms, 3);
      if (i < self->m_fixed.size()) {
        pinRows(self->m_fixed[i], image, &pos);
      }
      if (image.getCell() != cell) {
        image.setCell(cell);
      }
      moveTo(image, pos);
      // One image at a time: Matter keeps the stress that came with its
      // own force call, so the stress read below costs no second call.
      const AtomMatrix &force = image.getForces();
      req->energies[r] = image.getPotentialEnergy();
      double *grad = req->gradients + r * dof;
      const double *src = force.data();
      for (Eigen::Index k = 0; k < dof; ++k) {
        grad[k] = -src[k];
      }
      cellToRowMajor(image.cauchyStress(), req->stresses + 9 * r);
    }
    return RGSADDLE_OK;
  } catch (...) {
    return RGSADDLE_SURFACE_FAILED;
  }
#else
  return RGSADDLE_SURFACE_FAILED;
#endif
}

XtsciBand::XtsciBand(NudgedElasticBand &neb, const Parameters &params)
    : m_neb{&neb},
      m_maxMove{params.optimizer_options().max_move},
      m_threads{params.main_options().parallel} {
  rgsaddle_version_t stamp{};
  checkStatus(rgsaddle_abi_stamp(&stamp), "rgsaddle_abi_stamp");
  if (stamp.major != RGSADDLE_ABI_MAJOR) {
    throw std::runtime_error("incompatible rgsaddle ABI");
  }
  const auto &nebOpt = params.neb_options();
  checkBandOptions(neb, nebOpt, stamp);
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
  // eOn removes the mean force of an image only when every atom is free
  // (eonc::neb::zeroTranslation). With a fixed atom the mean stays, so the
  // session never pushes fixed atoms that the host pins back after the step.
  config.flags = 0;
  if (!neb.path.empty() &&
      neb.path[0]->numberOfFreeAtoms() < neb.path[0]->numberOfAtoms()) {
    config.flags |= RGSADDLE_BAND_KEEP_TRANSLATION;
  }
#if RGSADDLE_ABI_MINOR >= 21
  // From minor 21 a climbing session refines its climbing image alone
  // unless told not to. eOn's own controller runs OCI-NEB (ci_mmf), so
  // the session stays the plain CI-NEB, and every request keeps a band
  // shape the surface callback reads.
  config.flags |= RGSADDLE_BAND_NO_MIN_MODE;
#endif
  config.tangent = tangentKind(nebOpt);
  config.spring = springKind(nebOpt);
  config.projection = projectionKind(nebOpt);
  config.method = methodKind(params);
  config.spring_k = nebOpt.spring.om.enabled
                        ? nebOpt.spring.constant * nebOpt.spring.om.k_scale
                        : nebOpt.spring.constant;
  config.spring_ks = m_springKs.empty() ? nullptr : m_springKs.data();
  if (nebOpt.climbing_image.enabled) {
#if RGSADDLE_ABI_MINOR >= 17
    // eOn's ci_after and ci_after_rel as given; either may be 0.
    config.flags |= RGSADDLE_BAND_CI_AFTER;
    config.ci_trigger_factor = nebOpt.climbing_image.trigger_factor;
#else
    config.ci_trigger_factor = nebOpt.climbing_image.trigger_factor > 0.0
                                   ? nebOpt.climbing_image.trigger_factor
                                   : 1.0e-300;
#endif
    config.ci_trigger_force = nebOpt.climbing_image.trigger_force;
  } else {
    config.ci_trigger_factor = 0.0;
    config.ci_trigger_force = 0.0;
  }
  m_solid = neb.solidState();
  if (!m_solid && !neb.path.empty() && neb.path[0]->getPeriodic()) {
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
  if (m_solid) {
#if RGSADDLE_ABI_MINOR >= 17
    // Each image keeps its cell; prepareSolidState already put every cell
    // in lower-triangular form, which rgsaddle requires.
    std::vector<double> cells(static_cast<size_t>(nImages * 9), 0.0);
    for (int64_t i = 0; i < nImages; ++i) {
      cellToRowMajor(neb.path[i]->getCell(), cells.data() + 9 * i);
    }
    const auto &solid = nebOpt.solid_state;
    m_band =
        rgsaddle_band_create_solid(&config, nImages, nAtoms, positions.data(),
                                   cells.data(), solid.pressure, solid.weight);
    if (m_band == nullptr) {
      throw std::runtime_error("rgsaddle_band_create_solid failed");
    }
    return;
#else
    throw std::runtime_error(
        "solid_state NEB on the rgsaddle band needs rgsaddle ABI minor 17 "
        "or newer (built against minor " +
        std::to_string(RGSADDLE_ABI_MINOR) +
        "); the older band keeps every cell fixed. Use a newer rgsaddle or "
        "an eOn optimizer for opt_method");
#endif
  }
  m_band = rgsaddle_band_create(&config, nImages, nAtoms, positions.data());
  if (m_band == nullptr) {
    throw std::runtime_error("rgsaddle_band_create failed");
  }
#ifdef WITH_RGSADDLE_BAND_PRECON
  if (nebOpt.band.precon == "exp") {
    // a, c_stab and r_cut_factor at 0 take rgsaddle's 3, 0.1 and 2.
    const int rc = rgsaddle_band_set_exp_precon(m_band, 0.0, 0.0,
                                                nebOpt.band.precon_mu, 0.0);
    if (rc != RGSADDLE_OK) {
      rgsaddle_band_free(m_band);
      m_band = nullptr;
      checkStatus(rc, "rgsaddle_band_set_exp_precon");
    }
  }
#endif
}

XtsciBand::~XtsciBand() { rgsaddle_band_free(m_band); }

template <typename Work> void XtsciBand::forImages(int64_t n, Work &&work) {
  if (m_threads && n > 1) {
    eonc::forEachImage(static_cast<long>(n),
                       [&](long k) { work(static_cast<int64_t>(k - 1)); });
    return;
  }
  for (int64_t i = 0; i < n; ++i) {
    work(i);
  }
}

void XtsciBand::syncFromPath() {
  const auto nImages = static_cast<int64_t>(m_neb->numImages + 2);
  const auto dof = static_cast<std::ptrdiff_t>(3 * m_neb->atoms);
  std::vector<double> positions(static_cast<size_t>(nImages * dof), 0.0);
  forImages(nImages, [&](int64_t i) {
    const AtomMatrix &pos = m_neb->path[static_cast<size_t>(i)]->getPositions();
    std::copy(pos.data(), pos.data() + dof, positions.begin() + i * dof);
  });
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
  if (m_solid) {
    stepSolid();
    return;
  }
  const auto nImages = static_cast<int64_t>(m_neb->numImages + 2);
  m_fixed.resize(static_cast<size_t>(nImages));
  forImages(nImages, [&](int64_t i) {
    m_fixed[static_cast<size_t>(i)] =
        m_neb->path[static_cast<size_t>(i)]->getPositions();
  });
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
  const int rc = xts_band_step(m_band, surfaceCallback, this, &report);
#if RGSADDLE_ABI_MINOR >= 17
  // A stalled step leaves the band as it was. The NEB loop then runs into
  // max_iterations, as it does when a host optimizer stops making progress.
  if (rc == RGSADDLE_STALLED) {
    return;
  }
#endif
  checkStatus(rc, "xts_band_step");
  const auto dof = static_cast<std::ptrdiff_t>(m_fixed[0].size());
  std::vector<double> positions(m_fixed.size() * static_cast<size_t>(dof), 0.0);
  checkStatus(rgsaddle_band_positions(m_band, positions.data()),
              "rgsaddle_band_positions");
  // The session reports -1 until an image climbs; the host then holds no
  // climbing image either, so its convergence test reads the whole band.
  m_neb->climbingImage = report.ci_index > 0 ? report.ci_index : 0;
  const int64_t interior = nImages - 2;
#if RGSADDLE_ABI_MINOR >= 5
  // The session's evaluation of the accepted band and its projected force.
  // An image whose last evaluation was elsewhere (the solver backed off to
  // a point it evaluated earlier) takes that evaluation instead of a fresh
  // force call, and the host convergence check reads the session's force.
  std::vector<double> energies(m_fixed.size());
  std::vector<double> gradients(m_fixed.size() * static_cast<size_t>(dof));
  const bool evaluated =
      rgsaddle_band_evaluation(m_band, energies.data(), gradients.data(),
                               nullptr) == RGSADDLE_OK;
  std::vector<double> projected(static_cast<size_t>(interior * dof));
  const bool haveProjection =
      rgsaddle_band_evaluation(m_band, nullptr, nullptr, projected.data()) ==
      RGSADDLE_OK;
#endif
  // Endpoints are fixed by the session. Each interior image takes the
  // accepted point; interior fixed atoms stay put. The accepted point is
  // usually the band's last evaluation, which the image still holds.
  forImages(interior, [&](int64_t k) {
    const auto i = static_cast<size_t>(k + 1);
    Matter &image = *m_neb->path[i];
    AtomMatrix pos =
        AtomMatrix::Map(positions.data() + static_cast<std::ptrdiff_t>(i) * dof,
                        m_neb->atoms, 3);
    pinRows(m_fixed[i], image, &pos);
    moveTo(image, pos);
#if RGSADDLE_ABI_MINOR >= 5
    if (evaluated && image.needsForceUpdate()) {
      const AtomMatrix grad = AtomMatrix::Map(
          gradients.data() + static_cast<std::ptrdiff_t>(i) * dof, m_neb->atoms,
          3);
      image.setEvaluation(-grad, energies[i]);
    }
    if (haveProjection) {
      *m_neb->projectedForce[i] = AtomMatrix::Map(
          projected.data() + static_cast<std::ptrdiff_t>(k) * dof, m_neb->atoms,
          3);
    }
#endif
  });
#if RGSADDLE_ABI_MINOR >= 5
  if (haveProjection) {
    // The session already projected the force. The host convergence
    // check reads that force and does not rebuild it, so the highest
    // image (sp.con), the climbing image and the tangents (neb.dat) are
    // refreshed here from the accepted band.
    m_neb->adoptSessionProjection(report.ci_index);
    m_neb->movedAfterForceCall = false;
    return;
  }
#endif
  m_neb->movedAfterForceCall = true;
}

void XtsciBand::stepSolid() {
#if RGSADDLE_ABI_MINOR >= 17
  const auto nImages = static_cast<int64_t>(m_neb->numImages + 2);
  const auto nAtoms = static_cast<int64_t>(m_neb->atoms);
  m_fixed.resize(static_cast<size_t>(nImages));
  for (int64_t i = 0; i < nImages; ++i) {
    m_fixed[static_cast<size_t>(i)] = m_neb->path[i]->getPositions();
  }
  rgsaddle_report_t report{};
  const int rc =
      rgsaddle_band_step_solid(m_band, solidSurfaceCallback, this, &report);
  if (rc == RGSADDLE_STALLED) {
    return;
  }
  checkStatus(rc, "rgsaddle_band_step_solid");
  const auto dof = static_cast<std::ptrdiff_t>(3 * nAtoms);
  std::vector<double> positions(static_cast<size_t>(nImages * dof));
  std::vector<double> cells(static_cast<size_t>(nImages * 9));
  checkStatus(rgsaddle_band_positions(m_band, positions.data()),
              "rgsaddle_band_positions");
  checkStatus(rgsaddle_band_cells(m_band, cells.data()), "rgsaddle_band_cells");
  for (int64_t i = 1; i + 1 < nImages; ++i) {
    Matter &image = *m_neb->path[static_cast<size_t>(i)];
    AtomMatrix pos =
        AtomMatrix::Map(positions.data() + i * dof, m_neb->atoms, 3);
    pinRows(m_fixed[static_cast<size_t>(i)], image, &pos);
    const Matrix3d cell = cellFromRowMajor(cells.data() + 9 * i);
    if (image.getCell() != cell) {
      image.setCell(cell);
    }
    moveTo(image, pos);
  }
  if (report.ci_index >= 0) {
    m_neb->climbingImage = report.ci_index;
  }
  // eOn's convergence force reads its own joint projection of the moved
  // band (projectedForce and projectedCellForce); rebuild it there.
  m_neb->movedAfterForceCall = true;
#endif
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

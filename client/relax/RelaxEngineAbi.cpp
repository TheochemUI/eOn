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
#include "eon/relax/eon_relax_engine.h"

#include "relax/HostSurfacePotential.h"
#include "eon/Eigen.h"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/NudgedElasticBand.h"
#include "eon/Optimizer.h"
#include "eon/Parameters.h"
#include "eon_relax_engine.capnp.h"
#include "version.h"

#include <capnp/serialize.h>
#include <kj/array.h>
#include <kj/common.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

// Library TU. The test alias header is not included here.
using eonc::Matter;
using eonc::MinModeSaddleSearch;
using eonc::NEBObjectiveFunction;
using eonc::NudgedElasticBand;
using eonc::Optimizer;
using eonc::Parameters;
using eonc::ParametersLoadAccess;

extern "C" uint64_t eon_relax_version_hash(void);

namespace {

uint64_t fnv1a64(const char *s) {
  uint64_t h = 14695981039346656037ULL;
  if (!s) {
    return h;
  }
  for (const unsigned char *p = reinterpret_cast<const unsigned char *>(s); *p;
       ++p) {
    h ^= static_cast<uint64_t>(*p);
    h *= 1099511628211ULL;
  }
  return h;
}

void set_err(char *errbuf, size_t errlen, const char *msg) {
  if (!errbuf || errlen == 0) {
    return;
  }
  std::snprintf(errbuf, errlen, "%s", msg ? msg : "");
}

bool known_kind(eon_relax_kind_t kind) { return EON_RELAX_KIND_IS_KNOWN(kind); }

int neb_status(NudgedElasticBand::NEBStatus st) {
  return static_cast<int>(st);
}

void fill_matter(Matter &m, const eon_relax_band_t *band, long image,
                 uint64_t epoch) {
  m.resize(band->n_atoms);
  Eigen::VectorXi z(band->n_atoms);
  for (long a = 0; a < band->n_atoms; ++a) {
    z(a) = band->atomic_nrs[a];
  }
  m.setAtomicNrs(z);
  Eigen::VectorXd masses = Eigen::VectorXd::Ones(band->n_atoms);
  if (band->version.minor >= 1 && band->masses) {
    for (long a = 0; a < band->n_atoms; ++a) {
      masses(a) = band->masses[a];
    }
  }
  m.setMasses(masses);
  AtomMatrix pos(band->n_atoms, 3);
  const double *src = band->positions + image * 3 * band->n_atoms;
  // Match Potential::force / Matter storage (AtomMatrix row-major xyz).
  Eigen::Map<const AtomMatrix> mapped(src, band->n_atoms, 3);
  pos = mapped;
  const double *b = band->boxes + image * 9;
  const Matrix3d cell = Eigen::Map<const Matrix3d>(b);
  // Cell first. PBC stays off so dest optimizer writes
  // (setPositionsFree) do not wrap a host-frame x=-0.667 to 9.333
  // and hop the climbing image onto a wrapped shoulder.
  m.setPeriodic(false);
  m.setCell(cell);
  m.setPositions(pos);
  if (band->is_fixed) {
    for (long a = 0; a < band->n_atoms; ++a) {
      m.setFixed(a, band->is_fixed[a] ? 1 : 0);
    }
  } else {
    for (long a = 0; a < band->n_atoms; ++a) {
      m.setFixed(a, 0);
    }
  }
  m.setSurfaceEpoch(epoch);
}

void write_image(double *dst, const Matter &m, long n_atoms) {
  const AtomMatrix &pos = m.getPositions();
  Eigen::Map<AtomMatrix> mapped(dst, n_atoms, 3);
  mapped = pos;
}

int masses_ok(const eon_relax_band_t *band) {
  if (!band || band->version.minor < 1 || !band->masses) {
    return EON_RELAX_OK;
  }
  for (int64_t a = 0; a < band->n_atoms; ++a) {
    const double m = band->masses[a];
    if (!(m > 0.0) || !std::isfinite(m)) {
      return EON_RELAX_INVALID_PARAMETER;
    }
  }
  return EON_RELAX_OK;
}

void stamp_outcome(eon_relax_outcome_t *out, eon_relax_kind_t kind,
                   uint64_t epoch) {
  std::memset(out, 0, sizeof(*out));
  out->version.major = EON_RELAX_ABI_MAJOR;
  out->version.minor = EON_RELAX_ABI_MINOR;
  out->kind = static_cast<int32_t>(kind);
  out->version_hash = eon_relax_version_hash();
  out->surface_epoch = epoch;
}

void dirty_keep_frame(Matter &m) {
  const bool pbc = m.getPeriodic();
  m.setPeriodic(false);
  m.setPositions(m.getPositions());
  m.setPeriodic(pbc);
}

void dirty_endpoints(NudgedElasticBand &neb) {
  dirty_keep_frame(*neb.path[0]);
  dirty_keep_frame(*neb.path[static_cast<size_t>(neb.numImages + 1)]);
}

int apply_relax_params(Parameters &params, eon_relax_kind_t *kind,
                       uint64_t *epoch, const void *config, size_t config_len,
                       char *errbuf, size_t errlen) {
  if (config_len % sizeof(::capnp::word) != 0) {
    set_err(errbuf, errlen, "-13 RelaxEngineParams capnp root");
    return EON_RELAX_CAPNP_ROOT;
  }
  try {
    const auto nwords = config_len / sizeof(::capnp::word);
    std::vector<::capnp::word> words(nwords);
    std::memcpy(words.data(), config, config_len);
    ::kj::ArrayPtr<const ::capnp::word> view(words.data(), nwords);
    ::capnp::FlatArrayMessageReader reader(view);
    auto root = reader.getRoot<eonc::params_ssot::RelaxEngineParams>();
    const auto kind_txt = root.getKind();
    if (kind_txt == "neb") {
      *kind = EON_RELAX_KIND_NEB;
    } else if (kind_txt == "saddle") {
      *kind = EON_RELAX_KIND_SADDLE;
    } else {
      set_err(errbuf, errlen, "-11 unknown RelaxEngineParams.kind");
      return EON_RELAX_UNKNOWN_KIND;
    }
    auto neb = root.getNeb();
    auto &neb_o = ParametersLoadAccess::neb_options(params);
    neb_o.image_count = static_cast<long>(neb.getImageCount());
    neb_o.max_iterations = static_cast<long>(neb.getMaxIterations());
    neb_o.force_tolerance = neb.getForceTolerance();
    neb_o.climbing_image.enabled = neb.getClimbingImage();
    neb_o.endpoints.minimize = neb.getMinimizeEndpoints();
    auto saddle = root.getSaddle();
    auto &saddle_o = ParametersLoadAccess::saddle_search_options(params);
    saddle_o.max_iterations = static_cast<long>(saddle.getMaxIterations());
    saddle_o.converged_force = saddle.getConvergedForce();
    *epoch = root.getSurfaceEpoch();
    ParametersLoadAccess::main_options(params).randomSeed =
        static_cast<long>(root.getRandomSeed());
    ParametersLoadAccess::gp_surrogate_options(params).uncertainty =
        root.getUncertainty();
    return EON_RELAX_OK;
  } catch (const std::exception &) {
    set_err(errbuf, errlen, "-13 RelaxEngineParams capnp root");
    return EON_RELAX_CAPNP_ROOT;
  } catch (...) {
    set_err(errbuf, errlen, "-13 RelaxEngineParams capnp root");
    return EON_RELAX_CAPNP_ROOT;
  }
}

} // namespace

struct EonRelaxEngine {
  Parameters params;
  eon_relax_kind_t kind{EON_RELAX_KIND_NEB};
  uint64_t epoch{0};
  int last_rc{EON_RELAX_OK};
  bool unusable{false};

  // Stepper state (eon_relax_step): band, objective, and optimizer
  // persist across calls so quasi-Newton history accumulates on one
  // surface epoch; eon_relax_reset drops all of it.
  std::shared_ptr<eonc::HostSurfacePotential> step_pot;
  std::unique_ptr<NudgedElasticBand> step_neb;
  std::shared_ptr<NEBObjectiveFunction> step_objf;
  std::unique_ptr<Optimizer> step_opt;
  eon_relax_surface_fn step_surface{nullptr};
  void *step_user{nullptr};
  double step_baseline{0.0};
  long step_iteration{0};
};

static int stamp_rc(EonRelaxEngine *eng, int rc) {
  if (eng) {
    eng->last_rc = rc;
  }
  return rc;
}

static int fail_out(EonRelaxEngine *eng, eon_relax_outcome_t *out, int rc,
                    eon_relax_kind_t kind, uint64_t epoch) {
  if (out) {
    stamp_outcome(out, kind, epoch);
    out->status = -1;
  }
  return stamp_rc(eng, rc);
}

static void clear_stepper(EonRelaxEngine *eng) {
  eng->step_opt.reset();
  eng->step_objf.reset();
  eng->step_neb.reset();
  eng->step_pot.reset();
  eng->step_surface = nullptr;
  eng->step_user = nullptr;
  eng->step_baseline = 0.0;
  eng->step_iteration = 0;
}

extern "C" {

int eon_relax_abi_version(void) {
  return static_cast<int>((EON_RELAX_ABI_MAJOR << 16) | EON_RELAX_ABI_MINOR);
}

int eon_relax_available(void) { return 1; }

int eon_relax_abi_stamp(eon_relax_version_t *out) {
  if (!out) {
    return EON_RELAX_INVALID_PARAMETER;
  }
  out->major = EON_RELAX_ABI_MAJOR;
  out->minor = EON_RELAX_ABI_MINOR;
  return EON_RELAX_OK;
}

const char *eon_relax_version_hash_str(void) {
  static const std::string id = VERSION + "+git." + GIT_HASH_FULL + "+feat." +
                                std::to_string(fnv1a64(FEATURES_STRING.c_str()));
  return id.c_str();
}

uint64_t eon_relax_version_hash(void) {
  return fnv1a64(eon_relax_version_hash_str());
}

EonRelaxEngine *eon_relax_create(const void *config, size_t config_len,
                                 char *errbuf, size_t errlen) {
  if (config == nullptr && config_len == 0) {
    auto *eng = new (std::nothrow) EonRelaxEngine();
    if (!eng) {
      set_err(errbuf, errlen, "-19 out of memory");
      return nullptr;
    }
    eng->kind = EON_RELAX_KIND_NEB;
    return eng;
  }
  if (config == nullptr || config_len == 0) {
    set_err(errbuf, errlen, "-13 RelaxEngineParams capnp root");
    return nullptr;
  }
  auto *eng = new (std::nothrow) EonRelaxEngine();
  if (!eng) {
    set_err(errbuf, errlen, "-19 out of memory");
    return nullptr;
  }
  const int prc = apply_relax_params(eng->params, &eng->kind, &eng->epoch,
                                     config, config_len, errbuf, errlen);
  if (prc != EON_RELAX_OK) {
    delete eng;
    return nullptr;
  }
  return eng;
}

int eon_relax_set_surface_epoch(EonRelaxEngine *eng, uint64_t epoch) {
  if (!eng) {
    return EON_RELAX_NULL_ENGINE;
  }
  if (eng->unusable) {
    return stamp_rc(eng, EON_RELAX_SURFACE_FATAL);
  }
  eng->epoch = epoch;
  if (eng->step_pot) {
    eng->step_pot->setEpoch(epoch);
  }
  eng->last_rc = EON_RELAX_OK;
  return EON_RELAX_OK;
}

int eon_relax_last_error(const EonRelaxEngine *eng) {
  if (!eng) {
    return EON_RELAX_NULL_ENGINE;
  }
  return eng->last_rc;
}

int eon_relax_run(EonRelaxEngine *eng, eon_relax_band_t *band,
                  eon_relax_surface_fn surface, void *user,
                  eon_relax_outcome_t *out) {
  if (!eng) {
    return EON_RELAX_NULL_ENGINE;
  }
  if (eng->unusable) {
    return fail_out(eng, out, EON_RELAX_SURFACE_FATAL, eng->kind, eng->epoch);
  }
  if (!band) {
    return fail_out(eng, out, EON_RELAX_NULL_BAND, eng->kind, eng->epoch);
  }
  if (!surface) {
    return fail_out(eng, out, EON_RELAX_NULL_SURFACE, eng->kind, eng->epoch);
  }
  if (!out) {
    return stamp_rc(eng, EON_RELAX_NULL_OUTCOME);
  }
  clear_stepper(eng);
  if (band->n_atoms <= 0) {
    return fail_out(eng, out, EON_RELAX_NATOMS, eng->kind, eng->epoch);
  }
  if (!band->positions) {
    return fail_out(eng, out, EON_RELAX_NULL_POSITIONS, eng->kind, eng->epoch);
  }
  if (!band->atomic_nrs) {
    return fail_out(eng, out, EON_RELAX_NULL_ATOMIC_NRS, eng->kind, eng->epoch);
  }
  if (!band->boxes) {
    return fail_out(eng, out, EON_RELAX_NULL_BOXES, eng->kind, eng->epoch);
  }
  if (band->version.major != EON_RELAX_ABI_MAJOR) {
    return fail_out(eng, out, EON_RELAX_ABI_MISMATCH, eng->kind, eng->epoch);
  }
  if (!known_kind(eng->kind)) {
    return fail_out(eng, out, EON_RELAX_UNKNOWN_KIND, eng->kind, eng->epoch);
  }
  const int mass_rc = masses_ok(band);
  if (mass_rc != EON_RELAX_OK) {
    return fail_out(eng, out, mass_rc, eng->kind, eng->epoch);
  }
  stamp_outcome(out, eng->kind, eng->epoch);

  auto pot = std::make_shared<eonc::HostSurfacePotential>(surface, user,
                                                          eng->epoch);
  try {
    if (eng->kind == EON_RELAX_KIND_NEB) {
      if (band->n_images < 3) {
        out->status = -1;
        return stamp_rc(eng, EON_RELAX_BAND_TOO_SHORT);
      }
      const int64_t want =
          static_cast<int64_t>(eng->params.neb_options().image_count) + 2;
      if (band->n_images != want) {
        out->status = -1;
        return stamp_rc(eng, EON_RELAX_BAND_SIZE);
      }
      std::vector<Matter> path;
      path.reserve(static_cast<size_t>(band->n_images));
      for (int64_t i = 0; i < band->n_images; ++i) {
        Matter img(pot, eng->params);
        fill_matter(img, band, static_cast<long>(i), eng->epoch);
        path.push_back(std::move(img));
      }
      auto neb = std::make_unique<NudgedElasticBand>(std::move(path),
                                                     eng->params, pot);
      pot->bindPath(neb->path);
      dirty_endpoints(*neb);
      const auto st = neb->compute();
      if (st == NudgedElasticBand::NEBStatus::INIT ||
          st == NudgedElasticBand::NEBStatus::RUNNING) {
        out->status = neb_status(st);
        return stamp_rc(eng, EON_RELAX_UNKNOWN_STATUS);
      }
      out->status = neb_status(st);
      NEBObjectiveFunction ranked(neb.get(), eng->params);
      if (ranked.isUncertain()) {
        out->status = EON_RELAX_NEB_MAX_UNCERTAINTY;
      }
      out->iterations = neb->lastIteration();
      out->climbing_image = neb->climbingImage;
      out->max_force = neb->convergenceForce();
      out->surface_epoch = pot->surfaceEpoch();
      for (int64_t i = 1; i + 1 < band->n_images; ++i) {
        write_image(band->positions + i * 3 * band->n_atoms,
                    *neb->path[static_cast<size_t>(i)],
                    static_cast<long>(band->n_atoms));
      }
      return stamp_rc(eng, EON_RELAX_OK);
    }

    if (band->n_images != 1) {
      out->status = -1;
      return stamp_rc(eng, EON_RELAX_SADDLE_NIMAGES);
    }
    if (!band->mode) {
      out->status = -1;
      return stamp_rc(eng, EON_RELAX_NULL_MODE);
    }
    auto matter = std::make_shared<Matter>(pot, eng->params);
    fill_matter(*matter, band, 0, eng->epoch);
    AtomMatrix mode(band->n_atoms, 3);
    Eigen::Map<const AtomMatrix> mmap(band->mode, band->n_atoms, 3);
    mode = mmap;
    const double reactant_e = matter->getPotentialEnergy();
    MinModeSaddleSearch search(matter, mode, reactant_e, eng->params, pot);
    const int sst = search.run();
    write_image(band->positions, *matter,
                static_cast<long>(band->n_atoms));
    out->iterations = search.getIterationCount();
    out->surface_epoch = pot->surfaceEpoch();
    out->status = sst;
    return stamp_rc(eng, EON_RELAX_OK);
  } catch (const eonc::SurfaceRecoverable &rec) {
    out->status = -1;
    return stamp_rc(eng, rec.rc);
  } catch (const std::exception &) {
    out->status = -1;
    eng->unusable = true;
    return stamp_rc(eng, EON_RELAX_SURFACE_FATAL);
  }
}

int eon_relax_step(EonRelaxEngine *eng, eon_relax_band_t *band,
                   eon_relax_surface_fn surface, void *user,
                   eon_relax_outcome_t *out) {
  if (!eng) {
    return EON_RELAX_NULL_ENGINE;
  }
  if (eng->unusable) {
    return fail_out(eng, out, EON_RELAX_SURFACE_FATAL, eng->kind, eng->epoch);
  }
  if (!band) {
    return fail_out(eng, out, EON_RELAX_NULL_BAND, eng->kind, eng->epoch);
  }
  if (!surface) {
    return fail_out(eng, out, EON_RELAX_NULL_SURFACE, eng->kind, eng->epoch);
  }
  if (!out) {
    return stamp_rc(eng, EON_RELAX_NULL_OUTCOME);
  }
  if (eng->kind == EON_RELAX_KIND_SADDLE) {
    return fail_out(eng, out, EON_RELAX_INVALID_PARAMETER, eng->kind,
                    eng->epoch);
  }
  if (eng->kind != EON_RELAX_KIND_NEB) {
    return fail_out(eng, out, EON_RELAX_UNKNOWN_KIND, eng->kind, eng->epoch);
  }
  if (band->n_atoms <= 0) {
    return fail_out(eng, out, EON_RELAX_NATOMS, EON_RELAX_KIND_NEB, eng->epoch);
  }
  if (!band->positions) {
    return fail_out(eng, out, EON_RELAX_NULL_POSITIONS, EON_RELAX_KIND_NEB,
                    eng->epoch);
  }
  if (!band->atomic_nrs) {
    return fail_out(eng, out, EON_RELAX_NULL_ATOMIC_NRS, EON_RELAX_KIND_NEB,
                    eng->epoch);
  }
  if (!band->boxes) {
    return fail_out(eng, out, EON_RELAX_NULL_BOXES, EON_RELAX_KIND_NEB,
                    eng->epoch);
  }
  if (band->version.major != EON_RELAX_ABI_MAJOR) {
    return fail_out(eng, out, EON_RELAX_ABI_MISMATCH, EON_RELAX_KIND_NEB,
                    eng->epoch);
  }
  const int mass_rc = masses_ok(band);
  if (mass_rc != EON_RELAX_OK) {
    return fail_out(eng, out, mass_rc, EON_RELAX_KIND_NEB, eng->epoch);
  }
  stamp_outcome(out, EON_RELAX_KIND_NEB, eng->epoch);
  if (eng->step_neb &&
      (band->n_atoms != eng->step_neb->atoms ||
       band->n_images != eng->step_neb->numImages + 2)) {
    out->status = -1;
    return stamp_rc(eng, EON_RELAX_INVALID_PARAMETER);
  }
  if (band->n_images < 3) {
    out->status = -1;
    return stamp_rc(eng, EON_RELAX_BAND_TOO_SHORT);
  }
  const int64_t want =
      static_cast<int64_t>(eng->params.neb_options().image_count) + 2;
  if (band->n_images != want) {
    out->status = -1;
    return stamp_rc(eng, EON_RELAX_BAND_SIZE);
  }
  if (eng->step_neb &&
      (eng->step_surface != surface || eng->step_user != user)) {
    out->status = -1;
    return stamp_rc(eng, EON_RELAX_INVALID_PARAMETER);
  }

  try {
    if (!eng->step_neb) {
      eng->step_pot = std::make_shared<eonc::HostSurfacePotential>(
          surface, user, eng->epoch);
      std::vector<Matter> path;
      path.reserve(static_cast<size_t>(band->n_images));
      for (int64_t i = 0; i < band->n_images; ++i) {
        Matter img(eng->step_pot, eng->params);
        fill_matter(img, band, static_cast<long>(i), eng->epoch);
        path.push_back(std::move(img));
      }
      eng->step_neb = std::make_unique<NudgedElasticBand>(std::move(path),
                                                          eng->params,
                                                          eng->step_pot);
      eng->step_pot->bindPath(eng->step_neb->path);
      dirty_endpoints(*eng->step_neb);
      eng->step_neb->E_ref =
          std::min(eng->step_neb->path[0]->getPotentialEnergy(),
                   eng->step_neb
                       ->path[eng->step_neb->numImages + 1]
                       ->getPotentialEnergy());
      eng->step_neb->updateForces();
      eng->step_objf = std::make_shared<NEBObjectiveFunction>(
          eng->step_neb.get(), eng->params);
      eng->step_opt = eonc::helpers::create::mkOptim(
          eng->step_objf, eng->params.neb_options().opt_method, eng->params);
      eng->step_surface = surface;
      eng->step_user = user;
      eng->step_baseline = eng->step_neb->convergenceForce();
      eng->step_iteration = 0;
    } else {
      // The host may have moved images between steps (acquisition,
      // reparameterization); resync and mark the band dirty.
      for (int64_t i = 0; i < band->n_images; ++i) {
        fill_matter(*eng->step_neb->path[static_cast<size_t>(i)], band,
                    static_cast<long>(i), eng->epoch);
      }
      eng->step_neb->movedAfterForceCall = true;
      eng->step_neb->updateForces();
    }

    NudgedElasticBand &neb = *eng->step_neb;
    const double convForce = neb.convergenceForce();
    const auto &ci_opt = eng->params.neb_options().climbing_image;
    const bool ci_active =
        ci_opt.enabled &&
        (convForce < eng->step_baseline * ci_opt.trigger_factor ||
         convForce < ci_opt.trigger_force);
    neb.setCIEnabled(ci_active);

    const double force_tol = eng->params.neb_options().force_tolerance;
    const long max_iter = eng->params.neb_options().max_iterations;
    int status = EON_RELAX_NEB_RUNNING;
    if (eng->step_objf && eng->step_objf->isUncertain()) {
      status = EON_RELAX_NEB_MAX_UNCERTAINTY;
    } else if (convForce <= force_tol) {
      status = EON_RELAX_NEB_GOOD;
    } else if (eng->step_iteration >= max_iter) {
      status = EON_RELAX_NEB_BAD_MAX_ITERATIONS;
    } else {
      eng->step_opt->step(eng->params.optimizer_options().max_move);
      ++eng->step_iteration;
      neb.updateForces();
      if (eng->step_objf && eng->step_objf->isUncertain()) {
        status = EON_RELAX_NEB_MAX_UNCERTAINTY;
      } else if (neb.convergenceForce() <= force_tol) {
        status = EON_RELAX_NEB_GOOD;
      } else if (eng->step_iteration >= max_iter) {
        status = EON_RELAX_NEB_BAD_MAX_ITERATIONS;
      }
    }

    // Interior images only: eOn never moves the endpoints, and Matter
    // wraps positions into the cell, so writing them back would hand
    // the caller a different frame than the one it supplied.
    for (int64_t i = 1; i + 1 < band->n_images; ++i) {
      write_image(band->positions + i * 3 * band->n_atoms,
                  *neb.path[static_cast<size_t>(i)],
                  static_cast<long>(band->n_atoms));
    }
    out->status = status;
    out->iterations = eng->step_iteration;
    out->climbing_image = neb.climbingImage;
    out->max_force = neb.convergenceForce();
    out->surface_epoch = eng->step_pot->surfaceEpoch();
    return stamp_rc(eng, EON_RELAX_OK);
  } catch (const eonc::SurfaceRecoverable &rec) {
    out->status = -1;
    clear_stepper(eng);
    return stamp_rc(eng, rec.rc);
  } catch (const std::exception &) {
    out->status = -1;
    eng->unusable = true;
    return stamp_rc(eng, EON_RELAX_SURFACE_FATAL);
  }
}

int eon_relax_reset(EonRelaxEngine *eng) {
  if (!eng) {
    return EON_RELAX_NULL_ENGINE;
  }
  if (eng->unusable) {
    return stamp_rc(eng, EON_RELAX_SURFACE_FATAL);
  }
  clear_stepper(eng);
  eng->last_rc = EON_RELAX_OK;
  return EON_RELAX_OK;
}

void eon_relax_destroy(EonRelaxEngine *eng) {
  delete eng;
}

const char *eon_relax_status_name(eon_relax_kind_t kind, int status) {
  if (kind == EON_RELAX_KIND_NEB) {
    switch (status) {
    case EON_RELAX_NEB_GOOD:
      return "GOOD";
    case EON_RELAX_NEB_INIT:
      return "INIT";
    case EON_RELAX_NEB_BAD_MAX_ITERATIONS:
      return "BAD_MAX_ITERATIONS";
    case EON_RELAX_NEB_RUNNING:
      return "RUNNING";
    case EON_RELAX_NEB_MAX_UNCERTAINTY:
      return "MAX_UNCERTAINTY";
    default:
      return nullptr;
    }
  }
  if (kind == EON_RELAX_KIND_SADDLE) {
    const auto msg = MinModeSaddleSearch::statusMessage(status);
    if (msg == "Unknown status") {
      return nullptr;
    }
    return msg.data();
  }
  return nullptr;
}

} // extern "C"

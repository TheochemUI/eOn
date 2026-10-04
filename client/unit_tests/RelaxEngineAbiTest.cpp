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
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/ConFileIO.h"
#include "eon/DynLib.h"
#include "eon/Eigen.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/NudgedElasticBand.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eon/relax/eon_relax_engine.h"
#include "relax/HostSurfacePotential.h"
#include "eon_relax_engine.capnp.h"

#include <capnp/message.h>
#include <capnp/serialize.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

struct SurfaceCtx {
  std::shared_ptr<Potential> pot;
  double variance{0.0};
  uint64_t epoch{0};
  long calls{0};
  int64_t product_id{-1};
  uint64_t seen_epoch{0};
  long recover_on_call{0};
};

static int surface_forward(void *user, eon_relax_surface_request_t *req) {
  auto *ctx = static_cast<SurfaceCtx *>(user);
  ++ctx->calls;
  if (req->epoch_out) {
    ctx->seen_epoch = *req->epoch_out;
  }
  if (ctx->recover_on_call > 0 && ctx->calls == ctx->recover_on_call) {
    return 5;
  }
  const long n_images = static_cast<long>(req->n_images);
  const long n_atoms = static_cast<long>(req->n_atoms);
  const double *positions = req->positions;
  const int32_t *atomic_nrs_i32 = req->atomic_nrs;
  std::vector<int> nrs(atomic_nrs_i32, atomic_nrs_i32 + n_atoms);
  const int *atomic_nrs = nrs.data();
  const double *boxes = req->boxes;
  double *energies = req->energies;
  double *forces = req->forces;
  double *variances = req->variances;
  uint64_t *epoch_out = req->epoch_out;
  for (long s = 0; s < n_images; ++s) {
    if (ctx->pot) {
      double var = 0.0;
      ctx->pot->force(n_atoms, positions + s * 3 * n_atoms, atomic_nrs,
                      forces + s * 3 * n_atoms, &energies[s], &var,
                      boxes + s * 9);
    } else {
      const double *p = positions + s * 3 * n_atoms;
      double *f = forces + s * 3 * n_atoms;
      double e = 0.0;
      for (long a = 0; a < n_atoms; ++a) {
        const double x = p[3 * a];
        const double y = p[3 * a + 1];
        const double z = p[3 * a + 2];
        e += 0.5 * (x * x + y * y + z * z);
        f[3 * a] = -x;
        f[3 * a + 1] = -y;
        f[3 * a + 2] = -z;
      }
      energies[s] = e;
    }
    if (variances) {
      variances[s] = ctx->variance;
    }
  }
  if (epoch_out && ctx->epoch != 0) {
    *epoch_out = ctx->epoch;
  }
  return 0;
}

static void pack_harmonic_band(std::vector<double> &pos,
                               std::vector<double> &boxes, std::vector<int32_t> &z,
                               long n_images) {
  const long n_atoms = 1;
  pos.assign(static_cast<size_t>(n_images * 3 * n_atoms), 0.0);
  boxes.assign(static_cast<size_t>(n_images * 9), 0.0);
  z.assign(1, 1);
  for (long im = 0; im < n_images; ++im) {
    const double t =
        static_cast<double>(im) / static_cast<double>(n_images - 1);
    pos[static_cast<size_t>(im * 3)] = 0.3 - 0.6 * t;
    double *b = boxes.data() + im * 9;
    b[0] = 10.0;
    b[4] = 10.0;
    b[8] = 10.0;
  }
}

static int surface_fail(void *, eon_relax_surface_request_t *) { return -3; }

static std::vector<uint8_t> pack_relax_params(const char *kind,
                                              int64_t saddle_iters = 20,
                                              bool climbing = true) {
  ::capnp::MallocMessageBuilder msg;
  auto root = msg.initRoot<eonc::params_ssot::RelaxEngineParams>();
  root.setKind(kind);
  auto neb = root.getNeb();
  neb.setMinimizeEndpoints(true);
  neb.setClimbingImage(climbing);
  auto saddle = root.getSaddle();
  saddle.setMaxIterations(saddle_iters);
  const auto words = ::capnp::messageToFlatArray(msg);
  const auto bytes = words.asBytes();
  return {bytes.begin(), bytes.end()};
}

TEST_CASE("HostSurfacePotential forceBatch counts each system",
          "[relax][pot]") {
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  eonc::HostSurfacePotential pot(surface_forward, &ctx, 0);
  double pos0[3]{0.1, 0.0, 0.0};
  double pos1[3]{0.2, 0.0, 0.0};
  double frc0[3]{};
  double frc1[3]{};
  double box0[9]{10.0, 0, 0, 0, 10.0, 0, 0, 0, 10.0};
  double box1[9]{10.0, 0, 0, 0, 10.0, 0, 0, 0, 10.0};
  int nrs0[1]{1};
  int nrs1[1]{1};
  const double *positions[] = {pos0, pos1};
  const int *atomicNrs[] = {nrs0, nrs1};
  double *forces[] = {frc0, frc1};
  const double *boxes[] = {box0, box1};
  double energies[2]{};
  double variances[2]{};
  pot.forceBatch(2, 1, positions, atomicNrs, forces, energies, variances,
                 boxes);
  REQUIRE(pot.forceCallCounter.load() == 2);
  REQUIRE(ctx.calls == 1);
}

TEST_CASE("HostSurfacePotential force plus Matter counts one call",
          "[relax][pot][epoch]") {
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  auto pot = std::make_shared<eonc::HostSurfacePotential>(surface_forward,
                                                          &ctx, 0);
  Parameters params;
  Matter m(pot, params);
  m.resize(1);
  Eigen::VectorXi z(1);
  z(0) = 1;
  m.setAtomicNrs(z);
  m.setMasses(Eigen::VectorXd::Ones(1));
  Matrix3d cell = Matrix3d::Zero();
  cell(0, 0) = 10.0;
  cell(1, 1) = 10.0;
  cell(2, 2) = 10.0;
  m.setPeriodic(false);
  m.setCell(cell);
  AtomMatrix p(1, 3);
  p << 0.2, 0.0, 0.0;
  m.setPositions(p);
  ctx.calls = 0;
  const double e1 = m.getPotentialEnergy();
  REQUIRE(ctx.calls == 1);
  REQUIRE(pot->forceCallCounter.load() == 1);
  REQUIRE(m.getPotentialEnergy() == e1);
  REQUIRE(ctx.calls == 1);
  pot->setEpoch(7);
  REQUIRE(m.getPotentialEnergy() == e1);
  REQUIRE(ctx.calls == 2);
}

TEST_CASE("NEBObjectiveFunction returns interior host masses",
          "[relax][neb][masses]") {
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  auto pot = std::make_shared<eonc::HostSurfacePotential>(surface_forward,
                                                          &ctx, 0);
  Parameters params;
  ParametersLoadAccess::neb_options(params).image_count = 5;
  std::vector<Matter> path;
  path.reserve(7);
  for (int im = 0; im < 7; ++im) {
    Matter m(pot, params);
    m.resize(2);
    Eigen::VectorXi z(2);
    z << 1, 1;
    m.setAtomicNrs(z);
    Eigen::VectorXd masses(2);
    masses << 1.0, 16.0;
    m.setMasses(masses);
    Matrix3d cell = Matrix3d::Zero();
    cell(0, 0) = 10.0;
    cell(1, 1) = 10.0;
    cell(2, 2) = 10.0;
    m.setPeriodic(false);
    m.setCell(cell);
    AtomMatrix p(2, 3);
    p.setZero();
    p(0, 0) = 0.3 - 0.1 * im;
    p(1, 0) = 0.1;
    m.setPositions(p);
    path.push_back(std::move(m));
  }
  auto neb = std::make_unique<NudgedElasticBand>(std::move(path), params, pot);
  NEBObjectiveFunction objf(neb.get(), params);
  const VectorXd masses = objf.getMasses();
  REQUIRE(masses.size() == 10);
  REQUIRE(masses(0) == 1.0);
  REQUIRE(masses(1) == 16.0);
  REQUIRE(masses(8) == 1.0);
  REQUIRE(masses(9) == 16.0);
}

TEST_CASE("relax engine sheared box keeps RowMajor data order",
          "[relax][neb][box]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  for (long im = 0; im < n_images; ++im) {
    boxes[static_cast<size_t>(im * 9 + 1)] = 0.25;
  }
  struct BoxCtx {
    SurfaceCtx fwd;
    double seen{0.0};
  } bctx{};
  auto record = [](void *user, eon_relax_surface_request_t *req) -> int {
    auto *c = static_cast<BoxCtx *>(user);
    if (req->boxes) {
      c->seen = req->boxes[1];
    }
    return surface_forward(&c->fwd, req);
  };
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_step(eng, &band, record, &bctx, &out) == EON_RELAX_OK);
  REQUIRE(bctx.seen == 0.25);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine first surface keeps host-frame negatives",
          "[relax][neb][pbc]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  double min_x = 0.0;
  auto record = [](void *user, eon_relax_surface_request_t *req) -> int {
    double *mn = static_cast<double *>(user);
    const long n = static_cast<long>(req->n_images * req->n_atoms);
    for (long i = 0; i < n; ++i) {
      *mn = std::min(*mn, req->positions[static_cast<size_t>(3 * i)]);
    }
    SurfaceCtx ctx{nullptr, 0.0, 0, 0};
    return surface_forward(&ctx, req);
  };
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_step(eng, &band, record, &min_x, &out) == EON_RELAX_OK);
  REQUIRE(min_x < -0.25);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine ABI stamp is 2.1", "[relax][abi]") {
  REQUIRE(eon_relax_abi_version() ==
          static_cast<int>((EON_RELAX_ABI_MAJOR << 16) | EON_RELAX_ABI_MINOR));
  REQUIRE(eon_relax_available() == 1);
  eon_relax_version_t stamp{};
  REQUIRE(eon_relax_abi_stamp(&stamp) == EON_RELAX_OK);
  REQUIRE(stamp.major == 2);
  REQUIRE(stamp.minor == 1);
  REQUIRE(eon_relax_abi_stamp(nullptr) == EON_RELAX_INVALID_PARAMETER);
  REQUIRE(eon_relax_version_hash() != 0);
  const char *id = eon_relax_version_hash_str();
  REQUIRE(id != nullptr);
  REQUIRE(std::string(id).find("+git.") != std::string::npos);
  REQUIRE(std::string(id).find("+feat.") != std::string::npos);
  const auto git_pos = std::string(id).find("+git.");
  REQUIRE(git_pos != std::string::npos);
  const auto feat_pos = std::string(id).find("+feat.", git_pos);
  REQUIRE(feat_pos != std::string::npos);
  REQUIRE(feat_pos - (git_pos + 5) == 40);
}

TEST_CASE("relax engine destroy NULL is a no-op", "[relax][abi]") {
  eon_relax_destroy(nullptr);
  REQUIRE(eon_relax_last_error(nullptr) == EON_RELAX_NULL_ENGINE);
  REQUIRE(eon_relax_set_surface_epoch(nullptr, 1) == EON_RELAX_NULL_ENGINE);
  REQUIRE(eon_relax_reset(nullptr) == EON_RELAX_NULL_ENGINE);
}

TEST_CASE("relax engine status_name is fail-closed", "[relax][abi]") {
  REQUIRE(std::string(eon_relax_status_name(EON_RELAX_KIND_NEB,
                                            EON_RELAX_NEB_GOOD)) == "GOOD");
  REQUIRE(std::string(eon_relax_status_name(
              EON_RELAX_KIND_NEB, EON_RELAX_NEB_MAX_UNCERTAINTY)) ==
          "MAX_UNCERTAINTY");
  REQUIRE(eon_relax_status_name(EON_RELAX_KIND_NEB, 99) == nullptr);
  REQUIRE(eon_relax_status_name(EON_RELAX_KIND_INVALID, 0) == nullptr);
}

TEST_CASE("relax engine create NULL config and reject unknown capnp",
          "[relax][abi]") {
  char err[128]{};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, err, sizeof(err));
  REQUIRE(eng != nullptr);
  const char blob[] = "not-capnp";
  REQUIRE(eon_relax_create(blob, sizeof(blob), err, sizeof(err)) == nullptr);
  REQUIRE(std::string(err).find("-13") != std::string::npos);
  REQUIRE(eon_relax_create(blob, 0, err, sizeof(err)) == nullptr);
  REQUIRE(eon_relax_run(nullptr, nullptr, nullptr, nullptr, nullptr) ==
          EON_RELAX_NULL_ENGINE);
  eon_relax_band_t short_band{};
  short_band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  short_band.n_images = 2;
  short_band.n_atoms = 1;
  double pos1[6]{};
  int32_t z1[1]{1};
  double box1[18]{};
  short_band.positions = pos1;
  short_band.atomic_nrs = z1;
  short_band.boxes = box1;
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_BAND_TOO_SHORT);
  short_band.n_images = 4;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_BAND_SIZE);
  REQUIRE(eon_relax_last_error(eng) == EON_RELAX_BAND_SIZE);
  short_band.n_images = 7;
  short_band.n_atoms = 0;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_NATOMS);
  short_band.n_atoms = 1;
  short_band.positions = nullptr;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_NULL_POSITIONS);
  short_band.positions = pos1;
  short_band.atomic_nrs = nullptr;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_NULL_ATOMIC_NRS);
  short_band.atomic_nrs = z1;
  short_band.boxes = nullptr;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_NULL_BOXES);
  short_band.boxes = box1;
  short_band.version.major = 0;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_ABI_MISMATCH);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine NEB on a harmonic surface converges",
          "[relax][neb]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);

  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  const int rc = eon_relax_run(eng, &band, surface_forward, &ctx, &out);
  REQUIRE(rc == EON_RELAX_OK);
  REQUIRE(out.status == EON_RELAX_NEB_GOOD);
  REQUIRE(out.kind == EON_RELAX_KIND_NEB);
  REQUIRE(std::isfinite(out.max_force));
  REQUIRE(out.max_force <= 0.01 + 1e-6);
  REQUIRE(ctx.calls > 0);
  REQUIRE(out.version_hash != 0);
  // Endpoints stay put: NudgedElasticBandJob owns endpoint min, not compute().
  REQUIRE(std::abs(pos[0] - 0.3) < 1e-6);
  REQUIRE(std::abs(pos[3]) <= 0.2 + 1e-6);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine MAX_UNCERTAINTY surfaces from host variance",
          "[relax][neb][uncertainty]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);

  SurfaceCtx ctx{nullptr, 10.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  const int rc = eon_relax_run(eng, &band, surface_forward, &ctx, &out);
  REQUIRE(rc == EON_RELAX_OK);
  REQUIRE(out.status == EON_RELAX_NEB_MAX_UNCERTAINTY);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine surface failure is fail-closed", "[relax][abi]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);

  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  const int rc = eon_relax_run(eng, &band, surface_fail, nullptr, &out);
  REQUIRE(rc == EON_RELAX_SURFACE_FATAL);
  eon_relax_destroy(eng);
}

TEST_CASE("libeon_relax_engine exports the C waist", "[relax][dlopen]") {
  const char *env = std::getenv("EON_RELAX_ENGINE");
  std::vector<std::string> candidates;
  if (env && *env) {
    candidates.emplace_back(env);
  }
  candidates.emplace_back("libeon_relax_engine.so");
  eonc::dynlib::Handle h{};
  for (const auto &path : candidates) {
    h = eonc::dynlib::open(path.c_str());
    if (h) {
      break;
    }
  }
  REQUIRE(h);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_create") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_run") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_destroy") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_status_name") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_last_error") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_version_hash_str") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_abi_version") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_step") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_reset") != nullptr);
  auto abi = eonc::dynlib::loadSym<int (*)()>(h, "eon_relax_abi_version");
  REQUIRE(abi != nullptr);
  REQUIRE(abi() == static_cast<int>((EON_RELAX_ABI_MAJOR << 16) |
                                    EON_RELAX_ABI_MINOR));
  eonc::dynlib::close(h);
}


TEST_CASE("relax engine stepper converges the harmonic band one step at a "
          "time",
          "[relax][neb][stepper]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);

  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();

  eon_relax_outcome_t out{};
  int steps = 0;
  int rc = EON_RELAX_OK;
  do {
    rc = eon_relax_step(eng, &band, surface_forward, &ctx, &out);
    REQUIRE(rc == EON_RELAX_OK);
    ++steps;
  } while (out.status == EON_RELAX_NEB_RUNNING && steps < 500);
  REQUIRE(out.status == EON_RELAX_NEB_GOOD);
  REQUIRE(out.iterations <= steps);
  REQUIRE(std::isfinite(out.max_force));
  // Endpoints are the caller's; the stepper never moves them.
  REQUIRE(std::abs(pos[0] - 0.3) < 1e-12);
  REQUIRE(std::abs(pos[(n_images - 1) * 3] + 0.3) < 1e-12);

  // A different surface pointer without a reset is refused.
  REQUIRE(eon_relax_step(eng, &band, surface_fail, &ctx, &out) ==
          EON_RELAX_INVALID_PARAMETER);
  // Reset drops the band state; the next step reinitializes.
  REQUIRE(eon_relax_reset(eng) == EON_RELAX_OK);
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine stepper surfaces MAX_UNCERTAINTY",
          "[relax][neb][stepper][uncertainty]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  SurfaceCtx ctx{nullptr, 10.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  REQUIRE(out.status == EON_RELAX_NEB_MAX_UNCERTAINTY);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine create parses RelaxEngineParams kind tokens",
          "[relax][abi][capnp]") {
  char err[128]{};
  const auto neb = pack_relax_params("neb");
  EonRelaxEngine *eng =
      eon_relax_create(neb.data(), neb.size(), err, sizeof(err));
  REQUIRE(eng != nullptr);
  eon_relax_destroy(eng);

  const auto sad = pack_relax_params("saddle");
  eng = eon_relax_create(sad.data(), sad.size(), err, sizeof(err));
  REQUIRE(eng != nullptr);
  eon_relax_destroy(eng);

  const auto bad = pack_relax_params("dimer");
  REQUIRE(eon_relax_create(bad.data(), bad.size(), err, sizeof(err)) ==
          nullptr);
  REQUIRE(std::string(err).find("-11") != std::string::npos);
}

TEST_CASE("relax engine saddle kind requires one image and a mode",
          "[relax][saddle]") {
  const auto sad = pack_relax_params("saddle");
  EonRelaxEngine *eng = eon_relax_create(sad.data(), sad.size(), nullptr, 0);
  REQUIRE(eng != nullptr);

  double pos[3]{0.1, 0.0, 0.0};
  double boxes[9]{10, 0, 0, 0, 10, 0, 0, 0, 10};
  int32_t z[1]{1};
  double mode[3]{1.0, 0.0, 0.0};
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = 2;
  band.n_atoms = 1;
  band.positions = pos;
  band.atomic_nrs = z;
  band.boxes = boxes;
  band.mode = mode;
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_run(eng, &band, surface_forward, nullptr, &out) ==
          EON_RELAX_SADDLE_NIMAGES);
  REQUIRE(out.status == -1);
  REQUIRE(out.version.major == 2);
  band.n_images = 1;
  band.mode = nullptr;
  REQUIRE(eon_relax_run(eng, &band, surface_forward, nullptr, &out) ==
          EON_RELAX_NULL_MODE);
  band.mode = mode;
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  const int rc = eon_relax_run(eng, &band, surface_forward, &ctx, &out);
  REQUIRE(rc == EON_RELAX_OK);
  REQUIRE(out.kind == EON_RELAX_KIND_SADDLE);
  REQUIRE(eon_relax_status_name(EON_RELAX_KIND_SADDLE, out.status) != nullptr);
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_INVALID_PARAMETER);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine step BAND_* stamps a zero leftover outcome",
          "[relax][abi][stepper]") {
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  double pos[6]{};
  int32_t z[1]{1};
  double boxes[18]{};
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = 2;
  band.n_atoms = 1;
  band.positions = pos;
  band.atomic_nrs = z;
  band.boxes = boxes;
  eon_relax_outcome_t out{};
  out.iterations = 99;
  out.climbing_image = 7;
  out.max_force = 3.14;
  REQUIRE(eon_relax_step(eng, &band, surface_fail, nullptr, &out) ==
          EON_RELAX_BAND_TOO_SHORT);
  REQUIRE(out.version.major == 2);
  REQUIRE(out.version.minor == 1);
  REQUIRE(out.status == -1);
  REQUIRE(out.iterations == 0);
  REQUIRE(out.climbing_image == 0);
  REQUIRE(out.max_force == 0.0);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine SURFACE_FATAL makes the instance unusable",
          "[relax][abi]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  long calls = 0;
  auto count_fail = [](void *user, eon_relax_surface_request_t *) -> int {
    ++*static_cast<long *>(user);
    return -3;
  };
  REQUIRE(eon_relax_run(eng, &band, count_fail, &calls, &out) ==
          EON_RELAX_SURFACE_FATAL);
  const long after_fatal = calls;
  REQUIRE(after_fatal > 0);
  REQUIRE(eon_relax_run(eng, &band, count_fail, &calls, &out) ==
          EON_RELAX_SURFACE_FATAL);
  REQUIRE(calls == after_fatal);
  REQUIRE(eon_relax_step(eng, &band, count_fail, &calls, &out) ==
          EON_RELAX_SURFACE_FATAL);
  REQUIRE(calls == after_fatal);
  REQUIRE(eon_relax_set_surface_epoch(eng, 9) == EON_RELAX_SURFACE_FATAL);
  REQUIRE(eon_relax_last_error(eng) == EON_RELAX_SURFACE_FATAL);
  REQUIRE(eon_relax_reset(eng) == EON_RELAX_SURFACE_FATAL);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine epoch_out busts Matter energy and variance",
          "[relax][neb][epoch]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  const long first = ctx.calls;
  REQUIRE(first > 0);
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  const long second = ctx.calls;
  ctx.epoch = 7;
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  REQUIRE(ctx.calls > second);
  REQUIRE(out.surface_epoch == 7);
  REQUIRE(first > 0);
  eon_relax_destroy(eng);
}

static int surface_product_var(void *user, eon_relax_surface_request_t *req) {
  auto *ctx = static_cast<SurfaceCtx *>(user);
  ++ctx->calls;
  const long n_images = static_cast<long>(req->n_images);
  const long n_atoms = static_cast<long>(req->n_atoms);
  for (long s = 0; s < n_images; ++s) {
    const double *p = req->positions + s * 3 * n_atoms;
    double *f = req->forces + s * 3 * n_atoms;
    double e = 0.0;
    for (long a = 0; a < n_atoms; ++a) {
      const double x = p[3 * a];
      const double y = p[3 * a + 1];
      const double z = p[3 * a + 2];
      e += 0.5 * (x * x + y * y + z * z);
      f[3 * a] = -x;
      f[3 * a + 1] = -y;
      f[3 * a + 2] = -z;
    }
    req->energies[s] = e;
    if (req->variances) {
      const int64_t id = req->image_ids ? req->image_ids[s] : EON_RELAX_IMAGE_NONE;
      req->variances[s] = (id == ctx->product_id) ? 10.0 : 0.0;
    }
  }
  if (req->epoch_out) {
    *req->epoch_out = ctx->epoch;
  }
  return 0;
}

TEST_CASE("relax engine MAX_UNCERTAINTY sees the product image",
          "[relax][neb][uncertainty]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  SurfaceCtx ctx{};
  ctx.product_id = n_images - 1;
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  const int rc = eon_relax_run(eng, &band, surface_product_var, &ctx, &out);
  REQUIRE(rc == EON_RELAX_OK);
  REQUIRE(out.status == EON_RELAX_NEB_MAX_UNCERTAINTY);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine copies host masses and refuses non-positive",
          "[relax][abi][masses]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  double masses[1]{-1.0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  band.masses = masses;
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_run(eng, &band, surface_forward, nullptr, &out) ==
          EON_RELAX_INVALID_PARAMETER);
  masses[0] = 16.0;
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  REQUIRE(eon_relax_run(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  REQUIRE(out.status == EON_RELAX_NEB_GOOD);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine NULL is_fixed frees a previously fixed atom",
          "[relax][abi][fixed]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  int32_t fixed[1]{1};
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  band.is_fixed = fixed;
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  const long calls0 = ctx.calls;
  REQUIRE(calls0 > 0);
  band.is_fixed = nullptr;
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  REQUIRE(ctx.calls > calls0);
  eon_relax_destroy(eng);
}

static int surface_recover(void *, eon_relax_surface_request_t *) { return 5; }

TEST_CASE("relax engine NULL_BAND stamps leftover outcome", "[relax][abi]") {
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_outcome_t out{};
  out.status = EON_RELAX_NEB_GOOD;
  out.iterations = 9;
  REQUIRE(eon_relax_run(eng, nullptr, surface_fail, nullptr, &out) ==
          EON_RELAX_NULL_BAND);
  REQUIRE(out.status == -1);
  REQUIRE(out.version.major == 2);
  REQUIRE(out.iterations == 0);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine first-step recoverable does not leave a dead stepper",
          "[relax][neb][stepper]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_step(eng, &band, surface_recover, nullptr, &out) == 5);
  REQUIRE(out.status == -1);
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine live-step recoverable drops QN history",
          "[relax][neb][stepper]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  ctx.recover_on_call = ctx.calls + 1;
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) == 5);
  REQUIRE(out.status == -1);
  ctx.recover_on_call = 0;
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine set_surface_epoch updates a live pot",
          "[relax][neb][epoch]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  REQUIRE(eon_relax_set_surface_epoch(eng, 9) == EON_RELAX_OK);
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  REQUIRE(ctx.seen_epoch == 9);
  REQUIRE(out.surface_epoch == 9);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine live step refuses a changed atom count",
          "[relax][neb][stepper]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  band.n_atoms = 2;
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_INVALID_PARAMETER);
  REQUIRE(out.status == -1);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine run leaves caller endpoints in the caller frame",
          "[relax][neb]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);
  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_run(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  REQUIRE(std::abs(pos[0] - 0.3) < 1e-12);
  REQUIRE(std::abs(pos[(n_images - 1) * 3] + 0.3) < 1e-12);
  eon_relax_destroy(eng);
}

static int two_atom_saddle_surface(void *, eon_relax_surface_request_t *req) {
  const long n_images = static_cast<long>(req->n_images);
  const long n_atoms = static_cast<long>(req->n_atoms);
  for (long s = 0; s < n_images; ++s) {
    const double *p = req->positions + s * 3 * n_atoms;
    double *f = req->forces + s * 3 * n_atoms;
    for (long k = 0; k < 3 * n_atoms; ++k) {
      f[k] = 0.0;
    }
    const double x0 = p[0];
    const double y0 = p[1] - 1.0;
    const double x1 = n_atoms > 1 ? p[3] : 0.0;
    const double y1 = n_atoms > 1 ? p[4] : 0.0;
    const double xs0 = x0 * x0 - 1.0;
    const double xs1 = x1 * x1 - 1.0;
    req->energies[s] =
        0.5 * xs0 * xs0 + 0.5 * xs1 * xs1 + 0.5 * (y0 * y0 + y1 * y1);
    f[0] = -2.0 * x0 * xs0;
    f[1] = -y0;
    if (n_atoms > 1) {
      f[3] = -2.0 * x1 * xs1;
      f[4] = -y1;
    }
    if (req->variances) {
      req->variances[s] = 0.0;
    }
  }
  return 0;
}

TEST_CASE("relax engine stepper keeps the climbing image on a local max",
          "[relax][neb][stepper][ci]") {
  const long n_images = 7;
  const long n_atoms = 2;
  std::vector<double> pos(static_cast<size_t>(n_images * 3 * n_atoms), 0.0);
  std::vector<double> boxes(static_cast<size_t>(n_images * 9), 0.0);
  std::vector<int32_t> z = {1, 1};
  for (long im = 0; im < n_images; ++im) {
    const double t =
        static_cast<double>(im) / static_cast<double>(n_images - 1);
    const double x = 2.0 - 4.0 * t;
    pos[static_cast<size_t>(im * 6)] = x;
    pos[static_cast<size_t>(im * 6 + 3)] = x;
    double *b = boxes.data() + im * 9;
    b[0] = b[4] = b[8] = 10.0;
  }
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = n_atoms;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_step(eng, &band, two_atom_saddle_surface, nullptr, &out) ==
          EON_RELAX_OK);
  const int64_t first_ci = out.climbing_image;
  REQUIRE(first_ci >= 1);
  REQUIRE(first_ci <= 5);
  REQUIRE(eon_relax_step(eng, &band, two_atom_saddle_surface, nullptr, &out) ==
          EON_RELAX_OK);
  REQUIRE(out.climbing_image == first_ci);
  REQUIRE(eon_relax_step(eng, &band, two_atom_saddle_surface, nullptr, &out) ==
          EON_RELAX_OK);
  REQUIRE(out.climbing_image == first_ci);
  for (long im = 1; im + 1 < n_images; ++im) {
    const double x = pos[static_cast<size_t>(im * 6)];
    REQUIRE(x < 2.5);
    REQUIRE(x > -2.5);
  }
  eon_relax_destroy(eng);
}

} // namespace tests

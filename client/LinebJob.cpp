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
#include "eon/LinebJob.h"

#include "eon/ConFileIO.h"
#include "eon/EonLogger.h"
#include "eon/JobResult.h"
#include "eon/LinebSeed.h"
#include "eon/Matter.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"
#include "eon/potentials/RgpotAdapter/ResultCache.h"

#include <cmath>
#include <limits>
#include <map>
#include <nlohmann/json.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef WITH_RGSADDLE_LINEB
#include <rgsaddle.h>
#endif

namespace eonc {

bool linebAvailable() {
#ifdef WITH_RGSADDLE_LINEB
  rgsaddle_version_t abi{};
  return rgsaddle_abi_stamp(&abi) == RGSADDLE_OK &&
         abi.major == RGSADDLE_ABI_MAJOR;
#else
  return false;
#endif
}

namespace {

#ifdef WITH_RGSADDLE_LINEB

struct LinebSurface {
  Matter *probe = nullptr;
  std::vector<long> freeAtoms;
  int64_t surfaceRows = 0;
};

/// Map free-atom positions into the full Matter and write energy plus
/// energy gradients (minus forces) for those atoms. Same sign convention
/// as XtsciIrc / Instanton.
rgsaddle_status_t surfaceCallback(void *user, rgsaddle_surface_request_t *req) {
  auto *ctx = static_cast<LinebSurface *>(user);
  if (ctx == nullptr || ctx->probe == nullptr || req == nullptr ||
      req->positions == nullptr || req->energies == nullptr ||
      req->gradients == nullptr) {
    return RGSADDLE_SURFACE_FAILED;
  }
  const auto nFree = static_cast<int64_t>(ctx->freeAtoms.size());
  if (req->n_atoms != nFree) {
    return RGSADDLE_SHAPE;
  }

  const int64_t nImages =
      (req->flags & RGSADDLE_REQ_ONE_IMAGE) ? 1 : req->n_images;
  if (nImages < 1) {
    return RGSADDLE_SHAPE;
  }

  try {
    AtomMatrix pos = ctx->probe->getPositions();
    for (int64_t image = 0; image < nImages; ++image) {
      const double *xyz = req->positions + image * nFree * 3;
      for (int64_t a = 0; a < nFree; ++a) {
        const long atom = ctx->freeAtoms[static_cast<size_t>(a)];
        pos(atom, 0) = xyz[3 * a + 0];
        pos(atom, 1) = xyz[3 * a + 1];
        pos(atom, 2) = xyz[3 * a + 2];
      }
      ctx->probe->setPositions(pos);
      req->energies[image] = ctx->probe->getPotentialEnergy();
      const AtomMatrix forces = ctx->probe->getForces();
      double *grad = req->gradients + image * nFree * 3;
      for (int64_t a = 0; a < nFree; ++a) {
        const long atom = ctx->freeAtoms[static_cast<size_t>(a)];
        grad[3 * a + 0] = -forces(atom, 0);
        grad[3 * a + 1] = -forces(atom, 1);
        grad[3 * a + 2] = -forces(atom, 2);
      }
    }
    ctx->surfaceRows += nImages;
    return RGSADDLE_OK;
  } catch (...) {
    return RGSADDLE_SURFACE_FAILED;
  }
}

#ifdef RGSADDLE_HAS_LINEB_ORBIT
void writeRateOrbit(const std::string &filename, const Matter &seed,
                    const std::vector<long> &freeAtoms,
                    const std::vector<double> &positions,
                    const std::vector<double> &gradients,
                    nlohmann::json metadata) {
  Matter bead(seed);
  // Preserve the continuous Cartesian path across periodic cell boundaries.
  // The copy is only serialized; the potential is never called on it.
  bead.setPeriodic(false);
  AtomMatrix xyz = bead.getPositions();
  const size_t width = 3 * freeAtoms.size();
  std::vector<readcon::ConFrame> frames;
  frames.reserve(positions.size() / width);
  auto &observed = metadata["lineb_orbit"];
  for (const auto atom : freeAtoms) {
    observed["mobile_atom_ids"].push_back(bead.getAtomIndex(atom));
  }
  for (size_t offset = 0; offset < positions.size(); offset += width) {
    auto &rows = observed["gradients_eV_per_angstrom"];
    rows = nlohmann::json::array();
    for (size_t a = 0; a < freeAtoms.size(); ++a) {
      const size_t k = offset + 3 * a;
      xyz.row(freeAtoms[a]) =
          Eigen::Map<const Eigen::RowVector3d>(positions.data() + k);
      rows.push_back({gradients[k], gradients[k + 1], gradients[k + 2]});
    }
    bead.setPositions(xyz);
    io::ConFrameMetadata frame;
    frame.frame_index = frames.size();
    frame.write_con_forces = false;
    frame.raw_json = metadata.dump();
    frames.push_back(io::matterToConFrame(bead, &frame));
  }
  if (!io::io_ok(io::writeConFrames(filename, frames))) {
    throw std::runtime_error("lineb: cannot write " + filename);
  }
}
#endif

void unpackPositions(const std::vector<double> &flat,
                     const std::vector<long> &freeAtoms,
                     std::vector<std::shared_ptr<Matter>> &band) {
  const auto nFree = freeAtoms.size();
  if (nFree == 0 || flat.size() % (nFree * 3) != 0) {
    throw std::runtime_error("lineb: optimized path has the wrong length");
  }
  const size_t nImages = flat.size() / (nFree * 3);
  if (nImages != band.size()) {
    throw std::runtime_error("lineb: optimized path image count changed");
  }
  for (size_t i = 0; i < nImages; ++i) {
    AtomMatrix pos = band[i]->getPositions();
    for (size_t a = 0; a < nFree; ++a) {
      const long atom = freeAtoms[a];
      const size_t base = (i * nFree + a) * 3;
      pos(atom, 0) = flat[base + 0];
      pos(atom, 1) = flat[base + 1];
      pos(atom, 2) = flat[base + 2];
    }
    band[i]->setPositions(pos);
  }
}

#endif // WITH_RGSADDLE_LINEB

} // namespace

std::vector<std::string> LinebJob::run(void) {
  if (!linebAvailable()) {
    throw std::runtime_error(
        "lineb: needs -Dwith_rgsaddle=true against a rgsaddle that exports "
        "rgsaddle_lineb_solve (lineb C API)");
  }

#ifdef WITH_RGSADDLE_LINEB
  const auto &o = params.lineb_options();
  if (!(o.energy < std::numeric_limits<double>::infinity()) ||
      !std::isfinite(o.energy)) {
    throw std::runtime_error("lineb: [Lineb] energy must be a finite tunneling "
                             "energy below the barrier");
  }
  if (o.max_steps < 1) {
    throw std::runtime_error("lineb: [Lineb] max_steps must be at least 1");
  }
  if (o.beads < 0) {
    throw std::runtime_error(
        "lineb: [Lineb] beads must be nonnegative (0 is path only)");
  }

#ifdef RGPOT_HAS_CACHE
  std::optional<ScopedRgpotCache> resultCache;
  if (!o.result_cache.empty()) {
    auto *provider = dynamic_cast<RgpotCacheProvider *>(pot.get());
    if (!provider) {
      throw std::runtime_error("lineb: result_cache requires an rgpot adapter");
    }
    resultCache.emplace(*provider, o.result_cache);
  }
  const auto cacheCounts = [&] {
    return resultCache ? resultCache->counts() : rgpot::cache::EvalCounts{};
  };
#else
  if (!o.result_cache.empty()) {
    throw std::runtime_error(
        "lineb: result_cache requires rgpot built with cache support");
  }
#endif
  const auto jobStart = pot->forceCallCounter.load();
  auto band = lineb::loadSeedBand(params, pot, o.initial_path);
  std::vector<long> freeAtoms;
  std::vector<double> masses;
  lineb::collectFreeAtoms(*band.front(), freeAtoms, masses);
  for (size_t i = 1; i < band.size(); ++i) {
    std::vector<long> otherFree;
    std::vector<double> otherMasses;
    lineb::collectFreeAtoms(*band[i], otherFree, otherMasses);
    if (otherFree != freeAtoms) {
      throw std::runtime_error(
          "lineb: every seed image must share the same free-atom mask");
    }
  }

  std::vector<double> positions;
  // Unwrapped Cartesian for the C API only; Matter stays on its wrapped seed.
  lineb::packSeedCartesian(band, freeAtoms, o.molecular, positions);
  std::vector<double> optimized(positions.size(), 0.0);

  rgsaddle_lineb_config_t cfg{};
  const rgsaddle_status_t cfgRc =
      rgsaddle_lineb_config_default(o.opt_tol, o.perp_tol, &cfg);
  if (cfgRc != RGSADDLE_OK) {
    throw std::runtime_error(
        std::string("lineb: rgsaddle_lineb_config_default: ") +
        rgsaddle_status_name(cfgRc));
  }
  cfg.spring = o.spring;
  cfg.endpoint_scale = o.endpoint_scale;
  cfg.time_step = o.time_step;
  cfg.max_move = o.max_move;
  cfg.turning_width = o.turning_width;
  cfg.hbar = o.hbar;
  cfg.boltzmann = o.boltzmann;
  cfg.quadrature_pairs = o.quadrature_pairs;

  Matter probe(*band.front());
  LinebSurface ctx{&probe, freeAtoms};

#ifdef RGPOT_HAS_CACHE
  const auto cacheAtSolve = cacheCounts();
#endif
  const auto solveStart = pot->forceCallCounter.load();
  rgsaddle_lineb_result_t result{};
#ifndef RGSADDLE_LINEB_PROJECTED_VERLET
  if (o.projected_verlet) {
    throw std::runtime_error(
        "lineb: projected_verlet requires rgsaddle with RGSADDLE_LINEB_PROJECTED_VERLET");
  }
#endif
#ifndef RGSADDLE_LINEB_LINEAR_SEED
  if (o.linear_seed || o.endpoint_convergence) {
    throw std::runtime_error("lineb: linear_seed and endpoint_convergence "
                             "require rgsaddle with fixed-energy seed options");
  }
#endif
#ifdef RGSADDLE_LINEB_LINEAR_SEED
  uint32_t options =
      (o.linear_seed ? RGSADDLE_LINEB_LINEAR_SEED : 0u) |
      (o.endpoint_convergence ? RGSADDLE_LINEB_ENDPOINT_CONVERGENCE : 0u);
#ifdef RGSADDLE_LINEB_PROJECTED_VERLET
  if (o.projected_verlet) {
    options |= RGSADDLE_LINEB_PROJECTED_VERLET;
  }
#endif
#ifdef WITH_RGSADDLE_LINEB_COUNTS
  rgsaddle_lineb_counts_t counts{};
#ifdef RGSADDLE_HAS_LINEB_RATE_DIAGNOSTICS
  int64_t negativeModes = -1;
  int64_t retainedModes = -1;
#ifdef RGSADDLE_HAS_LINEB_ORBIT
  const size_t width = 3 * freeAtoms.size();
  if (static_cast<size_t>(o.beads) > std::vector<double>().max_size() / width) {
    throw std::length_error("lineb: orbit exceeds addressable memory");
  }
  std::vector<double> orbitPositions(static_cast<size_t>(o.beads) * width);
  std::vector<double> orbitGradients(orbitPositions.size());
  rgsaddle_lineb_orbit_t orbit{};
  orbit.capacity = orbitPositions.size();
  orbit.positions = orbitPositions.data();
  orbit.gradients = orbitGradients.data();
  const rgsaddle_status_t rc = rgsaddle_lineb_solve_counted_orbit(
      &cfg, static_cast<int64_t>(band.size()),
      static_cast<int64_t>(freeAtoms.size()), positions.data(), masses.data(),
      o.molecular ? 1 : 0, o.energy, o.max_steps, o.beads, surfaceCallback,
      &ctx, optimized.data(), &result,
      o.endpoint_kappa > 0.0 ? &o.endpoint_kappa : nullptr, options, &counts,
      &negativeModes, &retainedModes, o.beads ? &orbit : nullptr);
#else
  const rgsaddle_status_t rc = rgsaddle_lineb_solve_counted_diagnostics(
      &cfg, static_cast<int64_t>(band.size()),
      static_cast<int64_t>(freeAtoms.size()), positions.data(), masses.data(),
      o.molecular ? 1 : 0, o.energy, o.max_steps, o.beads, surfaceCallback,
      &ctx, optimized.data(), &result,
      o.endpoint_kappa > 0.0 ? &o.endpoint_kappa : nullptr, options, &counts,
      &negativeModes, &retainedModes);
#endif
#else
  const rgsaddle_status_t rc = rgsaddle_lineb_solve_counted(
      &cfg, static_cast<int64_t>(band.size()),
      static_cast<int64_t>(freeAtoms.size()), positions.data(), masses.data(),
      o.molecular ? 1 : 0, o.energy, o.max_steps, o.beads, surfaceCallback,
      &ctx, optimized.data(), &result,
      o.endpoint_kappa > 0.0 ? &o.endpoint_kappa : nullptr, options, &counts);
#endif
#else
  const rgsaddle_status_t rc = rgsaddle_lineb_solve_options(
      &cfg, static_cast<int64_t>(band.size()),
      static_cast<int64_t>(freeAtoms.size()), positions.data(), masses.data(),
      o.molecular ? 1 : 0, o.energy, o.max_steps, o.beads, surfaceCallback,
      &ctx, optimized.data(), &result,
      o.endpoint_kappa > 0.0 ? &o.endpoint_kappa : nullptr, options);
#endif
#elif defined(WITH_RGSADDLE_LINEB_KAPPA)
  const rgsaddle_status_t rc = rgsaddle_lineb_solve_with_kappa(
      &cfg, static_cast<int64_t>(band.size()),
      static_cast<int64_t>(freeAtoms.size()), positions.data(), masses.data(),
      o.molecular ? 1 : 0, o.energy, o.max_steps, o.beads, surfaceCallback,
      &ctx, optimized.data(), &result,
      o.endpoint_kappa > 0.0 ? &o.endpoint_kappa : nullptr);
#else
  if (o.endpoint_kappa > 0.0) {
    throw std::runtime_error(
        "lineb: endpoint_kappa requires rgsaddle_lineb_solve_with_kappa");
  }
  const rgsaddle_status_t rc = rgsaddle_lineb_solve(
      &cfg, static_cast<int64_t>(band.size()),
      static_cast<int64_t>(freeAtoms.size()), positions.data(), masses.data(),
      o.molecular ? 1 : 0, o.energy, o.max_steps, o.beads, surfaceCallback,
      &ctx, optimized.data(), &result);
#endif

  if (rc != RGSADDLE_OK) {
    std::string message =
        std::string("lineb: rgsaddle_lineb_solve: ") + rgsaddle_status_name(rc);
    const char *detail = rgsaddle_last_error();
    if (detail && *detail) {
      message += std::string(" (") + detail + ")";
    }
    throw std::runtime_error(message);
  }

#ifdef RGPOT_HAS_CACHE
  const auto cacheAtOutput = cacheCounts();
#endif
  const auto solveEnd = pot->forceCallCounter.load();
  unpackPositions(optimized, freeAtoms, band);

  const std::string pathFile = "lineb.con";
  const std::string resultsFile = "results.dat";
  std::vector<io::ConFrameMetadata> metas(band.size());
  for (size_t i = 0; i < band.size(); ++i) {
    metas[i].energy = band[i]->getPotentialEnergy();
  }
#ifdef WITH_RGSADDLE_LINEB_COUNTS
  std::map<std::string, int64_t> solverCalls{
      {"seed", counts.seed},
      {"optimization", counts.optimizer},
      {"path_preparation", counts.preparation},
      {"quadrature", counts.quadrature},
      {"dynamics", counts.dynamics},
      {"rate", counts.rate}};
#else
  std::map<std::string, int64_t> solverCalls{
      {"optimizer_reported", result.optimizer_calls},
      {"measure_reported", result.measure_calls},
      {"dynamics", result.dynamics_calls},
      {"rate", result.rate_calls}};
#endif
  int64_t classified = 0;
  for (const auto &entry : solverCalls)
    classified += entry.second;
#ifdef WITH_RGSADDLE_LINEB_COUNTS
  if (classified != ctx.surfaceRows)
    throw std::runtime_error("lineb: solver counts disagree with observed "
                             "configuration evaluations");
#else
  solverCalls.emplace("unclassified", ctx.surfaceRows - classified);
#endif
  solverCalls.emplace("total", ctx.surfaceRows);
  const auto outputEnd = pot->forceCallCounter.load();
  const std::map<std::string, size_t> potentialCalls{
      {"setup", solveStart - jobStart},
      {"solve", solveEnd - solveStart},
      {"output", outputEnd - solveEnd},
      {"total", outputEnd - jobStart}};
  const auto provenance = provenanceForJob(params);
  nlohmann::json metadata = {
      {"evaluation_counts",
       {{"scope", "job"},
        {"boundary", "solver_callback_geometry"},
        {"surface_calls", solverCalls},
        {"potential_calls", potentialCalls},
        {"driver_rows_by_provider",
         {{provenance.rgpot_name, outputEnd - jobStart}}}}},
      {"lineb",
       {{"scope", "job"},
        {"energy", o.energy},
        {"beads", o.beads},
        {"iterations", result.iterations},
        {"action", result.action},
        {"temperature", result.temperature},
        {"temperature_dynamics", result.temperature_dynamics},
        {"half_period", result.half_period},
        {"half_period_dynamics", result.half_period_dynamics},
        {"log_rate_over_reactant",
         o.beads ? nlohmann::json(result.log_rate_over_reactant)
                 : nlohmann::json(nullptr)},
        {"optimizer", o.projected_verlet ? "projected_verlet" : "fire"},
        {"surface_call_unit", "configuration"},
        {"surface_calls", solverCalls},
        {"potential_calls", potentialCalls}}}};
  auto &spectrum = metadata["lineb"]["spectrum"];
  spectrum = nullptr;
#ifdef RGSADDLE_HAS_LINEB_RATE_DIAGNOSTICS
  if (negativeModes >= 0 && retainedModes >= 0) {
    spectrum = {{"negative_modes", negativeModes}, {"modes", retainedModes}};
  }
#endif
  auto &providerResults = metadata["evaluation_counts"]["rgpot_results"];
  providerResults = nullptr;
#ifdef RGPOT_HAS_CACHE
  if (resultCache) {
    const auto end = cacheCounts();
    const auto delta = [](auto first, auto last) {
      return nlohmann::json{{"computed", last.computed - first.computed},
                            {"served", last.served - first.served}};
    };
    const rgpot::cache::EvalCounts zero{};
    providerResults = delta(zero, end);
    providerResults["scope"] = "job";
    providerResults["phases"] = {{"setup", delta(zero, cacheAtSolve)},
                                 {"solve", delta(cacheAtSolve, cacheAtOutput)},
                                 {"output", delta(cacheAtOutput, end)}};
  }
#endif
  std::vector<std::string> files{resultsFile, pathFile};
#ifdef RGSADDLE_HAS_LINEB_ORBIT
  if (o.beads > 0) {
    const auto component = static_cast<size_t>(orbit.max_coordinate);
    const long atom = freeAtoms.at(component / 3);
    metadata["lineb"]["stationarity"] = {
        {"discretization", "ring_polymer"},
        {"max_abs_eV_per_angstrom", orbit.max_abs_residual},
        {"rms_eV_per_angstrom", orbit.rms_residual},
        {"max_bead", orbit.max_bead},
        {"max_atom_id", band.front()->getAtomIndex(atom)},
        {"max_axis", component % 3}};
    files.emplace_back("lineb-orbit.con");
    writeRateOrbit(files.back(), *band.front(), freeAtoms, orbitPositions,
                   orbitGradients, metadata);
  }
#endif
  const auto metadataJson = metadata.dump();
  for (size_t i = 0; i < metas.size(); ++i) {
    metas[i].neb_bead = i;
    metas[i].raw_json = metadataJson;
  }
  if (!io::io_ok(io::writeNebPath(pathFile, band, metas))) {
    throw std::runtime_error("lineb: cannot write " + pathFile);
  }

  auto env = JobResultEnvelope::fromMinimization(
      RunStatus::GOOD, params.potential_options().potential,
      pot->forceCallCounter.load(), true,
      band[band.size() / 2]->getPotentialEnergy());
  env.provenance = provenance;
  env.job_type = "lineb";
  env.extras.emplace_back("action", result.action);
  env.extras.emplace_back("half_period", result.half_period);
  env.extras.emplace_back("half_period_dynamics", result.half_period_dynamics);
  env.extras.emplace_back("temperature", result.temperature);
  env.extras.emplace_back("temperature_dynamics", result.temperature_dynamics);
  env.extras.emplace_back("log_rate_over_reactant",
                          result.log_rate_over_reactant);
  env.extras.emplace_back("iterations", static_cast<double>(result.iterations));
  env.extras.emplace_back("optimizer_calls",
                          static_cast<double>(result.optimizer_calls));
  env.extras.emplace_back("measure_calls",
                          static_cast<double>(result.measure_calls));
  env.extras.emplace_back("dynamics_calls",
                          static_cast<double>(result.dynamics_calls));
  env.extras.emplace_back("rate_calls", static_cast<double>(result.rate_calls));
  for (const auto &[phase, count] : solverCalls) {
    env.extras.emplace_back("surface_calls_" + phase,
                            static_cast<double>(count));
  }
  for (const auto &[phase, count] : potentialCalls) {
    env.extras.emplace_back("potential_calls_" + phase,
                            static_cast<double>(count));
  }
  env.extras.emplace_back("energy", o.energy);
  env.extras.emplace_back("beads", static_cast<double>(o.beads));
  env.extras.emplace_back("molecular", o.molecular ? 1.0 : 0.0);
  env.extras.emplace_back("n_images", static_cast<double>(band.size()));
  env.extras.emplace_back("n_free_atoms",
                          static_cast<double>(freeAtoms.size()));
  env.writeResultsDat(resultsFile);

  EONC_LOG_INFO(
      "[Lineb] action {:.6g}, T {:.4g} K (dynamics {:.4g} K), {} steps",
      result.action, result.temperature, result.temperature_dynamics,
      result.iterations);

  return files;
#else
  return {};
#endif
}

} // namespace eonc

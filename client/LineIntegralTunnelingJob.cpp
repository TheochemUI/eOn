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
#include "eon/LineIntegralTunnelingJob.h"

#include "eon/ConFileIO.h"
#include "eon/EonLogger.h"
#include "eon/JobResult.h"
#include "eon/LineIntegralTunneling.h"
#include "eon/PotRegistry.h"
#include "eon/TunnelingJobSupport.h"

#include "magic_enum/magic_enum.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace eonc {
namespace {

struct ForceCalls {
  size_t liNeb{0};
  size_t continuation{0};
  long hessianProducts{0};
};

std::string pathFilename(long index) {
  std::ostringstream name;
  name << "line_integral_path_" << std::setfill('0') << std::setw(3) << index
       << ".con";
  return name.str();
}

void writePath(const std::string &filename,
               const tunneling::LineIntegralPath &path,
               const tunneling::LineIntegralOrbit &orbit,
               const tunneling::RateInstanton &rate, const MassWeighted &mw,
               const Matter &reference, const ForceCalls &calls) {
  Matter frame(reference);
  const double logRateSeconds =
      rate.logRate - std::log(tunneling::kTimeUnitSeconds);
  for (size_t i = 0; i < path.images.size(); ++i) {
    mw.place(path.images[i], frame);
    io::ConFrameMetadata metadata;
    metadata.frame_index = static_cast<uint64_t>(i);
    metadata.energy = path.energies[i];
    metadata.write_con_forces = false;
    if (i == 0) {
      metadata.scalars = {
          {"line_integral_energy", path.energy},
          {"line_integral_action", path.action},
          {"line_integral_temperature_K", orbit.temperature},
          {"line_integral_period", orbit.period},
          {"line_integral_period_error", orbit.periodError},
          {"line_integral_log_rate_internal", rate.logRate},
          {"line_integral_log_rate_per_s", logRateSeconds},
          {"line_integral_force_calls_li_neb",
           static_cast<double>(calls.liNeb)},
          {"line_integral_log_rate_error", rate.logRateError},
          {"line_integral_hessian_vector_products",
           static_cast<double>(calls.hessianProducts)},
          {"line_integral_force_calls_continuation",
           static_cast<double>(calls.continuation)},
      };
    }
    if (!io::io_ok(frame.matter2con(filename, i > 0, &metadata))) {
      throw std::runtime_error("line-integral tunneling cannot write " +
                               filename);
    }
  }
}

OptType optimizerType(const std::string &name) {
  const auto parsed =
      magic_enum::enum_cast<OptType>(name, magic_enum::case_insensitive);
  if (!parsed || *parsed == OptType::None || *parsed == OptType::Unknown) {
    throw std::invalid_argument(
        "line-integral tunneling optimizer is not an eOn optimizer: " + name);
  }
  return *parsed;
}

} // namespace

std::vector<std::string> LineIntegralTunnelingJob::run() {
  const auto &o = params.line_integral_tunneling_options();
  const std::string tableFilename = "line_integral_tunneling.dat";
  const std::string resultsFilename = "results.dat";
  std::vector<std::string> returnFiles{resultsFilename, tableFilename};

  Matter reactant(pot, params);
  if (!io::io_ok(reactant.con2matter(o.reactant_filename))) {
    throw std::runtime_error("line-integral tunneling cannot read " +
                             o.reactant_filename);
  }
  Matter saddle(pot, params);
  if (!io::io_ok(saddle.con2matter(o.saddle_filename))) {
    throw std::runtime_error("line-integral tunneling cannot read " +
                             o.saddle_filename);
  }
  if (reactant.numberOfAtoms() != saddle.numberOfAtoms()) {
    throw std::runtime_error(
        "line-integral tunneling reactant and saddle atom counts differ");
  }
  alignRigid(reactant, saddle);

  TunnelingJobContext context(
      params, pot, reactant,
      {.requireWholeAtomConstraints = true, .forwardNonPeriodicCell = true});
  const MassWeighted &mw = context.coordinates();
  const VectorXd qReactant = VectorXd::Zero(mw.dimension());
  const VectorXd qSaddle = mw.toQ(saddle);
  auto potential = [&](const std::vector<VectorXd> &q,
                       std::vector<double> &energies,
                       std::vector<VectorXd> &gradients) {
    context.evaluate(q, energies, gradients);
  };
  auto align = [&](std::vector<VectorXd> &path) {
    Matter image(reactant);
    Matter previous(reactant);
    for (VectorXd &q : path) {
      mw.place(q, image);
      alignRigid(previous, image);
      q = mw.toQ(image);
      previous = image;
    }
  };

  std::vector<VectorXd> mep;
  const auto frames = readcon::read_all_frames(o.mep_filename);
  Matter previousMep(reactant);
  for (const auto &input : frames) {
    Matter image(reactant);
    if (!io::io_ok(io::con2matter(image, input))) {
      throw std::runtime_error("line-integral tunneling cannot read " +
                               o.mep_filename);
    }
    if (image.numberOfAtoms() != reactant.numberOfAtoms()) {
      throw std::runtime_error(
          "line-integral tunneling MEP atom count differs");
    }
    alignRigid(previousMep, image);
    mep.push_back(mw.toQ(image));
    previousMep = image;
  }
  if (mep.size() < 3) {
    throw std::runtime_error(
        "line-integral tunneling MEP has fewer than three images");
  }

  std::vector<double> mepEnergies;
  std::vector<VectorXd> mepGradients;
  context.evaluate(mep, mepEnergies, mepGradients);
  std::vector<double> stationaryEnergies;
  std::vector<VectorXd> stationaryGradients;
  context.evaluate({qReactant, qSaddle}, stationaryEnergies,
                   stationaryGradients);
  const double reactantEnergy = stationaryEnergies[0];
  const double saddleEnergy = stationaryEnergies[1];
  if (!std::isfinite(reactantEnergy) || !std::isfinite(saddleEnergy) ||
      !std::all_of(mepEnergies.begin(), mepEnergies.end(),
                   [](double value) { return std::isfinite(value); })) {
    throw std::runtime_error(
        "line-integral tunneling received a non-finite potential energy");
  }
  // The saddle's unstable mode starts from the MEP tangent through its
  // highest image.
  const auto highest = static_cast<size_t>(
      std::max_element(mepEnergies.begin(), mepEnergies.end()) -
      mepEnergies.begin());
  const VectorXd tangent = mep[std::min(highest + 1, mep.size() - 1)] -
                           mep[highest > 0 ? highest - 1 : 0];
  const TunnelingJobContext::LowestMode unstable =
      context.lowestModeAt(qSaddle, tangent);
  if (!(unstable.curvature < 0.0)) {
    throw std::runtime_error(
        "line-integral tunneling saddle has no unstable mode");
  }
  // T_c = hbar omega_b / (2 pi kB), the period at which the harmonic
  // barrier's imaginary-time oscillation first fits. The Ritz value is
  // within its residual of the eigenvalue, so orbits are checked against
  // the longest crossover period that allows.
  const double crossoverTemperature =
      tunneling::kHbar * std::sqrt(-unstable.curvature) /
      (2.0 * std::numbers::pi * tunneling::kBoltzmann);
  const double weakestCurvature = -unstable.curvature - unstable.residual;
  if (!(weakestCurvature > 0.0)) {
    throw std::runtime_error(
        "line-integral tunneling saddle mode residual exceeds its curvature");
  }
  const double crossoverPeriod =
      2.0 * std::numbers::pi / std::sqrt(weakestCurvature);
  context.findRotationZeroModes(qReactant);
  const MatrixXd reactantRigid =
      mw.rigidBasis(reactant, context.rotationZero());
  const tunneling::HessianProduct reactantProducts =
      context.hessianProductAt(qReactant);
  tunneling::RingRigidBodies rigidBodies;
  if (static_cast<long>(mw.sqrtMasses().size()) == reactant.numberOfAtoms()) {
    rigidBodies.sqrtMasses = mw.sqrtMasses();
    rigidBodies.reference = mw.referenceFree();
    rigidBodies.rotations = context.rotationZero();
  }

  tunneling::LineIntegralOptions options;
  options.images = o.images;
  options.spring = o.spring;
  options.contour = o.contour;
  options.forceTolerance = o.force_tolerance;
  options.maxIterations = o.max_iterations;
  options.timeStep = o.li_neb_time_step;

  std::ofstream table(tableFilename);
  if (!table) {
    throw std::runtime_error("line-integral tunneling cannot write " +
                             tableFilename);
  }
  table << std::setprecision(17)
        << "# index energy_eV action temperature_K period period_error "
           "ln_k_internal "
           "ln_k_per_s ln_k_error rate_per_s force_calls_li_neb "
           "force_calls_continuation hessian_vector_products "
           "interior_gradient endpoint_gradient negative_modes\n";

  size_t totalLiNeb = 0;
  size_t totalContinuation = 0;
  tunneling::LineIntegralPath previous;
  tunneling::RateInstanton lastRate;
  tunneling::LineIntegralOrbit lastOrbit;
  std::vector<tunneling::LineIntegralContinuationPoint> continuation;
  bool determinantEstimatedErrorTargetSatisfied = true;
  long level = 0;
  double energy = saddleEnergy - o.energy_offset;
  while (true) {
    std::vector<VectorXd> seed =
        level == 0
            ? tunneling::lineIntegralSeed(mep, mepEnergies, energy, o.images)
            : tunneling::lineIntegralContinuation(mep, mepEnergies, previous,
                                                  energy, o.images);

    const size_t beforeBand = PotRegistry::get().total_force_calls();
    tunneling::LineIntegralPath path;
    if (o.optimizer == "velocity_verlet") {
      path = tunneling::optimizeLineIntegralPath(std::move(seed), energy,
                                                 potential, options, align);
    } else {
      path = tunneling::optimizeLineIntegralPath(
          std::move(seed), energy, potential, options, params,
          optimizerType(o.optimizer), align);
    }
    if (!path.converged) {
      throw std::runtime_error(
          "line-integral tunneling path did not converge at energy " +
          std::to_string(energy));
    }
    const size_t afterBand = PotRegistry::get().total_force_calls();
    ForceCalls calls;
    if (level == 0) {
      calls.liNeb = afterBand - beforeBand;
      totalLiNeb += calls.liNeb;
    } else {
      calls.continuation = afterBand - beforeBand;
      totalContinuation += calls.continuation;
    }

    const tunneling::LineIntegralOrbit orbit = tunneling::lineIntegralOrbit(
        path, o.replicas, o.period_tolerance, crossoverPeriod);

    const long beforeProducts = context.hessianVectorProducts();
    const tunneling::RateInstanton rate = tunneling::hqtstRate(
        orbit,
        {.ring = context.ringHessianProducts(orbit.replicas),
         .reactant = reactantProducts},
        reactantEnergy, reactantRigid,
        {.estimatedLogRateErrorTarget = o.determinant_estimated_error,
         .seed = static_cast<std::uint64_t>(level)},
        rigidBodies);
    determinantEstimatedErrorTargetSatisfied &= rate.converged;
    calls.hessianProducts = context.hessianVectorProducts() - beforeProducts;
    if (rate.negativeModes != 1) {
      throw std::runtime_error(
          "line-integral tunneling ring Hessian does not have one negative "
          "mode");
    }
    const double logRateSeconds =
        rate.logRate - std::log(tunneling::kTimeUnitSeconds);
    const std::string pathFile = pathFilename(level);
    writePath(pathFile, path, orbit, rate, mw, reactant, calls);
    returnFiles.push_back(pathFile);

    table << level << ' ' << energy << ' ' << path.action << ' '
          << orbit.temperature << ' ' << orbit.period << ' '
          << orbit.periodError << ' ' << rate.logRate << ' ' << logRateSeconds
          << ' ' << rate.logRateError << ' ' << rate.rate << ' ' << calls.liNeb
          << ' ' << calls.continuation << ' ' << calls.hessianProducts << ' '
          << path.interiorGradient << ' ' << path.endpointGradient << ' '
          << rate.negativeModes << '\n';
    EONC_LOG_INFO("[Line Integral Tunneling] E={:.8g} eV W={:.8g} T={:.8g} K "
                  "ln(k/s^-1)={:.8g}+-{:.3g}; force calls LI-NEB={}, "
                  "continuation={}; Hessian-vector products={}",
                  energy, path.action, orbit.temperature, logRateSeconds,
                  rate.logRateError, calls.liNeb, calls.continuation,
                  calls.hessianProducts);
    previous = path;
    lastRate = rate;
    lastOrbit = orbit;
    continuation.push_back({energy, orbit.temperature, rate.logRate});
    if (orbit.temperature <= o.lowest_temperature) {
      break;
    }
    energy = tunneling::nextLineIntegralEnergy(continuation, o.energy_step,
                                               o.log_rate_tolerance,
                                               o.lowest_temperature);
    ++level;
  }
  table.close();

  const RunStatus status = determinantEstimatedErrorTargetSatisfied
                               ? RunStatus::GOOD
                               : RunStatus::FAIL_ESTIMATED_ERROR_TARGET;
  auto envelope = JobResultEnvelope::fromMinimization(
      status, params.potential_options().potential,
      PotRegistry::get().total_force_calls(), false, 0.0);
  envelope.job_type = "line_integral_tunneling";
  envelope.extras = {
      {"line_integral_levels", static_cast<double>(continuation.size())},
      {"line_integral_energy_eV", previous.energy},
      {"line_integral_action", previous.action},
      {"line_integral_temperature_K", lastOrbit.temperature},
      {"line_integral_crossover_temperature_K", crossoverTemperature},
      {"line_integral_period", lastOrbit.period},
      {"line_integral_period_error", lastOrbit.periodError},
      {"line_integral_log_rate_internal", lastRate.logRate},
      {"line_integral_log_rate_per_s",
       lastRate.logRate - std::log(tunneling::kTimeUnitSeconds)},
      {"line_integral_rate_per_s", lastRate.rate},
      {"line_integral_log_rate_error", lastRate.logRateError},
      {"line_integral_determinant_estimated_error_target_satisfied",
       determinantEstimatedErrorTargetSatisfied ? 1.0 : 0.0},
      {"line_integral_finite_difference_hessian_fallback",
       context.usesFiniteDifferenceHessianFallback() ? 1.0 : 0.0},
      {"force_calls_li_neb", static_cast<double>(totalLiNeb)},
      {"force_calls_continuation", static_cast<double>(totalContinuation)},
      {"hessian_vector_products",
       static_cast<double>(context.hessianVectorProducts())},
  };
  envelope.writeResultsDat(resultsFilename);
  return returnFiles;
}

} // namespace eonc

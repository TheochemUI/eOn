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
#include "eon/LineIntegralTunneling.h"

#include "catch2/catch_amalgamated.hpp"
#include "eon/ConFileIO.h"
#include "eon/LineIntegralTunnelingJob.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"
#include "eon/TunnelingJobSupport.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <sstream>
#include <thread>

using namespace eonc;
using namespace eonc::tunneling;
using Catch::Matchers::WithinAbs;

namespace {

class NonFiniteEnergyPotential final : public Potential {
public:
  explicit NonFiniteEnergyPotential(const Parameters &parameters)
      : Potential(PotType::LJCLUSTER, parameters) {}

  void force(long nAtoms, const double *, const int *, double *forces,
             double *energy, double *variance, const double *) override {
    std::fill(forces, forces + 3 * nAtoms, 0.0);
    *energy = std::numeric_limits<double>::quiet_NaN();
    if (variance != nullptr) {
      *variance = 0.0;
    }
  }
};

class CurrentPathGuard {
public:
  CurrentPathGuard()
      : original_(std::filesystem::current_path()) {}
  ~CurrentPathGuard() { std::filesystem::current_path(original_); }

private:
  std::filesystem::path original_;
};

double logPartitionRatio(const RateInstanton &rate) {
  const double bnh = rate.betaN * kHbar;
  const double zeroModeFactor =
      0.5 *
      std::log(rate.bN / (2.0 * std::numbers::pi * rate.betaN * kHbar * kHbar));
  return rate.logRate + std::log(bnh) - zeroModeFactor +
         rate.betaN * rate.ringPotential;
}

/// The documentation tables and the numbers quoted in their prose are both
/// generated from these files, so they are written at full precision.
void writeTable(const std::filesystem::path &file,
                const std::vector<std::string> &header,
                const std::vector<std::vector<double>> &rows) {
  std::ofstream out(file);
  REQUIRE(out.good());
  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  for (size_t column = 0; column < header.size(); ++column) {
    out << (column == 0 ? "" : ",") << header[column];
  }
  out << '\n';
  for (const auto &row : rows) {
    REQUIRE(row.size() == header.size());
    for (size_t column = 0; column < row.size(); ++column) {
      out << (column == 0 ? "" : ",") << row[column];
    }
    out << '\n';
  }
  REQUIRE(out.good());
}

std::vector<VectorXd> loadMep(const Matter &reactant,
                              const MassWeighted &coordinates) {
  std::vector<VectorXd> mep;
  std::ifstream input("cuh2_mep_h.dat");
  REQUIRE(input.good());
  std::array<double, 6> hydrogen{};
  while (input >> hydrogen.front()) {
    for (size_t i = 1; i < hydrogen.size(); ++i) {
      REQUIRE(input >> hydrogen[i]);
    }
    AtomMatrix positions = reactant.getPositions();
    const long mobileAtoms = static_cast<long>(hydrogen.size() / 3);
    const long firstMobileAtom = reactant.numberOfAtoms() - mobileAtoms;
    REQUIRE(firstMobileAtom >= 0);
    for (long atom = firstMobileAtom; atom < reactant.numberOfAtoms(); ++atom) {
      for (long axis = 0; axis < 3; ++axis) {
        positions(atom, axis) =
            hydrogen[static_cast<size_t>((atom - firstMobileAtom) * 3 + axis)];
      }
    }
    Matter image(reactant);
    image.setPositions(positions);
    mep.push_back(coordinates.toQ(image));
  }
  REQUIRE(mep.size() == 22);
  std::reverse(mep.begin(), mep.end());
  return mep;
}

void writeJobInputs(Matter &reactant, const MassWeighted &coordinates,
                    const std::vector<VectorXd> &mep, size_t saddle) {
  REQUIRE(io::io_ok(reactant.matter2con("reactant.con")));
  Matter saddleMatter(reactant);
  coordinates.place(mep[saddle], saddleMatter);
  REQUIRE(io::io_ok(saddleMatter.matter2con("saddle.con")));
  for (size_t i = 0; i < mep.size(); ++i) {
    Matter image(reactant);
    coordinates.place(mep[i], image);
    REQUIRE(io::io_ok(image.matter2con("neb.con", i > 0)));
  }
}

std::vector<std::vector<double>>
readNumericTable(const std::filesystem::path &file) {
  std::ifstream input(file);
  REQUIRE(input.good());
  std::string line;
  REQUIRE(static_cast<bool>(std::getline(input, line)));
  std::vector<std::vector<double>> rows;
  while (std::getline(input, line)) {
    std::istringstream stream(line);
    std::vector<double> row;
    double value = 0.0;
    while (stream >> value) {
      row.push_back(value);
    }
    REQUIRE_FALSE(row.empty());
    rows.push_back(std::move(row));
  }
  return rows;
}

std::map<std::string, std::string>
readResults(const std::filesystem::path &file) {
  std::ifstream input(file);
  REQUIRE(input.good());
  std::map<std::string, std::string> values;
  std::string value;
  std::string key;
  while (input >> value >> key) {
    values.emplace(key, value);
  }
  return values;
}

constexpr int documentedRatioDecimalPlaces = 4;

std::vector<double> documentedRefinedRatios() {
  const char *source = std::getenv("EON_SOURCE_DIR");
  REQUIRE(source != nullptr);
  std::ifstream input(std::filesystem::path(source) / "docs/source/fig/data" /
                      "cuh2_line_integral_convergence.csv");
  REQUIRE(input.good());
  std::string line;
  REQUIRE(static_cast<bool>(std::getline(input, line)));
  std::vector<double> ratios;
  while (std::getline(input, line)) {
    std::istringstream row(line);
    std::vector<double> columns;
    std::string value;
    while (std::getline(row, value, ',')) {
      columns.push_back(std::stod(value));
    }
    REQUIRE(columns.size() == 10);
    const double scale =
        std::pow(10.0, static_cast<double>(documentedRatioDecimalPlaces));
    ratios.push_back(
        std::round(std::exp(columns[7] - columns[9]) * scale) / scale);
  }
  REQUIRE_FALSE(ratios.empty());
  return ratios;
}

std::vector<std::array<double, 2>> documentedJobRates() {
  const char *source = std::getenv("EON_SOURCE_DIR");
  REQUIRE(source != nullptr);
  std::ifstream input(std::filesystem::path(source) / "docs/source/fig/data" /
                      "cuh2_line_integral_7level.csv");
  REQUIRE(input.good());
  std::string line;
  REQUIRE(static_cast<bool>(std::getline(input, line)));
  std::vector<std::array<double, 2>> rates;
  while (std::getline(input, line)) {
    std::istringstream row(line);
    std::string value;
    std::vector<double> columns;
    while (std::getline(row, value, ',')) {
      columns.push_back(std::stod(value));
    }
    REQUIRE(columns.size() == 10);
    rates.push_back(
        {columns[1], columns[2] + std::log(kTimeUnitSeconds)});
  }
  REQUIRE(rates.size() >= 2);
  return rates;
}

} // namespace

TEST_CASE("Cu(110) plus two H has an LI-NEB tunnelling path",
          "[LineIntegralTunneling][CuH2]") {
  Parameters parameters;
  ParametersLoadAccess::potential_options(parameters).potential = PotType::CUH2;
  auto potential = helpers::sharePotential(
      helpers::makePotential(PotType::CUH2, parameters));
  Matter reactant(potential, parameters);
  REQUIRE(io::io_ok(reactant.con2matter("min2_0Cu.con")));

  Matter partiallyConstrained(reactant);
  long firstMobileAtom = -1;
  for (long atom = 0; atom < reactant.numberOfAtoms(); ++atom) {
    if (!reactant.getFixed(atom)) {
      firstMobileAtom = atom;
      break;
    }
  }
  REQUIRE(firstMobileAtom >= 0);
  partiallyConstrained.setFixed(firstMobileAtom, 0, 1);
  REQUIRE_THROWS_AS(TunnelingJobContext(parameters, potential,
                                        partiallyConstrained,
                                        {.requireWholeAtomConstraints = true,
                                         .forwardNonPeriodicCell = true}),
                    std::invalid_argument);

  const unsigned detectedThreads = std::thread::hardware_concurrency();
  const long hessianProductThreads =
      static_cast<long>(detectedThreads == 0 ? 1 : detectedThreads);
  INFO("Hessian product threads: " << hessianProductThreads);
  TunnelingJobContext context(
      parameters, potential, reactant,
      {.requireWholeAtomConstraints = true,
       .forwardNonPeriodicCell = true,
       .hessianProductThreads = hessianProductThreads});
  const MassWeighted &coordinates = context.coordinates();
  REQUIRE(potential->supportsMassWeightedHessianVectorProduct());
  REQUIRE_FALSE(context.usesFiniteDifferenceHessianFallback());
  const std::vector<VectorXd> mep = loadMep(reactant, coordinates);
  std::vector<double> energies;
  std::vector<VectorXd> gradients;
  context.evaluate(mep, energies, gradients);
  const size_t saddle = static_cast<size_t>(
      std::max_element(energies.begin(), energies.end()) - energies.begin());
  REQUIRE(saddle > 0);
  REQUIRE(saddle + 1 < mep.size());

  // Section II E of Ásgeirsson et al. uses 5-10 LI images, k=20 A^-2,
  // kappa=1000 eV^-1 A^-1, |g|<0.01 A^-1, and P=150 for the independently
  // optimized ring reference. Its atomic-force criterion is 10^-8 eV/A.
  constexpr long paperImages = 10;
  constexpr long paperReplicas = 150;
  constexpr double paperSpring = 20.0;
  constexpr double paperContour = 1000.0;
  constexpr double paperPathTolerance = 0.01;
  const double paperReferenceTolerance = std::pow(10.0, -8.0);
  constexpr double maximumEnergyStep = 0.01;
  constexpr double lowestTemperature = 79.0;
  constexpr double refinementTemperature = 100.0;
  constexpr double maximumLiInstantonLogRateDifference = 0.25;
  constexpr double figureReadUncertaintyLog10 = 0.5;
  constexpr long rateStride = 5;
  const double energyOffset =
      maximumEnergyStep * static_cast<double>(rateStride);
  const double logRateTolerance = std::log(10.0) * figureReadUncertaintyLog10;
  LineIntegralOptions options;
  options.images = paperImages;
  options.spring = paperSpring;
  options.contour = paperContour;
  options.forceTolerance = paperPathTolerance;

  const double determinantEstimatedError =
      ParametersLoadAccess::line_integral_tunneling_options(parameters)
          .determinant_estimated_error;

  // Rates take the products alone, as the job does. The dense Hessians
  // below, one product per column, are oracles for the crossover and the
  // reference instanton only. eOn's forward-difference Hessian is off by
  // about 2% on the H modes, which would bias the reference.
  const VectorXd reactantQ = VectorXd::Zero(coordinates.dimension());
  context.findRotationZeroModes(reactantQ);
  const MatrixXd reactantRigid =
      coordinates.rigidBasis(reactant, context.rotationZero());
  const HessianProduct reactantProducts = context.hessianProductAt(reactantQ);
  RingRigidBodies rigidBodies;
  if (static_cast<long>(coordinates.sqrtMasses().size()) ==
      reactant.numberOfAtoms()) {
    rigidBodies.sqrtMasses = coordinates.sqrtMasses();
    rigidBodies.reference = coordinates.referenceFree();
    rigidBodies.rotations = context.rotationZero();
  }
  auto lineIntegralRate = [&](const LineIntegralOrbit &candidate,
                              std::uint64_t seed) {
    return hqtstRate(candidate,
                     {.ring = context.ringHessianProducts(candidate.replicas),
                      .reactant = reactantProducts},
                     energies.front(), reactantRigid,
                     {.estimatedLogRateErrorTarget = determinantEstimatedError,
                      .seed = seed},
                     rigidBodies);
  };

  const TunnelingJobContext::LowestMode unstable =
      context.lowestModeAt(mep[saddle], mep[saddle + 1] - mep[saddle - 1]);
  REQUIRE(unstable.curvature < 0.0);
  const double weakestCurvature = -unstable.curvature - unstable.residual;
  REQUIRE(weakestCurvature > 0.0);
  const double crossoverPeriod =
      2.0 * std::numbers::pi / std::sqrt(weakestCurvature);
  const double periodTolerance =
      std::sqrt(std::numeric_limits<double>::epsilon()) * crossoverPeriod +
      paperPathTolerance / paperContour;

  auto denseHessian = [&](const VectorXd &q) {
    const HessianProduct product = context.hessianProductAt(q);
    const long freedom = coordinates.dimension();
    MatrixXd columns(freedom, freedom);
    for (long i = 0; i < freedom; ++i) {
      columns.col(i) = product(VectorXd::Unit(freedom, i));
    }
    return MatrixXd(0.5 * (columns + columns.transpose()));
  };
  const MatrixXd reactantHessian = denseHessian(reactantQ);
  const MatrixXd saddleHessian = denseHessian(mep[saddle]);
  // For a unit u, the symmetric oracle has an eigenvalue within
  // |H_o u - theta u| <= |H_o u - H u| + |H u - theta u| of theta.
  {
    const MatrixXd &oracle = saddleHessian;
    const VectorXd u = unstable.vector.normalized();
    const VectorXd hu = context.hessianProductAt(mep[saddle])(u);
    const double reach =
        (oracle * u - hu).norm() + (hu - unstable.curvature * u).norm();
    const VectorXd saddleEigenvalues =
        Eigen::SelfAdjointEigenSolver<MatrixXd>(oracle).eigenvalues();
    INFO("Lanczos curvature " << unstable.curvature << ", residual "
                              << unstable.residual << ", oracle lowest "
                              << saddleEigenvalues(0) << ", reach " << reach);
    REQUIRE(std::abs(saddleEigenvalues(0) - unstable.curvature) <= reach);
    REQUIRE(saddleEigenvalues(1) > 0.0);
  }
  const VectorXd reactantEigenvalues =
      Eigen::SelfAdjointEigenSolver<MatrixXd>(reactantHessian).eigenvalues();
  INFO(reactantEigenvalues.transpose());
  REQUIRE(reactantEigenvalues.minCoeff() > 0.0);

  RateInstantonOptions referenceOptions;
  referenceOptions.beads = paperReplicas;
  referenceOptions.maxIterations = 5000;
  referenceOptions.forceTolerance = paperReferenceTolerance;
  referenceOptions.halfRing = true;

  std::vector<size_t> forceCalls;
  std::vector<double> allTemperatures;
  std::vector<double> allLogRates;
  std::vector<double> temperatures;
  std::vector<double> logRateDifferences;
  std::vector<double> rateRatios;
  std::vector<double> partitionLogDifferences;
  std::vector<double> replicaConvergenceErrors;
  std::vector<double> refinedRateRatios;
  std::vector<double> refinedPartitionDifferences;
  std::vector<double> refinedRateErrors;
  std::vector<long> hessianProducts;
  std::vector<std::vector<double>> callRows;
  std::vector<std::vector<double>> rateRows;
  std::vector<std::vector<double>> convergenceRows;
  std::vector<LineIntegralContinuationPoint> continuation;
  LineIntegralPath previous;
  long level = 0;
  double energy = energies[saddle] - energyOffset;
  while (true) {
    const LineIntegralPath higherEnergyPath = previous;
    REQUIRE(energies.front() < energy);
    REQUIRE(energies.back() < energy);
    const std::vector<VectorXd> seed =
        level == 0 ? lineIntegralSeed(mep, energies, energy, options.images)
                   : lineIntegralContinuation(mep, energies, previous, energy,
                                              options.images);
    const size_t beforeBand = PotRegistry::get().total_force_calls();
    const LineIntegralPath path = optimizeLineIntegralPath(
        seed, energy,
        [&](const std::vector<VectorXd> &q, std::vector<double> &v,
            std::vector<VectorXd> &g) { context.evaluate(q, v, g); },
        options);
    forceCalls.push_back(PotRegistry::get().total_force_calls() - beforeBand);
    INFO("E=" << energy << " eV, LI calls=" << forceCalls.back());
    REQUIRE(path.converged);
    REQUIRE(path.action > 0.0);
    previous = path;

    const LineIntegralOrbit orbit = lineIntegralOrbit(
        path, paperReplicas, periodTolerance, crossoverPeriod);
    const auto stream = static_cast<std::uint64_t>(level);
    const RateInstanton lineIntegral = lineIntegralRate(orbit, stream);
    REQUIRE(lineIntegral.negativeModes == 1);
    REQUIRE(lineIntegral.logRateError <= determinantEstimatedError);
    allTemperatures.push_back(orbit.temperature);
    allLogRates.push_back(lineIntegral.logRate);
    hessianProducts.push_back(lineIntegral.hessianProducts);
    continuation.push_back({energy, orbit.temperature, lineIntegral.logRate});
    callRows.push_back({energy, orbit.temperature,
                        static_cast<double>(forceCalls.back()),
                        static_cast<double>(lineIntegral.hessianProducts)});
    const bool reachedLowestTemperature =
        orbit.temperature <= lowestTemperature;

    if (level % rateStride != 0 && !reachedLowestTemperature) {
      energy = nextLineIntegralEnergy(continuation, maximumEnergyStep,
                                      logRateTolerance, lowestTemperature);
      ++level;
      continue;
    }
    temperatures.push_back(orbit.temperature);

    RateInstanton reference = optimizeRateInstanton(
        mep[saddle], saddleHessian, orbit.period / kHbar, orbit.replicas,
        [&](const std::vector<VectorXd> &q, std::vector<double> &v,
            std::vector<VectorXd> &g) { context.evaluate(q, v, g); },
        referenceOptions);
    REQUIRE(reference.converged);
    instantonRate(
        reference,
        [&](long, const VectorXd &q) { return context.beadHessianAt(q); },
        reactantHessian, energies.front(), saddleHessian, energies[saddle]);
    REQUIRE(reference.negativeModes == 1);
    const double difference = lineIntegral.logRate - reference.logRate;
    REQUIRE(std::abs(difference) < maximumLiInstantonLogRateDifference);
    const double rateRatio = std::exp(difference);
    logRateDifferences.push_back(difference);
    rateRatios.push_back(rateRatio);
    rateRows.push_back({orbit.temperature, lineIntegral.logRate,
                        lineIntegral.logRateError, reference.logRate});
    partitionLogDifferences.push_back(logPartitionRatio(lineIntegral) -
                                      logPartitionRatio(reference));
    INFO("E=" << energy << " eV, T=" << orbit.temperature
              << " K, W=" << path.action << ", LI calls=" << forceCalls.back()
              << ", ln(k)_LI=" << lineIntegral.logRate << ", ln(k)_reference="
              << reference.logRate << ", k_LI/k_reference=" << rateRatio);

    // Reported points at or below this temperature separate path
    // discretization (n) and ring discretization (P) errors.
    if (orbit.temperature <= refinementTemperature) {
      auto refinedPath = [&](long images) {
        LineIntegralOptions refinedOptions = options;
        refinedOptions.images = images;
        const std::vector<VectorXd> refinedSeed = lineIntegralContinuation(
            mep, energies, higherEnergyPath, energy, images);
        return optimizeLineIntegralPath(
            refinedSeed, energy,
            [&](const std::vector<VectorXd> &q, std::vector<double> &v,
                std::vector<VectorXd> &g) { context.evaluate(q, v, g); },
            refinedOptions);
      };
      const LineIntegralPath path20 = refinedPath(20);
      const LineIntegralPath path40 = refinedPath(40);
      LineIntegralOptions options80 = options;
      options80.images = 80;
      const LineIntegralPath path80 = optimizeLineIntegralPath(
          resampleLineIntegralPath(path40.images, options80.images), energy,
          [&](const std::vector<VectorXd> &q, std::vector<double> &v,
              std::vector<VectorXd> &g) { context.evaluate(q, v, g); },
          options80);
      REQUIRE(path20.converged);
      REQUIRE(path40.converged);
      REQUIRE(path80.converged);

      auto refinedRate = [&](const LineIntegralPath &candidate, long replicas) {
        return lineIntegralRate(lineIntegralOrbit(candidate, replicas,
                                                  periodTolerance,
                                                  crossoverPeriod),
                                stream);
      };
      const RateInstanton rate20P150 = refinedRate(path20, 150);
      const RateInstanton rate10P600 = refinedRate(path, 600);
      const RateInstanton rate40P150 = refinedRate(path40, 150);
      const RateInstanton rate40P300 = refinedRate(path40, 300);
      const RateInstanton rate40P600 = refinedRate(path40, 600);
      const RateInstanton rate80P150 = refinedRate(path80, 150);
      const RateInstanton rate80P300 = refinedRate(path80, 300);
      const RateInstanton rate80P600 = refinedRate(path80, 600);

      auto referenceRate = [&](const LineIntegralPath &candidate,
                               long replicas) {
        const LineIntegralOrbit candidateOrbit = lineIntegralOrbit(
            candidate, replicas, periodTolerance, crossoverPeriod);
        RateInstantonOptions candidateOptions = referenceOptions;
        candidateOptions.beads = replicas;
        RateInstanton result = optimizeRateInstanton(
            mep[saddle], saddleHessian, candidateOrbit.period / kHbar,
            candidateOrbit.replicas,
            [&](const std::vector<VectorXd> &q, std::vector<double> &v,
                std::vector<VectorXd> &g) { context.evaluate(q, v, g); },
            candidateOptions);
        REQUIRE(result.converged);
        instantonRate(
            result,
            [&](long, const VectorXd &q) { return context.beadHessianAt(q); },
            reactantHessian, energies.front(), saddleHessian, energies[saddle]);
        return result;
      };
      const RateInstanton reference80P150 = referenceRate(path80, 150);
      const RateInstanton reference80P300 = referenceRate(path80, 300);
      const RateInstanton reference80P600 = referenceRate(path80, 600);
      const LineIntegralOrbit orbit80P600 =
          lineIntegralOrbit(path80, 600, periodTolerance, crossoverPeriod);
      INFO("convergence T="
           << orbit.temperature << " K: ln(k) n10/P150=" << lineIntegral.logRate
           << ", n20/P150=" << rate20P150.logRate << ", n10/P600="
           << rate10P600.logRate << ", n40/P150=" << rate40P150.logRate
           << ", n40/P300=" << rate40P300.logRate << ", n40/P600="
           << rate40P600.logRate << ", n80/P150=" << rate80P150.logRate
           << ", n80/P600=" << rate80P600.logRate << " +- "
           << rate80P600.logRateError
           << ", n80 temperature=" << orbit80P600.temperature
           << " K, matched instanton P150=" << reference80P150.logRate
           << ", P300=" << reference80P300.logRate
           << ", matched P600=" << reference80P600.logRate);
      for (const RateInstanton *rate :
           {&rate20P150, &rate10P600, &rate40P150, &rate40P300, &rate40P600,
            &rate80P150, &rate80P300, &rate80P600}) {
        REQUIRE(rate->negativeModes == 1);
        REQUIRE(rate->logRateError <= determinantEstimatedError);
      }
      // Halving the ring spacing shrinks the change in ln k, up to the
      // stochastic error bars of the three rates involved.
      auto bars = [](const RateInstanton &a, const RateInstanton &b) {
        return a.logRateError + b.logRateError;
      };
      const double replicaShift =
          std::abs(rate10P600.logRate - lineIntegral.logRate);
      const double refinedRateError =
          std::abs(rate80P600.logRate - rate80P300.logRate) +
          std::abs(reference80P600.logRate - reference80P300.logRate) +
          std::abs(rate80P600.logRate - rate40P600.logRate) +
          rate80P600.logRateError;
      REQUIRE(std::abs(rate80P600.logRate - rate80P300.logRate) <
              std::abs(rate80P300.logRate - rate80P150.logRate) +
                  bars(rate80P600, rate80P300) + bars(rate80P300, rate80P150));
      REQUIRE(std::abs(reference80P600.logRate - reference80P300.logRate) <
              std::abs(reference80P300.logRate - reference80P150.logRate));
      REQUIRE(std::abs(logPartitionRatio(rate80P600) -
                       logPartitionRatio(rate80P300)) <
              std::abs(logPartitionRatio(rate80P300) -
                       logPartitionRatio(rate80P150)) +
                  bars(rate80P600, rate80P300) + bars(rate80P300, rate80P150));
      REQUIRE(std::abs(logPartitionRatio(reference80P600) -
                       logPartitionRatio(reference80P300)) <
              std::abs(logPartitionRatio(reference80P300) -
                       logPartitionRatio(reference80P150)));
      replicaConvergenceErrors.push_back(replicaShift);
      refinedRateRatios.push_back(
          std::exp(rate80P600.logRate - reference80P600.logRate));
      refinedPartitionDifferences.push_back(logPartitionRatio(rate80P600) -
                                            logPartitionRatio(reference80P600));
      refinedRateErrors.push_back(refinedRateError);
      convergenceRows.push_back(
          {orbit.temperature, lineIntegral.logRate, rate20P150.logRate,
           rate40P150.logRate, rate80P150.logRate, rate10P600.logRate,
           orbit80P600.temperature, rate80P600.logRate, rate80P600.logRateError,
           reference80P600.logRate});
    }
    if (reachedLowestTemperature) {
      break;
    }
    energy = nextLineIntegralEnergy(continuation, maximumEnergyStep,
                                    logRateTolerance, lowestTemperature);
    ++level;
  }
  const std::filesystem::path tableDirectory =
      std::filesystem::temp_directory_path() / "eon_cuh2_line_integral";
  std::filesystem::remove_all(tableDirectory);
  std::filesystem::create_directories(tableDirectory);
  writeTable(
      tableDirectory / "cuh2_line_integral_calls.csv",
      {"energy_eV", "temperature_K", "force_calls", "hessian_vector_products"},
      callRows);
  writeTable(tableDirectory / "cuh2_line_integral_rates.csv",
             {"temperature_K", "log_rate_li", "log_rate_li_error",
              "log_rate_instanton"},
             rateRows);
  writeTable(tableDirectory / "cuh2_line_integral_convergence.csv",
             {"temperature_n10_K", "log_rate_n10_P150", "log_rate_n20_P150",
              "log_rate_n40_P150", "log_rate_n80_P150", "log_rate_n10_P600",
              "temperature_n80_K", "log_rate_n80_P600",
              "log_rate_n80_P600_error", "log_rate_instanton_P600"},
             convergenceRows);
  const double averageCalls =
      static_cast<double>(
          std::accumulate(forceCalls.begin(), forceCalls.end(), size_t{0})) /
      static_cast<double>(forceCalls.size());
  CAPTURE(forceCalls, averageCalls, temperatures, logRateDifferences,
          rateRatios, partitionLogDifferences, replicaConvergenceErrors,
          refinedRateRatios, refinedPartitionDifferences, refinedRateErrors,
          hessianProducts);
  REQUIRE(rateRatios.size() == temperatures.size());
  REQUIRE(partitionLogDifferences.size() == rateRatios.size());
  for (size_t i = 0; i < rateRatios.size(); ++i) {
    REQUIRE(std::isfinite(rateRatios[i]));
    REQUIRE(std::isfinite(partitionLogDifferences[i]));
  }
  const size_t expectedRefinementCount = static_cast<size_t>(std::count_if(
      temperatures.begin(), temperatures.end(), [](double temperature) {
        return temperature <= refinementTemperature;
      }));
  REQUIRE(expectedRefinementCount > 0);
  REQUIRE(replicaConvergenceErrors.size() == expectedRefinementCount);
  REQUIRE(convergenceRows.size() == expectedRefinementCount);
  REQUIRE(refinedRateRatios.size() == replicaConvergenceErrors.size());
  REQUIRE(refinedPartitionDifferences.size() ==
          replicaConvergenceErrors.size());
  REQUIRE(refinedRateErrors.size() == refinedRateRatios.size());
  const std::vector<double> documentedRatios = documentedRefinedRatios();
  REQUIRE(documentedRatios.size() == refinedRateRatios.size());
  const double documentedRatioScale =
      std::pow(10.0, static_cast<double>(documentedRatioDecimalPlaces));
  const double documentedRatioRounding = 0.5 / documentedRatioScale;
  for (size_t i = 0; i < refinedRateRatios.size(); ++i) {
    // Lock the four decimal places printed by the generated guide table. The
    // half-unit rounding interval is wide enough for differences below the
    // displayed resolution without requiring bit-reproducible linear algebra.
    REQUIRE_THAT(refinedRateRatios[i],
                 WithinAbs(documentedRatios[i],
                           documentedRatioRounding));
    REQUIRE(std::isfinite(refinedPartitionDifferences[i]));
  }
  REQUIRE_FALSE(forceCalls.empty());
  REQUIRE(forceCalls.size() == callRows.size());
  REQUIRE(hessianProducts.size() == callRows.size());
  REQUIRE(std::all_of(forceCalls.begin(), forceCalls.end(),
                      [](size_t calls) { return calls > 0; }));
  REQUIRE(std::all_of(hessianProducts.begin(), hessianProducts.end(),
                      [](long products) { return products > 0; }));

}

TEST_CASE("CuH2 line-integral job runs end to end",
          "[LineIntegralTunneling][CuH2]") {
  Parameters parameters;
  ParametersLoadAccess::potential_options(parameters).potential = PotType::CUH2;
  auto potential = helpers::sharePotential(
      helpers::makePotential(PotType::CUH2, parameters));
  Matter reactant(potential, parameters);
  REQUIRE(io::io_ok(reactant.con2matter("min2_0Cu.con")));
  TunnelingJobContext context(
      parameters, potential, reactant,
      {.requireWholeAtomConstraints = true, .forwardNonPeriodicCell = true});
  const MassWeighted &coordinates = context.coordinates();
  const std::vector<VectorXd> mep = loadMep(reactant, coordinates);
  std::vector<double> energies;
  std::vector<VectorXd> gradients;
  context.evaluate(mep, energies, gradients);
  const size_t saddle = static_cast<size_t>(
      std::max_element(energies.begin(), energies.end()) - energies.begin());
  REQUIRE(saddle > 0);
  REQUIRE(saddle + 1 < mep.size());

  const std::vector<std::array<double, 2>> expected = documentedJobRates();
  const std::filesystem::path jobDirectory =
      std::filesystem::temp_directory_path() / "eon_line_integral_job_test";
  std::filesystem::remove_all(jobDirectory);
  std::filesystem::create_directories(jobDirectory);
  CurrentPathGuard currentPath;
  std::filesystem::current_path(jobDirectory);
  writeJobInputs(reactant, coordinates, mep, saddle);

  auto jobParameters = std::make_unique<Parameters>();
  ParametersLoadAccess::potential_options(*jobParameters).potential =
      PotType::CUH2;
  auto &jobOptions =
      ParametersLoadAccess::line_integral_tunneling_options(*jobParameters);
  jobOptions.images = 10;
  jobOptions.replicas = 150;
  jobOptions.energy_offset = 0.05;
  jobOptions.lowest_temperature = 0.5 * (expected[0][0] + expected[1][0]);
  jobOptions.energy_step = 0.01;
  jobOptions.log_rate_tolerance = 0.5 * std::log(10.0);
  jobOptions.spring = 20.0;
  jobOptions.contour = 1000.0;
  jobOptions.force_tolerance = 0.01;
  const TunnelingJobContext::LowestMode unstable =
      context.lowestModeAt(mep[saddle], mep[saddle + 1] - mep[saddle - 1]);
  REQUIRE(unstable.curvature < 0.0);
  const double weakestCurvature = -unstable.curvature - unstable.residual;
  REQUIRE(weakestCurvature > 0.0);
  const double crossoverPeriod =
      2.0 * std::numbers::pi / std::sqrt(weakestCurvature);
  jobOptions.period_tolerance =
      std::sqrt(std::numeric_limits<double>::epsilon()) * crossoverPeriod +
      jobOptions.force_tolerance / jobOptions.contour;
  const auto firstAtOrBelowTarget =
      std::find_if(expected.begin(), expected.end(), [&](const auto &row) {
        return row[0] <= jobOptions.lowest_temperature;
      });
  REQUIRE(firstAtOrBelowTarget != expected.end());
  const size_t expectedJobLevels =
      static_cast<size_t>(std::distance(expected.begin(), firstAtOrBelowTarget)) +
      1;
  REQUIRE(expectedJobLevels >= 2);
  constexpr double jobTemperatureTolerance = 0.05;
  constexpr double jobLogRateTolerance = 0.05;
  Runtime runtime;
  LineIntegralTunnelingJob job(std::move(jobParameters), runtime);
  const std::vector<std::string> outputs = job.run();
  REQUIRE(std::find(outputs.begin(), outputs.end(),
                    "line_integral_tunneling.dat") != outputs.end());
  REQUIRE(std::filesystem::exists("line_integral_path_000.con"));
  REQUIRE(std::filesystem::exists("line_integral_path_001.con"));
  REQUIRE(std::filesystem::exists("results.dat"));
  const std::vector<std::vector<double>> jobRows =
      readNumericTable("line_integral_tunneling.dat");
  REQUIRE(jobRows.size() == expectedJobLevels);
  for (size_t i = 0; i < jobRows.size(); ++i) {
    REQUIRE(jobRows[i].size() == 16);
    REQUIRE_THAT(jobRows[i][3],
                 WithinAbs(expected[i][0], jobTemperatureTolerance));
    REQUIRE_THAT(jobRows[i][6],
                 WithinAbs(expected[i][1], jobLogRateTolerance));
  }
  const std::map<std::string, std::string> results =
      readResults("results.dat");
  REQUIRE(results.at("termination_reason") == "0");
  REQUIRE(results.at("termination_reason_text") == "GOOD");
  REQUIRE(static_cast<size_t>(
              std::stol(results.at("line_integral_levels"))) ==
          expectedJobLevels);
  REQUIRE_THAT(std::stod(results.at("line_integral_log_rate_internal")),
               WithinAbs(expected[expectedJobLevels - 1][1],
                         jobLogRateTolerance));
}

TEST_CASE("Line-integral tunnelling declares its finite-difference fallback",
          "[LineIntegralTunneling][CuH2]") {
  Parameters parameters;
  ParametersLoadAccess::potential_options(parameters).potential =
      PotType::LJCLUSTER;
  auto potential = helpers::sharePotential(
      helpers::makePotential(PotType::LJCLUSTER, parameters));
  REQUIRE_FALSE(potential->supportsMassWeightedHessianVectorProduct());
  Matter reactant(potential, parameters);
  REQUIRE(io::io_ok(reactant.con2matter("min2_0Cu.con")));
  TunnelingJobContext context(
      parameters, potential, reactant,
      {.requireWholeAtomConstraints = true, .forwardNonPeriodicCell = true});
  REQUIRE(context.usesFiniteDifferenceHessianFallback());
  const MassWeighted &coordinates = context.coordinates();
  const std::vector<VectorXd> mep = loadMep(reactant, coordinates);
  const size_t middle = mep.size() / 2;
  const HessianProduct product = context.hessianProductAt(mep[middle]);
  const VectorXd direction = mep[middle + 1] - mep[middle - 1];
  REQUIRE(product(direction).allFinite());
  const RingHessianProducts ring = context.ringHessianProducts(
      {mep[middle - 1], mep[middle], mep[middle + 1]});
  std::vector<VectorXd> applied;
  ring({direction, direction, direction}, applied);
  REQUIRE(applied.size() == 3);
  for (const VectorXd &value : applied) {
    REQUIRE(value.allFinite());
  }
}

TEST_CASE("tunnelling context rejects non-finite potential energies",
          "[LineIntegralTunneling]") {
  Parameters parameters;
  auto potential = std::make_shared<NonFiniteEnergyPotential>(parameters);
  Matter reactant(potential, parameters);
  REQUIRE(io::io_ok(reactant.con2matter("min2_0Cu.con")));
  TunnelingJobContext context(parameters, potential, reactant);
  std::vector<double> energies;
  std::vector<VectorXd> gradients;
  REQUIRE_THROWS_AS(
      context.evaluate({VectorXd::Zero(context.coordinates().dimension())},
                       energies, gradients),
      std::runtime_error);
}

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

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <numbers>
#include <vector>

using namespace eonc;
using namespace eonc::tunneling;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double kTransverse = 1.5;
constexpr double kModifiedMullerBrownUnitsPerEv = 600.0;

void quarticPotential(const std::vector<VectorXd> &points,
                      std::vector<double> &energies,
                      std::vector<VectorXd> &gradients) {
  energies.resize(points.size());
  gradients.resize(points.size());
  for (size_t i = 0; i < points.size(); ++i) {
    const double x = points[i](0);
    const double y = points[i](1);
    const double x2 = x * x;
    energies[i] = (1.0 - x2) * (1.0 - x2) + 0.5 * kTransverse * y * y;
    gradients[i].resize(2);
    gradients[i] << 4.0 * x * (x2 - 1.0), kTransverse * y;
  }
}

MatrixXd quarticHessian(const VectorXd &point) {
  MatrixXd hessian = MatrixXd::Zero(2, 2);
  hessian(0, 0) = 12.0 * point(0) * point(0) - 4.0;
  hessian(1, 1) = kTransverse;
  return hessian;
}

struct Gaussian {
  double amplitude;
  double xx;
  double xy;
  double yy;
  double x;
  double y;
};

// The first four rows are the canonical Müller-Brown surface. The JCP paper
// says that Fig. 2 adds three Gaussians to shift and deepen the right minimum,
// but does not tabulate them. The final three rows are a reproducible
// reconstruction from that figure: one Gaussian removes the canonical basin
// at (0.623, 0.028), and two centred at the digitized shifted basin reproduce
// its displacement and depth. They are deliberately identified as a figure
// reconstruction rather than unpublished author parameters.
constexpr Gaussian kModifiedMullerBrown[] = {
    {-200.0, -1.0, 0.0, -10.0, 1.0, 0.0},
    {-100.0, -1.0, 0.0, -10.0, 0.0, 0.5},
    {-170.0, -6.5, 11.0, -6.5, -0.5, 1.5},
    {15.0, 0.7, 0.6, 0.7, -1.0, 1.0},
    {100.0, -10.0, 0.0, -10.0, 0.623, 0.028},
    {-100.0, -10.0, 0.0, -10.0, 1.0, 0.0},
    {-50.0, -5.0, 0.0, -5.0, 1.2, 0.0},
};

void modifiedMullerBrown(const std::vector<VectorXd> &points,
                         std::vector<double> &energies,
                         std::vector<VectorXd> &gradients) {
  energies.assign(points.size(), 0.0);
  gradients.assign(points.size(), VectorXd::Zero(2));
  for (size_t j = 0; j < points.size(); ++j) {
    for (const Gaussian &term : kModifiedMullerBrown) {
      const double dx = points[j](0) - term.x;
      const double dy = points[j](1) - term.y;
      const double exponent =
          term.xx * dx * dx + term.xy * dx * dy + term.yy * dy * dy;
      const double value =
          term.amplitude * std::exp(exponent) / kModifiedMullerBrownUnitsPerEv;
      energies[j] += value;
      gradients[j](0) += value * (2.0 * term.xx * dx + term.xy * dy);
      gradients[j](1) += value * (term.xy * dx + 2.0 * term.yy * dy);
    }
  }
}

std::pair<std::vector<VectorXd>, std::vector<double>> modifiedMullerBrownMep() {
  std::vector<VectorXd> mep;
  std::ifstream input("mep.dat");
  REQUIRE(input.good());
  double x = 0.0;
  double y = 0.0;
  while (input >> x >> y) {
    mep.push_back((VectorXd(2) << x, y).finished());
  }
  REQUIRE(mep.size() == 81);
  std::vector<double> energies;
  std::vector<VectorXd> gradients;
  modifiedMullerBrown(mep, energies, gradients);
  return {mep, energies};
}

std::pair<std::vector<VectorXd>, std::vector<double>> quarticMep(long images) {
  std::vector<VectorXd> path;
  std::vector<double> energies;
  for (long i = 0; i < images; ++i) {
    VectorXd point(2);
    point << -1.0 +
                 2.0 * static_cast<double>(i) / static_cast<double>(images - 1),
        0.0;
    path.push_back(point);
    energies.push_back(std::pow(1.0 - point(0) * point(0), 2));
  }
  return {path, energies};
}

LineIntegralPath exactQuarticPath(double energy, long images) {
  const double turningPoint = std::sqrt(1.0 - std::sqrt(energy));
  LineIntegralPath path;
  path.energy = energy;
  path.converged = true;
  for (long i = 0; i < images; ++i) {
    const double fraction =
        static_cast<double>(i) / static_cast<double>(images - 1);
    path.images.push_back(
        (VectorXd(2) << (2.0 * fraction - 1.0) * turningPoint, 0.0).finished());
  }
  quarticPotential(path.images, path.energies, path.gradients);
  const long intervals = 2 * (images - 1);
  const double step = (std::numbers::pi / 2.0) / static_cast<double>(intervals);
  auto transformedAction = [&](double theta) {
    const double sine = std::sin(theta);
    const double cosine = std::cos(theta);
    const double a2 = turningPoint * turningPoint;
    const double excess =
        a2 * cosine * cosine * (1.0 + std::sqrt(energy) - a2 * sine * sine);
    return std::sqrt(2.0 * excess) * turningPoint * cosine;
  };
  double halfAction =
      transformedAction(0.0) + transformedAction(std::numbers::pi / 2.0);
  for (long i = 1; i < intervals; ++i) {
    halfAction += (i % 2 == 0 ? 2.0 : 4.0) *
                  transformedAction(step * static_cast<double>(i));
  }
  path.action = 2.0 * step * halfAction / 3.0;
  return path;
}

double logPartitionRatio(const RateInstanton &rate) {
  const double bnh = rate.betaN * kHbar;
  const double zeroModeFactor =
      0.5 *
      std::log(rate.bN / (2.0 * std::numbers::pi * rate.betaN * kHbar * kHbar));
  return rate.logRate + std::log(bnh) - zeroModeFactor +
         rate.betaN * rate.ringPotential;
}

} // namespace

TEST_CASE("LI-NEB corner-cuts the reconstructed modified Muller-Brown surface",
          "[LineIntegralTunneling][MullerBrown]") {
  // Figure 3 was digitized at 0.1 displacement intervals. Integrating
  // sqrt(2 Delta V) over its red curve for unit mass gives W=0.81; reading
  // the 0.02 eV and 0.5 displacement tick spacing gives uncertainties of
  // 0.01 eV and 0.05, which propagate to 0.04 in W. The plotted endpoint
  // displacement is 2.05 with a 0.05 per-image reading uncertainty; independent
  // image errors are propagated in quadrature along the path.
  constexpr double figureAction = 0.81;
  constexpr double figureActionUncertainty = 0.04;
  constexpr double figureDisplacement = 2.05;
  constexpr double figureCoordinateUncertainty = 0.05;
  LineIntegralOptions options;
  const double figureDisplacementUncertainty =
      std::sqrt(static_cast<double>(options.images - 1)) *
      figureCoordinateUncertainty;
  const auto [mep, mepEnergies] = modifiedMullerBrownMep();
  const double energy = -110.0 / kModifiedMullerBrownUnitsPerEv;
  const std::vector<VectorXd> seed =
      lineIntegralSeed(mep, mepEnergies, energy, options.images);
  const LineIntegralPath path =
      optimizeLineIntegralPath(seed, energy, modifiedMullerBrown, options);
  CAPTURE(path.action, path.iterations, path.interiorGradient,
          path.endpointGradient);
  double pathLength = 0.0;
  for (size_t i = 1; i < path.images.size(); ++i) {
    pathLength += (path.images[i] - path.images[i - 1]).norm();
  }
  std::vector<double> relativeEnergies = path.energies;
  for (double &value : relativeEnergies) {
    value -= energy;
  }
  CAPTURE(pathLength, path.images, relativeEnergies);
  REQUIRE(path.converged);
  REQUIRE_THAT(path.action, WithinAbs(figureAction, figureActionUncertainty));
  REQUIRE_THAT(pathLength,
               WithinAbs(figureDisplacement, figureDisplacementUncertainty));

  const VectorXd intermediate =
      (VectorXd(2) << -0.05860977, 0.46837729).finished();
  auto closest = [&](const std::vector<VectorXd> &points) {
    double distance = (points.front() - intermediate).norm();
    for (const VectorXd &point : points) {
      distance = std::min(distance, (point - intermediate).norm());
    }
    return distance;
  };
  CAPTURE(closest(seed), closest(path.images));
  REQUIRE(closest(path.images) > closest(seed));
}

TEST_CASE("LI-NEB continuation follows exact quartic turning contours",
          "[LineIntegralTunneling]") {
  const auto [mep, mepEnergies] = quarticMep(101);
  LineIntegralOptions options;
  options.images = 10;
  options.spring = 20.0;
  options.contour = 1000.0;
  options.forceTolerance = 0.01;
  options.maxIterations = 5000;
  options.timeStep = 0.01;

  const double firstEnergy = 0.36;
  const auto firstSeed =
      lineIntegralSeed(mep, mepEnergies, firstEnergy, options.images);
  const LineIntegralPath first = optimizeLineIntegralPath(
      firstSeed, firstEnergy, quarticPotential, options);
  REQUIRE(first.converged);
  const double firstContourEnergyError =
      first.endpointGradient / options.contour;
  REQUIRE_THAT(first.energies.front(),
               WithinAbs(firstEnergy, firstContourEnergyError));
  REQUIRE_THAT(first.energies.back(),
               WithinAbs(firstEnergy, firstContourEnergyError));

  const double nextEnergy = 0.25;
  const auto nextSeed = lineIntegralContinuation(mep, mepEnergies, first,
                                                 nextEnergy, options.images);
  const LineIntegralPath next =
      optimizeLineIntegralPath(nextSeed, nextEnergy, quarticPotential, options);
  REQUIRE(next.converged);
  REQUIRE(next.action > first.action);
  const double turningPoint = std::sqrt(1.0 - std::sqrt(nextEnergy));
  const double endpointEnergyError = next.endpointGradient / options.contour;
  const double outerTurningPoint =
      std::sqrt(1.0 - std::sqrt(nextEnergy - endpointEnergyError));
  const double innerTurningPoint =
      std::sqrt(1.0 - std::sqrt(nextEnergy + endpointEnergyError));
  const double endpointPositionError = std::max(
      outerTurningPoint - turningPoint, turningPoint - innerTurningPoint);
  REQUIRE_THAT(next.images.front()(0),
               WithinAbs(-turningPoint, endpointPositionError));
  REQUIRE_THAT(next.images.back()(0),
               WithinAbs(turningPoint, endpointPositionError));
}

TEST_CASE("LI-NEB accepts an eOn optimizer", "[LineIntegralTunneling]") {
  const auto [mep, mepEnergies] = quarticMep(101);
  LineIntegralOptions options;
  const double energy = 0.36;
  Parameters parameters;
  const LineIntegralPath path = optimizeLineIntegralPath(
      lineIntegralSeed(mep, mepEnergies, energy, options.images), energy,
      quarticPotential, options, parameters, OptType::FIRE);
  REQUIRE(path.converged);
  REQUIRE(path.iterations < options.maxIterations);
}

TEST_CASE("LI-NEB rejects a multi-barrier MEP", "[LineIntegralTunneling]") {
  std::vector<VectorXd> mep;
  for (long i = 0; i < 5; ++i) {
    mep.push_back((VectorXd(2) << static_cast<double>(i), 0.0).finished());
  }
  const std::vector<double> energies{0.0, 2.0, 0.0, 1.5, 0.0};
  REQUIRE_THROWS_AS(lineIntegralSeed(mep, energies, 1.0, 10),
                    std::invalid_argument);
}

TEST_CASE("energy continuation limits the predicted log-rate change",
          "[LineIntegralTunneling]") {
  constexpr double initialEnergy = 0.5;
  constexpr double maximumStep = 0.02;
  constexpr double logRateTolerance = 0.5;
  constexpr double lowestTemperature = 80.0;
  std::vector<LineIntegralContinuationPoint> points{
      {initialEnergy, 100.0, 0.0}};
  REQUIRE_THAT(nextLineIntegralEnergy(points, maximumStep, logRateTolerance,
                                      lowestTemperature),
               WithinAbs(initialEnergy - maximumStep,
                         std::numeric_limits<double>::epsilon()));
  const double secondEnergy = initialEnergy - maximumStep / 2.0;
  const double secondTemperature =
      (points.front().temperature + lowestTemperature) / 2.0;
  const double secondLogRate = 2.0 * logRateTolerance;
  points.push_back({secondEnergy, secondTemperature, secondLogRate});
  const double measuredEnergyStep = points[0].energy - points[1].energy;
  const double measuredLogRateStep =
      std::abs(points[1].logRate - points[0].logRate);
  const double controlledStep =
      logRateTolerance * measuredEnergyStep / measuredLogRateStep;
  REQUIRE_THAT(nextLineIntegralEnergy(points, maximumStep, logRateTolerance,
                                      lowestTemperature),
               WithinAbs(points.back().energy - controlledStep,
                         std::numeric_limits<double>::epsilon()));
}

TEST_CASE("HQTST agrees with the ring instanton on a separable barrier",
          "[LineIntegralTunneling]") {
  const double energy = 0.25;
  const double rootEnergy = std::sqrt(energy);
  const double ellipticModulus =
      std::sqrt((1.0 - rootEnergy) / (1.0 + rootEnergy));
  const double exactPeriod = 4.0 * std::comp_ellint_1(ellipticModulus) /
                             std::sqrt(2.0 * (1.0 + rootEnergy));
  const double periodTolerance =
      std::sqrt(std::numeric_limits<double>::epsilon()) * exactPeriod;
  constexpr long pathImages = 80;
  constexpr long replicas = 128;
  const LineIntegralPath path = exactQuarticPath(energy, pathImages);
  const LineIntegralOrbit orbit =
      lineIntegralOrbit(path, replicas, periodTolerance, std::numbers::pi);
  CAPTURE(orbit.period, orbit.periodError, exactPeriod);
  REQUIRE_THAT(orbit.period, WithinAbs(exactPeriod, orbit.periodError));
  const VectorXd reactant = (VectorXd(2) << -1.0, 0.0).finished();
  const VectorXd saddle = VectorXd::Zero(2);
  const MatrixXd hReactant = quarticHessian(reactant);
  const MatrixXd hSaddle = quarticHessian(saddle);
  const HqtstOptions hqtst{.estimatedLogRateErrorTarget = 0.01, .seed = 3};
  auto rateOn = [&](const LineIntegralOrbit &candidate) {
    const std::vector<VectorXd> replicaPoints = candidate.replicas;
    long products = 0;
    const HqtstProducts exact{
        .ring =
            [replicaPoints, &products](const std::vector<VectorXd> &v,
                                       std::vector<VectorXd> &out) {
              out.resize(v.size());
              for (size_t j = 0; j < v.size(); ++j) {
                out[j] = quarticHessian(replicaPoints[j]) * v[j];
                ++products;
              }
            },
        .reactant =
            [&](const VectorXd &v) {
              ++products;
              return VectorXd(hReactant * v);
            }};
    const RateInstanton rate =
        hqtstRate(candidate, exact, 0.0, MatrixXd(2, 0), hqtst);
    REQUIRE(rate.hessianProducts == products);
    REQUIRE(rate.logRateError <= hqtst.estimatedLogRateErrorTarget);
    return rate;
  };
  const RateInstanton direct = rateOn(orbit);
  REQUIRE(direct.negativeModes == 1);

  RateInstantonOptions instantonOptions;
  instantonOptions.beads = replicas;
  instantonOptions.maxIterations = 5000;
  // The ring discretization is second order in 1/P, so the optimizer target
  // is scaled from that same resolution rather than a fixed literal.
  const double instantonForceTolerance =
      4.0 / (static_cast<double>(replicas) * static_cast<double>(replicas));
  REQUIRE(instantonForceTolerance > 0.0);
  instantonOptions.forceTolerance = instantonForceTolerance;
  instantonOptions.halfRing = true;
  const double beta = orbit.period / kHbar;
  RateInstanton instanton =
      optimizeRateInstanton(saddle, hSaddle, beta, orbit.replicas,
                            quarticPotential, instantonOptions);
  REQUIRE(instanton.converged);
  instantonRate(
      instanton,
      [](long, const VectorXd &point) { return quarticHessian(point); },
      hReactant, 0.0, hSaddle, 1.0);
  REQUIRE(instanton.negativeModes == 1);
  constexpr long refinedReplicas = 256;
  const LineIntegralOrbit refinedOrbit = lineIntegralOrbit(
      path, refinedReplicas, periodTolerance, std::numbers::pi);
  const RateInstanton refinedDirect = rateOn(refinedOrbit);
  REQUIRE(refinedDirect.negativeModes == 1);
  instantonOptions.beads = refinedReplicas;
  RateInstanton refinedInstanton =
      optimizeRateInstanton(saddle, hSaddle, beta, refinedOrbit.replicas,
                            quarticPotential, instantonOptions);
  REQUIRE(refinedInstanton.converged);
  instantonRate(
      refinedInstanton,
      [](long, const VectorXd &point) { return quarticHessian(point); },
      hReactant, 0.0, hSaddle, 1.0);
  REQUIRE(refinedInstanton.negativeModes == 1);

  // Each HQTST value carries its determinant error bar.
  const double determinantEstimatedError =
      direct.logRateError + refinedDirect.logRateError;
  const double rateDiscretizationChange =
      std::abs(refinedInstanton.logRate - instanton.logRate) +
      std::abs(refinedDirect.logRate - direct.logRate);
  const double coarseRateGap = std::abs(direct.logRate - instanton.logRate);
  const double refinedRateGap =
      std::abs(refinedDirect.logRate - refinedInstanton.logRate);
  const double partitionDiscretizationChange =
      std::abs(logPartitionRatio(refinedInstanton) -
               logPartitionRatio(instanton)) +
      std::abs(logPartitionRatio(refinedDirect) - logPartitionRatio(direct));
  const double coarsePartitionGap =
      std::abs(logPartitionRatio(direct) - logPartitionRatio(instanton));
  const double refinedPartitionGap = std::abs(
      logPartitionRatio(refinedDirect) - logPartitionRatio(refinedInstanton));
  CAPTURE(orbit.action, orbit.period, orbit.temperature, direct.ringPotential,
          instanton.ringPotential, direct.bN, instanton.bN, direct.logRate,
          instanton.logRate, refinedDirect.logRate, refinedInstanton.logRate,
          determinantEstimatedError, rateDiscretizationChange,
          coarsePartitionGap, refinedPartitionGap,
          partitionDiscretizationChange);
  REQUIRE(refinedRateGap < coarseRateGap + determinantEstimatedError);
  REQUIRE(refinedRateGap <
          rateDiscretizationChange + determinantEstimatedError);
  REQUIRE(refinedPartitionGap < coarsePartitionGap + determinantEstimatedError);
  REQUIRE(refinedPartitionGap <
          partitionDiscretizationChange + determinantEstimatedError);
}

TEST_CASE("HQTST needs every product it is given", "[LineIntegralTunneling]") {
  const LineIntegralPath path = exactQuarticPath(0.25, 40);
  const LineIntegralOrbit orbit = lineIntegralOrbit(
      path, 32, std::sqrt(std::numeric_limits<double>::epsilon()),
      std::numbers::pi);
  const HqtstProducts missing{
      .ring = {}, .reactant = [](const VectorXd &v) { return VectorXd(v); }};
  REQUIRE_THROWS_AS(hqtstRate(orbit, missing, 0.0, MatrixXd(2, 0), {}),
                    std::invalid_argument);
}

TEST_CASE("line-integral inputs fail closed on non-finite values",
          "[LineIntegralTunneling]") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  LineIntegralPath path = exactQuarticPath(0.25, 40);
  const double periodTolerance =
      std::sqrt(std::numeric_limits<double>::epsilon());

  path.action = nan;
  REQUIRE_THROWS_AS(
      lineIntegralOrbit(path, 32, periodTolerance, std::numbers::pi),
      std::invalid_argument);
  path = exactQuarticPath(0.25, 40);
  path.action = -1.0;
  REQUIRE_THROWS_AS(
      lineIntegralOrbit(path, 32, periodTolerance, std::numbers::pi),
      std::invalid_argument);
  path = exactQuarticPath(0.25, 40);
  path.energy = infinity;
  REQUIRE_THROWS_AS(
      lineIntegralOrbit(path, 32, periodTolerance, std::numbers::pi),
      std::invalid_argument);
  path = exactQuarticPath(0.25, 40);
  REQUIRE_THROWS_AS(lineIntegralOrbit(path, 32, nan, std::numbers::pi),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(lineIntegralOrbit(path, 32, periodTolerance, nan),
                    std::invalid_argument);

  const auto [mep, mepEnergies] = quarticMep(101);
  LineIntegralOptions options;
  options.spring = infinity;
  REQUIRE_THROWS_AS(
      optimizeLineIntegralPath(
          lineIntegralSeed(mep, mepEnergies, 0.36, options.images), 0.36,
          quarticPotential, options),
      std::invalid_argument);
  options = {};
  REQUIRE_THROWS_AS(
      optimizeLineIntegralPath(
          lineIntegralSeed(mep, mepEnergies, 0.36, options.images), nan,
          quarticPotential, options),
      std::invalid_argument);

  options = {};
  const auto nonFiniteGradient =
      [nan](const std::vector<VectorXd> &points,
            std::vector<double> &energies,
            std::vector<VectorXd> &gradients) {
        quarticPotential(points, energies, gradients);
        gradients[points.size() / 2](0) = nan;
      };
  REQUIRE_THROWS_AS(
      optimizeLineIntegralPath(
          lineIntegralSeed(mep, mepEnergies, 0.36, options.images), 0.36,
          nonFiniteGradient, options),
      std::runtime_error);
}

TEST_CASE("HQTST rejects every non-finite scalar input",
          "[LineIntegralTunneling]") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  const LineIntegralPath path = exactQuarticPath(0.25, 40);
  const LineIntegralOrbit valid = lineIntegralOrbit(
      path, 32, std::sqrt(std::numeric_limits<double>::epsilon()),
      std::numbers::pi);
  const HqtstProducts products{
      .ring = [](const std::vector<VectorXd> &v,
                 std::vector<VectorXd> &out) { out = v; },
      .reactant = [](const VectorXd &v) { return VectorXd(v); }};
  const MatrixXd rigid(2, 0);

  for (const auto member : {&LineIntegralOrbit::energy,
                            &LineIntegralOrbit::action,
                            &LineIntegralOrbit::period}) {
    LineIntegralOrbit orbit = valid;
    orbit.*member = nan;
    REQUIRE_THROWS_AS(hqtstRate(orbit, products, 0.0, rigid, {}),
                      std::invalid_argument);
  }
  LineIntegralOrbit orbit = valid;
  orbit.period = infinity;
  REQUIRE_THROWS_AS(hqtstRate(orbit, products, 0.0, rigid, {}),
                    std::invalid_argument);
  orbit = valid;
  orbit.energies.front() = nan;
  REQUIRE_THROWS_AS(hqtstRate(orbit, products, 0.0, rigid, {}),
                    std::invalid_argument);
  orbit = valid;
  orbit.replicas.front().setConstant(
      std::sqrt(std::numeric_limits<double>::max()));
  REQUIRE_THROWS_AS(hqtstRate(orbit, products, 0.0, rigid, {}),
                    std::runtime_error);
  REQUIRE_THROWS_AS(hqtstRate(valid, products, nan, rigid, {}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(
      hqtstRate(valid, products, std::numeric_limits<double>::max(), rigid, {}),
      std::runtime_error);
  HqtstOptions options;
  options.estimatedLogRateErrorTarget = infinity;
  REQUIRE_THROWS_AS(hqtstRate(valid, products, 0.0, rigid, options),
                    std::invalid_argument);
}

TEST_CASE("line-integral option validation rejects infinity",
          "[LineIntegralTunneling]") {
  const double infinity = std::numeric_limits<double>::infinity();
  line_integral_tunneling_options_t options;
  const std::array<double line_integral_tunneling_options_t::*, 10> fields{
      &line_integral_tunneling_options_t::energy_offset,
      &line_integral_tunneling_options_t::lowest_temperature,
      &line_integral_tunneling_options_t::energy_step,
      &line_integral_tunneling_options_t::log_rate_tolerance,
      &line_integral_tunneling_options_t::spring,
      &line_integral_tunneling_options_t::contour,
      &line_integral_tunneling_options_t::force_tolerance,
      &line_integral_tunneling_options_t::li_neb_time_step,
      &line_integral_tunneling_options_t::period_tolerance,
      &line_integral_tunneling_options_t::determinant_estimated_error};
  for (const auto field : fields) {
    auto invalid = options;
    invalid.*field = infinity;
    REQUIRE_THROWS_AS(validateOptions(invalid), std::invalid_argument);
  }
}

TEST_CASE("period uncertainty rejects quartic paths at crossover",
          "[LineIntegralTunneling]") {
  constexpr long images = 10;
  constexpr long replicas = 128;
  const double crossoverPeriod = std::numbers::pi;
  const double periodTolerance =
      std::sqrt(std::numeric_limits<double>::epsilon()) * crossoverPeriod;
  // Regression energies from the reported crossover failure: both lie inside
  // the independently estimated path-discretization uncertainty.
  const std::array<double, 2> barrierDeficits{0.003, 0.001};
  for (double deficit : barrierDeficits) {
    const LineIntegralPath path = exactQuarticPath(1.0 - deficit, images);
    REQUIRE_THROWS_AS(
        lineIntegralOrbit(path, replicas, periodTolerance, crossoverPeriod),
        CrossoverTemperatureError);
  }
}

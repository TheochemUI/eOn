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
#pragma once

#include "BaseStructures.h"
#include "Eigen.h"
#include "Parameters.h"
#include "RingDeterminant.h"
#include "Tunneling.h"

#include <cstdint>
#include <functional>
#include <stdexcept>
#include <vector>

namespace eonc::tunneling {

/// Validate the shared INI/JSON line-integral job options. Errors name the
/// rejected key so configuration failures are actionable.
void validateOptions(const line_integral_tunneling_options_t &options);

struct LineIntegralOptions {
  long images{10};
  double spring{20.0};
  double contour{1000.0};
  double forceTolerance{0.01};
  long maxIterations{5000};
  double timeStep{0.01};
};

using PathAlignment = std::function<void(std::vector<VectorXd> &)>;

struct LineIntegralPath {
  std::vector<VectorXd> images;
  std::vector<double> energies;
  std::vector<VectorXd> gradients;
  double energy{0.0};
  double action{0.0};
  double interiorGradient{0.0};
  double endpointGradient{0.0};
  long iterations{0};
  bool converged{false};
};

/// Points on the MEP between its two crossings of energy, resampled to an
/// evenly spaced path. The highest-energy MEP image separates the crossings.
std::vector<VectorXd> lineIntegralSeed(const std::vector<VectorXd> &mep,
                                       const std::vector<double> &energies,
                                       double energy, long images);

/// Lower-energy continuation seed: attach the newly exposed MEP pieces to the
/// previous optimal path, then resample to the requested image count.
std::vector<VectorXd> lineIntegralContinuation(
    const std::vector<VectorXd> &mep, const std::vector<double> &energies,
    const LineIntegralPath &previous, double energy, long images);

/// Resample an existing path to a different evenly spaced image count.
std::vector<VectorXd>
resampleLineIntegralPath(const std::vector<VectorXd> &path, long images);

struct LineIntegralContinuationPoint {
  double energy;
  double temperature;
  double logRate;
};

/// Select the next lower energy from the local changes in temperature and
/// log rate. The configured energy step is an upper bound.
double
nextLineIntegralEnergy(const std::vector<LineIntegralContinuationPoint> &points,
                       double maximumEnergyStep, double logRateTolerance,
                       double lowestTemperature);

/// LI-NEB minimization with velocity-projection velocity Verlet.
LineIntegralPath optimizeLineIntegralPath(std::vector<VectorXd> guess,
                                          double energy,
                                          const BatchPotential &potential,
                                          const LineIntegralOptions &options,
                                          const PathAlignment &align = {});

/// LI-NEB minimization through one of eOn's ordinary optimizers.
LineIntegralPath optimizeLineIntegralPath(
    std::vector<VectorXd> guess, double energy, const BatchPotential &potential,
    const LineIntegralOptions &options, const Parameters &parameters,
    OptType optimizer, const PathAlignment &align = {});

struct LineIntegralOrbit {
  std::vector<VectorXd> replicas;
  std::vector<double> energies;
  std::vector<double> coordinate;
  double period{0.0};
  double periodError{0.0};
  double temperature{0.0};
  double energy{0.0};
  double action{0.0};
};

class CrossoverTemperatureError : public std::domain_error {
public:
  using std::domain_error::domain_error;
};

/// Integrate imaginary time along the cubic path and mirror that half orbit
/// onto P equally spaced replicas.
LineIntegralOrbit lineIntegralOrbit(const LineIntegralPath &path, long replicas,
                                    double periodTolerance,
                                    double crossoverPeriod);

/// The exact Hessian-vector products HQTST is evaluated from.
struct HqtstProducts {
  /// H(q_j) v_j at every orbit replica.
  RingHessianProducts ring;
  /// H v at the reactant minimum.
  HessianProduct reactant;
};

struct HqtstOptions {
  /// Stopping threshold for the estimated ln k error from the two stochastic
  /// log-determinants; ln k carries half of each diagnostic estimate.
  double estimatedLogRateErrorTarget{0.1};
  std::uint64_t seed{0};
};

/// Ásgeirsson-Arnaldsson-Jónsson harmonic quantum TST with the partition
/// functions of InstantonJob, from Hessian-vector products alone: the ring
/// determinant is ringLogDetPrime on the replicas, the reactant's is
/// wellLogDet over the ring's free-particle frequencies. `reactantRigid`
/// holds the reactant's orthonormal rigid directions as columns; on the
/// ring they are the beads' own rigid motions when `rigidBodies` has
/// masses, and otherwise those columns copied to every bead.
RateInstanton hqtstRate(const LineIntegralOrbit &orbit,
                        const HqtstProducts &products, double reactantEnergy,
                        const MatrixXd &reactantRigid,
                        const HqtstOptions &options,
                        const RingRigidBodies &rigidBodies = {});

} // namespace eonc::tunneling

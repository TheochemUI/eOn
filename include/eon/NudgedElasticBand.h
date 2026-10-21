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

#include <filesystem>

#include "Eigen.h"
#include "EonLogger.h"

#include "EigenmodeStrategy.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "NEBProjection.h"
#include "NEBTangent.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"

#include <optional>
#include <readcon-core.hpp>
#include <vector>

namespace eonc::neb {
class OCINEBController;
}

namespace eonc {

// NEB method for determining a minimum energy path between two matter objects
class NudgedElasticBand {
  friend class eonc::neb::OCINEBController;

public:
  enum class NEBStatus {
    GOOD = 0,
    INIT = 1,
    BAD_MAX_ITERATIONS = 2,
    RUNNING,
    MAX_UNCERTAINTY
  };
  NudgedElasticBand(std::shared_ptr<Matter> initialPassed,
                    std::shared_ptr<Matter> finalPassed,
                    const Parameters &parametersPassed,
                    std::shared_ptr<Potential> potPassed);
  NudgedElasticBand(std::vector<Matter> initPath,
                    const Parameters &parametersPassed,
                    std::shared_ptr<Potential> potPassed);
  ~NudgedElasticBand() = default;

  NudgedElasticBand::NEBStatus compute(void);
  NudgedElasticBand::NEBStatus getStatus() { return this->status; };
  /// Iteration count from the last compute(). Zero before compute() runs.
  [[nodiscard]] long lastIteration() const noexcept { return iteration_; }
  [[nodiscard]] bool solidState() const noexcept { return solidState_; }
  [[nodiscard]] double solidJacobian() const noexcept { return solidJacobian_; }
  [[nodiscard]] const Matrix3d &cellForce(long image) const {
    return projectedCellForce.at(static_cast<size_t>(image));
  }
  // An explicit projection does not change the runtime CI flag.
  void updateForces(bool ci_active);
  void updateForces(void) { updateForces(ci_enabled_); }
  void setCIEnabled(bool enabled) {
    if (ci_enabled_ != enabled) {
      ci_enabled_ = enabled;
      // CI changes the projection while physical force caches remain valid.
      movedAfterForceCall = true;
    }
  }
  /// After a dimer has put this image on the saddle, band steps leave
  /// it there. A later dimer clears the pin before it walks again.
  void pinImage(long index) { pinnedImage_ = index; }
  void unpinImage() { pinnedImage_ = -1; }
  [[nodiscard]] long pinnedImage() const { return pinnedImage_; }
  double convergenceForce(void);
  /// Host bookkeeping for a band whose projected force came from the
  /// rgsaddle session: the highest interior image, the climbing image
  /// (`ciIndex`, or none when it is not an interior image) and the
  /// tangents, from the images' current energies and positions. The
  /// projected force stays as the session assembled it.
  void adoptSessionProjection(long ciIndex);
  void findExtrema(void);
  void printImageData(bool writeToFile = false, size_t idx = 0);

  /// In-memory ConFrames with the same NEB stamps as writePathCon / neb.con.
  /// Empty if path is incomplete. Does not write to disk.
  [[nodiscard]] std::vector<readcon::ConFrame>
  pathFrames(std::optional<size_t> bandIndex = std::nullopt);

  std::vector<std::shared_ptr<EigenmodeStrategy>> eigenmode_solvers;

  /// The identity band image i carries in a batch (forceBatchOwned):
  /// interior image i is i - 1, the reactant numImages and the product
  /// numImages + 1, so no image shares a key with an endpoint.
  [[nodiscard]] long imageOwner(std::size_t i) const noexcept {
    const auto idx = static_cast<long>(i);
    if (idx == 0)
      return numImages;
    return idx == numImages + 1 ? numImages + 1 : idx - 1;
  }

  int atoms{0};
  long numImages{0}, climbingImage{0}, numExtrema{0};
  // Reactant energy at construction. Zoom moves path[0] to the start of its
  // window, so reported energies are measured from this value instead.
  double reactantEnergy{0.0};
  std::vector<std::shared_ptr<Matter>> path; // NEB images
  std::vector<std::shared_ptr<AtomMatrix>> tangent;
  std::vector<std::shared_ptr<AtomMatrix>> projectedForce;
  std::vector<double> extremumEnergy;
  std::vector<double> extremumPosition;
  std::vector<double> extremumCurvature;

  std::size_t maxEnergyImage{0};
  // Geometry or runtime CI state differs from the cached band projection.
  bool movedAfterForceCall{false};
  bool perImagePotentials_{
      false}; ///< Whether per-image potential instances exist
  double ksp{0.0};
  double k_u{0.0}; // Upper-bound value for the spring constant
  double k_l{0.0}; // Lower-bound value for the spring constant
  double E_ref; // Reference energy chosen to be equal to the max energy of the
                // reactant or product energy minimum

private:
  bool ci_enabled_{false};
  // Image left fixed after a converged dimer. -1 pins nothing.
  long pinnedImage_{-1};
  double baseline_force{-1.0};
  Parameters params;
  std::shared_ptr<Potential> pot;
  NEBStatus status;
  long iteration_{0};
  eonc::log::Scoped log;

  // Cached strategies (constant across iterations)
  neb::TangentStrategy tangentStrat_;
  neb::ProjectionStrategy projectionStrat_;
  bool solidState_{false};
  double solidJacobian_{1.0};
  std::vector<Matrix3d> projectedCellForce;
  void prepareSolidState();
  void projectSolidState(bool ci_active);
};

class NEBObjectiveFunction : public ObjectiveFunction {
public:
  NEBObjectiveFunction(NudgedElasticBand *nebPassed,
                       const Parameters &parametersPassed)
      : ObjectiveFunction(parametersPassed), neb{nebPassed} {}
  // This is the odd one out, doesn't take a Matter so we null it

  ~NEBObjectiveFunction(void) {};

  VectorXd getGradient(bool fdstep = false);
  double getEnergy();
  void setPositions(const VectorXd &x);
  VectorXd getPositions();
  int degreesOfFreedom();
  bool isConverged();
  bool isUncertain();
  double getConvergence();
  VectorXd difference(const VectorXd &a, const VectorXd &b);
  // One mass per free atom on each interior image, in band order.
  VectorXd getMasses() const override;
  // The band residual rotates with its tangent and is not a PES gradient.
  bool supportsFiniteDifferenceCurvature() const override { return false; }
  // getEnergy() sums the image energies; the springs are not in it and the
  // climbing image moves up it.
  bool gradientMatchesEnergy() const override { return false; }
  bool reuseDisplacementCurvature() const override { return true; }
  NudgedElasticBand::NEBStatus status;

private:
  NudgedElasticBand *neb;
};

} // namespace eonc

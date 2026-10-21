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
#include "EonLogger.h"
#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

namespace eonc {

class BasinHoppingDisplaceAccess;

class BasinHoppingJob : public Job {
public:
  BasinHoppingJob(std::unique_ptr<Parameters> parameters, Runtime &rt)
      : Job(std::move(parameters), rt),
        current{std::make_shared<Matter>(pot, params)},
        trial{std::make_shared<Matter>(pot, params)}, fcalls{0} {}
  ~BasinHoppingJob(void) = default;

  std::vector<std::string> run(void) override;
  /// Metropolis weight for energy change `de`.
  /// Divides by `kB * temperature` only for an uphill hop when both are
  /// positive. Non-positive temperature or `kB` refuses that hop.
  static double metropolisProbability(double de, double kB, double temperature);
  /// Angle-axis map. `omega` is the rotation vector in radians.
  /// A zero vector returns `point`.
  static Vector3d rotateAngleAxis(const Vector3d &point, const Vector3d &omega);
  /// Column k is d(rotateAngleAxis)/d(omega_k).
  static Matrix3d rotateAngleAxisJacobian(const Vector3d &point,
                                         const Vector3d &omega);

  /// One rigid water molecule: oxygen and its two hydrogens.
  struct RigidBody {
    long atom[3]{};
  };

private:
  friend class BasinHoppingDisplaceAccess;
  VectorXd calculateDistanceFromCenter(Matter *matter);
  AtomMatrix displaceRandom(double maxDisplacement);
  void randomSwap(Matter *matter);
  /// True when the optimizer section is still the client default.
  /// Switches that default quench to L-BFGS.
  bool configureStockQuench();
  bool prepareRigidWater();
  void displaceRigid(Matter *matter, double translation, double angle);
  bool quench(Matter &matter, bool wantTight, double bestEnergy);
  bool quenchCartesian(Matter &matter, bool wantTight, double bestEnergy);
  bool quenchRigid(Matter &matter, bool wantTight, double bestEnergy);
  void adaptStockStep(double &step, int recentAccept, int period);
  void adaptRigidStep(double &step, int recentAccept, int period);
  void seedRigidTemperature(double energy);
  void seedAtomicTemperature(double energy);
  void liftTemperature(double deltaE);
  [[nodiscard]] double hopTemperature() const;
  [[nodiscard]] double swapProbability(const std::vector<long> &elements) const;
  [[nodiscard]] long jumpStepCount() const;
  void separateWaters(AtomMatrix &positions) const;
  /// Point-group projection. False when the geometry is already symmetric
  /// or the projection does not move the atoms.
  bool projectSymmetry(AtomMatrix &positions) const;
  /// Rotate atoms outside the median radius about a principal axis.
  bool rotateOuterShell(AtomMatrix &positions) const;
  /// One stock symmetry step written into `trial`. False leaves the
  /// caller to take an ordinary Cartesian displacement. `stale` is the
  /// number of steps since the best energy improved.
  bool stockSymmetryStep(int stale);
  /// True for a free cluster, including one stored in a large cell.
  [[nodiscard]] bool isolatedCluster() const;
  std::shared_ptr<Matter> current;
  std::shared_ptr<Matter> trial; // initial configuration
  std::vector<long> getElements(Matter *matter);
  friend struct BasinHoppingElementsTest;
  std::vector<std::string> returnFiles;
  int jump_count{0}; // count of jump moves
  int disp_count{0}; // count of displacement moves
  int swap_count{0}; // count of swap moves
  int symmetry_count{0}; // point-group and outer-shell steps
  int fcalls{0};
  bool stock_{false};
  bool rigidWater_{false};
  /// Recent acceptance is already at the target, so a clash must not
  /// keep raising the temperature.
  bool temperatureHeld_{false};
  double temperature_{300.0};
  double angularDisplacement_{0.5};
  double initialDisplacement_{0.5};
  double initialAngular_{0.5};
  /// kT seeded from the quenched well of a rigid water cluster.
  double rigidScale_{0.0};
  std::vector<RigidBody> rigidBodies_;

  std::vector<std::shared_ptr<Matter>> uniqueStructures;
  std::vector<double> uniqueEnergies;
  eonc::log::Scoped log;
};

} // namespace eonc

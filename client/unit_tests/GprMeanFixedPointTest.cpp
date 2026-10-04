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
// gpr_optim-ls8r gate: eonclient minimization on a frozen GP surface
// reaches the same fixed point as the gpr-internal band L-BFGS engine.
// The surface is a fitted-Cholesky GprMeanPotential; both drivers see
// the identical frozen model (the factor epoch must not move).

#include "catch2/catch_amalgamated.hpp"

#include <memory>

#include "eon/Eigen.h"
#include "eon/LBFGS.h"
#include "eon/ObjectiveFunction.h"
#include "eon/Parameters.h"

using eonc::LBFGS;
using eonc::Parameters;
using eonc::ParametersLoadAccess;

#include "gpr/ml/GaussianProcessRegression.h"
#include "gpr/neb/GPNEB.h"
#include "gpr/potentials/GprMeanPotential.h"

namespace {

gpr::Observation twoPoint() {
  gpr::Observation train;
  train.R.resize(2, 6);
  train.E.resize(1, 2);
  train.G.resize(2, 6);
  train.R << 8.98237316483057, 9.93723083577204, 7.89441632385049,
      7.65248322727496, 9.95590549457398, 7.87787958998366, 8.97856277303058,
      9.93211628067229, 7.89882761414426, 7.64888749663556, 9.95517051512886,
      7.87274215046670;
  train.E << 0.0, 370.498489473903e-6;
  train.G << 22.2174753484639e-3, -22.2402530044595e-3, 18.5353722102027e-3,
      -17.3837548235882e-3, 21.9137885404279e-3, -25.8161071576406e-3,
      21.4452456338877e-3, -33.1231900862804e-3, 28.3157675762483e-3,
      -15.8030003032378e-3, 20.8628464314946e-3, -38.3251019392285e-3;
  return train;
}

gpr::AtomsConfiguration twoAtomConfig() {
  gpr::AtomsConfiguration cfg;
  const gpr::Index_t n_atoms = 2;
  cfg.positions.resize(1, 3 * n_atoms);
  cfg.positions << 0.0, 0.0, 0.0, 0.0, 0.0, 3.5;
  cfg.is_frozen.resize(1, n_atoms);
  cfg.is_frozen.setConstant(0);
  cfg.id.resize(1, n_atoms);
  cfg.id << 0, 1;
  cfg.atomicNrs.resize(1, n_atoms);
  cfg.atomicNrs.setConstant(1);
  cfg.type.resize(1, n_atoms);
  cfg.type.setZero();
  cfg.atoms_mov.resize(n_atoms);
  cfg.atoms_mov.type.resize(1, n_atoms);
  cfg.atoms_mov.type.setZero();
  cfg.pairtype.resize(1, 1);
  cfg.pairtype(0, 0) = 0;
  cfg.n_pt = 1;
  cfg.atoms_mov.positions.resize(1, 6);
  cfg.atoms_mov.positions = cfg.positions;
  return cfg;
}

void fitCholesky(gpr::GaussianProcessRegression &model,
                 const gpr::Observation &train) {
  gpr::InputParameters params;
  gpr::AtomsConfiguration cfg = twoAtomConfig();
  cfg.atoms_mov.positions = train.R.topRows(1);
  model.initialize(params, cfg);
  model.setHssvdBackend(false);
  model.setDisableOutputNormalization(true);
  gpr::GPRSetup setup;
  setup.jitter_sigma2 = 0.0;
  model.setParameters(setup);
  model.getLikGaussian().setSigma2(1.0e-8);
  model.getConstantCovarianceFunction().setConstSigma2(1.0);
  model.getSexpAtCovarianceFunction().getLengthScaleRef().resize(1, 2);
  model.getSexpAtCovarianceFunction().getLengthScaleRef()(0, 0) = 1.0155;
  model.getSexpAtCovarianceFunction().getLengthScaleRef()(0, 1) = 0.03182;
  model.getSexpAtCovarianceFunction().setMagnSigma2(8.17e-6);
  model.getSexpAtCovarianceFunction().setConfInfo(cfg);
  model.setJitterSigma2(0.0);
  model.updateModelWithFullData(train);
}

constexpr long kAtoms = 2;
const int kAtomicNrs[kAtoms] = {1, 1};
const double kBox[9] = {25.0, 0.0, 0.0, 0.0, 25.0, 0.0, 0.0, 0.0, 25.0};

/// eonclient side: the frozen GP surface as an eOn ObjectiveFunction.
class GprSurfaceObjective final : public eonc::ObjectiveFunction {
  pot::GprMeanPotential &m_surface;
  VectorXd m_x;

public:
  GprSurfaceObjective(pot::GprMeanPotential &surface, const Parameters &params,
                      VectorXd start)
      : eonc::ObjectiveFunction(params),
        m_surface(surface),
        m_x(std::move(start)) {}

  double getEnergy() override {
    VectorXd F(m_x.size());
    double U = 0.0;
    m_surface.force(kAtoms, m_x.data(), kAtomicNrs, F.data(), &U, nullptr,
                    kBox);
    return U;
  }

  VectorXd getGradient(bool /*fdstep*/ = false) override {
    VectorXd F(m_x.size());
    double U = 0.0;
    m_surface.force(kAtoms, m_x.data(), kAtomicNrs, F.data(), &U, nullptr,
                    kBox);
    return -F;
  }

  void setPositions(const VectorXd &x) override { m_x = x; }
  VectorXd getPositions() override { return m_x; }
  int degreesOfFreedom() override { return static_cast<int>(m_x.size()); }
  bool isConverged() override {
    return getConvergence() < params.optimizer_options().converged_force;
  }
  double getConvergence() override { return getGradient().norm(); }
  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    return a - b;
  }
};

/// gpr side: the same surface behind the band engine's objective seam.
class FrozenSurfaceObjective final
    : public gpr::neb::detail::GlobalLBFGSObjective {
  pot::GprMeanPotential &surface_;
  Eigen::VectorXd x_;
  Eigen::VectorXd f_;

public:
  FrozenSurfaceObjective(pot::GprMeanPotential &surface, Eigen::VectorXd start)
      : surface_(surface),
        x_(std::move(start)),
        f_(Eigen::VectorXd::Zero(x_.size())) {}

  [[nodiscard]] Eigen::VectorXd positions() const override { return x_; }
  [[nodiscard]] Eigen::VectorXd force() const override { return f_; }
  void setPositions(const Eigen::VectorXd &positions) override {
    x_ = positions;
  }
  bool refresh(bool /*finite_difference_step*/) override {
    double U = 0.0;
    surface_.force(kAtoms, x_.data(), kAtomicNrs, f_.data(), &U, nullptr,
                   kBox);
    return f_.allFinite();
  }
  [[nodiscard]] Eigen::Index imageDof() const override { return x_.size(); }
};

} // namespace

TEST_CASE("eonclient LBFGS and the internal band engine share the frozen "
          "GP fixed point",
          "[gprd][gprmean][fixedpoint]") {
  gpr::GaussianProcessRegression model;
  const gpr::Observation train = twoPoint();
  fitCholesky(model, train);
  pot::GprMeanPotential surface(model);
  const auto epoch_before = surface.surfaceEpoch();

  // Start between the two training rows, slightly displaced.
  Eigen::VectorXd start = train.R.row(0).transpose();
  start(0) += 5.0e-4;
  start(3) -= 5.0e-4;

  Parameters params;
  auto &opt_o = ParametersLoadAccess::optimizer_options(params);
  opt_o.converged_force = 1e-9;
  opt_o.max_move = 0.05;
  opt_o.max_iterations = 2000;
  opt_o.lbfgs.memory = 20;
  opt_o.lbfgs.auto_scale = true;
  opt_o.lbfgs.inverse_curvature = 0.01;
  opt_o.lbfgs.angle_reset = true;
  opt_o.lbfgs.distance_reset = true;

  // eonclient path.
  auto objf =
      std::make_shared<GprSurfaceObjective>(surface, params, VectorXd(start));
  LBFGS eon_opt(objf, params);
  const int status = eon_opt.run(2000, opt_o.max_move);
  const VectorXd x_eon = objf->getPositions();
  REQUIRE(status == 1);

  // Internal band engine on the identical surface, identical start.
  FrozenSurfaceObjective frozen(surface, start);
  gpr::neb::detail::InternalLbfgsEngine engine;
  REQUIRE(frozen.refresh(false));
  int iters = 0;
  while (frozen.force().norm() > 1e-9 && iters < 2000) {
    if (!engine.step(frozen, opt_o.max_move)) {
      break;
    }
    REQUIRE(frozen.refresh(false));
    ++iters;
  }
  const Eigen::VectorXd x_int = frozen.positions();
  REQUIRE(frozen.force().norm() <= 1e-9);

  // Same fixed point: positions and energies agree, and the surface
  // never refit under either driver.
  REQUIRE((x_eon - x_int).cwiseAbs().maxCoeff() < 1e-4);
  GprSurfaceObjective probe_eon(surface, params, VectorXd(x_eon));
  GprSurfaceObjective probe_int(surface, params, VectorXd(x_int));
  REQUIRE(probe_eon.getEnergy() ==
          Catch::Approx(probe_int.getEnergy()).margin(1e-10));
  REQUIRE(surface.surfaceEpoch() == epoch_before);
}

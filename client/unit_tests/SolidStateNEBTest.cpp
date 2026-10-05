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
#include "eon/SolidStateNEB.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/NudgedElasticBand.h"

#include <cmath>
#include <memory>

namespace tests {
namespace {

static eonc::helpers::test::QuillTestLogger quillSetup;

// E = 0.5*kp*(sx-0.5)^2 + 0.5*kc*(Lx - L0 - gamma*4*sx*(1-sx))^2
// sx = x/Lx of atom 0. The preferred length is longest at sx = 0.5.
struct BreathPot final : eonc::Potential {
  double kp{20.0};
  double kc{5.0};
  double L0{5.0};
  double gamma{0.0};
  bool native{false};
  Matrix3d last{Matrix3d::Zero()};

  BreathPot()
      : Potential(eonc::PotType::LJ) {}

  [[nodiscard]] bool computesStress() const noexcept override { return native; }

  [[nodiscard]] Matrix3d cauchyStress() const override { return last; }

  void force(long nAtoms, const double *positions, const int * /*atomicNrs*/,
             double *forces, double *energy, double *variance,
             const double *box) override {
    if (nAtoms < 1) {
      throw std::invalid_argument("BreathPot needs one atom");
    }
    const Matrix3d cell = Matrix3d::Map(box);
    const Vector3d cartesian(positions[0], positions[1], positions[2]);
    const double sx = (cell.inverse().transpose() * cartesian)(0);
    const double length = cell(0, 0);
    const double g = 4.0 * sx * (1.0 - sx);
    const double u = sx - 0.5;
    const double delta = length - L0 - gamma * g;
    *energy = 0.5 * kp * u * u + 0.5 * kc * delta * delta;
    if (variance != nullptr) {
      *variance = 0.0;
    }
    for (long i = 0; i < nAtoms * 3; ++i) {
      forces[i] = 0.0;
    }
    const double dEdsx =
        kp * u + kc * delta * (-gamma) * 4.0 * (1.0 - 2.0 * sx);
    const Vector3d gradient = dEdsx * cell.inverse().col(0);
    forces[0] = -gradient(0);
    forces[1] = -gradient(1);
    forces[2] = -gradient(2);
    const double volume = std::abs(cell.determinant());
    last.setZero();
    last(0, 0) = (kc * delta * length) / volume;
  }
};

struct Band {
  Parameters params;
  std::shared_ptr<BreathPot> pot;
  std::shared_ptr<eonc::Matter> reactant;
  std::shared_ptr<eonc::Matter> product;

  explicit Band(bool solid, double gamma = 0.0, double pressure = 0.0)
      : pot{std::make_shared<BreathPot>()},
        reactant{nullptr},
        product{nullptr} {
    pot->gamma = gamma;
    ParametersLoadAccess::main_options(params).parallel = false;
    ParametersLoadAccess::optimizer_options(params).method = eonc::OptType::SD;
    ParametersLoadAccess::optimizer_options(params).max_move = 0.15;
    ParametersLoadAccess::optimizer_options(params).max_iterations = 200;
    ParametersLoadAccess::neb_options(params).opt_method = eonc::OptType::SD;
    ParametersLoadAccess::neb_options(params).image_count = 3;
    ParametersLoadAccess::neb_options(params).max_iterations = 80;
    ParametersLoadAccess::neb_options(params).force_tolerance = 0.05;
    ParametersLoadAccess::neb_options(params).climbing_image.enabled = false;
    ParametersLoadAccess::neb_options(params).solid_state.enabled = solid;
    ParametersLoadAccess::neb_options(params).solid_state.pressure = pressure;
    reactant = std::make_shared<eonc::Matter>(pot, params);
    product = std::make_shared<eonc::Matter>(pot, params);
    // 0.3 and 0.7 keep the short fractional arc through sx = 0.5.
    place(*reactant, 0.3);
    place(*product, 0.7);
  }

  void place(eonc::Matter &matter, double sx) const {
    matter.resize(2);
    matter.setAtomicNr(0, 1);
    matter.setAtomicNr(1, 1);
    Matrix3d cell = Matrix3d::Zero();
    cell(0, 0) = 5.0;
    cell(1, 1) = 4.0;
    cell(2, 2) = 4.0;
    matter.setCell(cell);
    matter.setPeriodic(true);
    AtomMatrix positions(2, 3);
    positions.setZero();
    positions(0, 0) = sx * 5.0;
    positions(0, 1) = 2.0;
    positions(0, 2) = 2.0;
    matter.setPositions(positions);
    matter.setFixed(1, 1);
  }

  std::unique_ptr<eonc::NudgedElasticBand> neb() const {
    return std::make_unique<eonc::NudgedElasticBand>(reactant, product, params,
                                                     pot);
  }
};

TEST_CASE("solid-state Jacobian matches the cell-volume scaling",
          "[neb][solid_state]") {
  const double volume = 64.0;
  const int nAtoms = 8;
  const double weight = 2.0;
  const double expected =
      std::cbrt(volume / static_cast<double>(nAtoms)) * std::sqrt(8.0) * weight;
  REQUIRE(eonc::neb::solidStateJacobian(volume, nAtoms, weight) ==
          Catch::Approx(expected));
  REQUIRE_THROWS_AS(eonc::neb::solidStateJacobian(0.0, 8, 1.0),
                    std::invalid_argument);
}

TEST_CASE("lower-triangular orientation keeps fractional coordinates",
          "[neb][solid_state]") {
  Matrix3d cell;
  cell << 0.0, 0.0, 5.0, 0.0, 4.0, 0.0, 3.0, 1.0, 0.0;
  AtomMatrix positions(2, 3);
  positions << 1.0, 2.0, 3.0, 0.5, 0.2, 0.1;
  const double volume = cell.determinant();
  const double distance = (positions.row(0) - positions.row(1)).norm();
  const AtomMatrix fractional = positions * cell.inverse();
  REQUIRE(eonc::neb::orientCellLowerTriangular(cell, positions));
  REQUIRE(cell(0, 1) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(cell(0, 2) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(cell(1, 2) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(cell(0, 0) > 0.0);
  REQUIRE(cell.determinant() == Catch::Approx(volume));
  REQUIRE((positions.row(0) - positions.row(1)).norm() ==
          Catch::Approx(distance));
  const AtomMatrix after = positions * cell.inverse();
  REQUIRE(after.isApprox(fractional, 1e-10));

  Matrix3d diagonal = Matrix3d::Zero();
  diagonal(0, 0) = 5.0;
  diagonal(1, 1) = 4.0;
  diagonal(2, 2) = 6.0;
  AtomMatrix held(2, 3);
  held << 1.0, 2.0, 3.0, 0.5, 0.2, 0.1;
  const AtomMatrix expect = held;
  REQUIRE(eonc::neb::orientCellLowerTriangular(diagonal, held));
  REQUIRE(held.isApprox(expect, 1e-12));
  REQUIRE(diagonal(0, 0) == Catch::Approx(5.0));
  REQUIRE(diagonal(1, 1) == Catch::Approx(4.0));
  REQUIRE(diagonal(2, 2) == Catch::Approx(6.0));
}

TEST_CASE("finite-difference stress matches the analytic Cauchy stress",
          "[neb][solid_state]") {
  auto pot = std::make_shared<BreathPot>();
  pot->gamma = 0.4;
  Parameters params;
  eonc::Matter image(pot, params);
  image.resize(1);
  image.setAtomicNr(0, 1);
  Matrix3d cell = Matrix3d::Zero();
  cell(0, 0) = 5.0;
  cell(1, 1) = 4.0;
  cell(2, 2) = 4.0;
  image.setCell(cell);
  image.setPeriodic(true);
  AtomMatrix positions(1, 3);
  positions << 2.5, 2.0, 2.0;
  image.setPositions(positions);
  image.getForces();
  const Matrix3d stress = eonc::neb::finiteDifferenceCauchyStress(image, 1e-6);
  REQUIRE(stress(0, 0) == Catch::Approx(pot->last(0, 0)).margin(1e-8));
  REQUIRE(stress(1, 1) == Catch::Approx(0.0).margin(1e-8));
  REQUIRE(stress(0, 1) == Catch::Approx(0.0).margin(1e-10));
}

namespace {

// BreathPot energies through the batch interface, counting the rounds.
struct BatchedBreath final : eonc::Potential {
  BreathPot inner;
  long rounds{0};
  long systems{0};

  BatchedBreath()
      : Potential(eonc::PotType::LJ) {}

  using Potential::force;
  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    rounds++;
    systems++;
    inner.force(nAtoms, positions, atomicNrs, forces, energy, variance, box);
  }

  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return true;
  }

  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const *atomicNrs, double *const *forces,
                  double *energies, double *variances,
                  const double *const *boxes) override {
    rounds++;
    systems += nSystems;
    for (long s = 0; s < nSystems; ++s) {
      double var = 0.0;
      inner.force(nAtoms, positions[s], atomicNrs[s], forces[s], &energies[s],
                  &var, boxes[s]);
      if (variances != nullptr) {
        variances[s] = var;
      }
    }
  }
};

} // namespace

TEST_CASE("finite-difference stress strains a whole band in one batch",
          "[neb][solid_state][batch]") {
  auto batched = std::make_shared<BatchedBreath>();
  batched->inner.gamma = 0.4;
  auto serial = std::make_shared<BreathPot>();
  serial->gamma = 0.4;
  Parameters params;
  std::vector<eonc::Matter> images;
  for (double sx : {0.35, 0.5, 0.65}) {
    eonc::Matter &image = images.emplace_back(batched, params);
    image.resize(1);
    image.setAtomicNr(0, 1);
    Matrix3d cell = Matrix3d::Zero();
    cell(0, 0) = 5.0 + sx;
    cell(1, 1) = 4.0;
    cell(2, 2) = 4.0;
    image.setCell(cell);
    image.setPeriodic(true);
    AtomMatrix positions(1, 3);
    positions << sx * cell(0, 0), 2.0, 2.0;
    image.setPositions(positions);
  }
  std::vector<const eonc::Matter *> band;
  for (const auto &image : images) {
    band.push_back(&image);
  }
  const auto stresses = eonc::neb::finiteDifferenceCauchyStresses(band, 1e-6);
  REQUIRE(batched->rounds == 1);
  REQUIRE(batched->systems == 36);
  REQUIRE(stresses.size() == 3);
  for (size_t k = 0; k < images.size(); ++k) {
    eonc::Matter alone(images[k]);
    alone.setPotential(serial);
    const Matrix3d one = eonc::neb::finiteDifferenceCauchyStress(alone, 1e-6);
    REQUIRE((stresses[k] - one).cwiseAbs().maxCoeff() < 1e-12);
  }
}

TEST_CASE("an unchanged cell leaves the solid-state atomic force unchanged",
          "[neb][solid_state]") {
  Band plain(false);
  Band solid(true);
  auto plainNeb = plain.neb();
  auto solidNeb = solid.neb();
  plainNeb->updateForces();
  solidNeb->updateForces();
  REQUIRE(solidNeb->solidState());
  eonc::NEBObjectiveFunction objective(solidNeb.get(), solid.params);
  REQUIRE(objective.degreesOfFreedom() == 3 * (3 * 2 + 9));
  for (long i = 1; i <= solidNeb->numImages; ++i) {
    REQUIRE(solidNeb->cellForce(i).norm() == Catch::Approx(0.0).margin(1e-8));
    const AtomMatrix difference =
        *solidNeb->projectedForce[static_cast<size_t>(i)] -
        *plainNeb->projectedForce[static_cast<size_t>(i)];
    // Near-zero forces make Eigen's relative isApprox reject a 1e-14 gap.
    REQUIRE(difference.norm() == Catch::Approx(0.0).margin(1e-8));
  }
}

TEST_CASE("a native stress comes with the band's own force call",
          "[neb][solid_state]") {
  Band analytic(true, 0.4);
  analytic.pot->native = true;
  Band difference(true, 0.4);
  auto nativeNeb = analytic.neb();
  auto fdNeb = difference.neb();
  size_t dirty = 0;
  for (long i = 1; i <= nativeNeb->numImages; ++i) {
    dirty += nativeNeb->path[static_cast<size_t>(i)]->needsForceUpdate();
  }
  REQUIRE(dirty > 0);
  const size_t before = analytic.pot->forceCallCounter;
  nativeNeb->updateForces();
  // One call per moved image: the stress the potential reported with it
  // is the one the cell force uses.
  REQUIRE(analytic.pot->forceCallCounter - before == dirty);
  fdNeb->updateForces();
  for (long i = 1; i <= nativeNeb->numImages; ++i) {
    REQUIRE(nativeNeb->cellForce(i)(0, 0) ==
            Catch::Approx(fdNeb->cellForce(i)(0, 0)).margin(1e-6));
  }
  // A moved image drops the cached stress with its forces.
  AtomMatrix moved = nativeNeb->path[2]->getPositions();
  moved(0, 0) += 0.05;
  nativeNeb->path[2]->setPositions(moved);
  const Matrix3d fresh = nativeNeb->path[2]->cauchyStress();
  REQUIRE(fresh(0, 0) == Catch::Approx(analytic.pot->last(0, 0)));
}

TEST_CASE("positive pressure pushes the cell inward", "[neb][solid_state]") {
  Band solid(true, 0.0, 0.02);
  auto neb = solid.neb();
  neb->updateForces();
  REQUIRE(neb->cellForce(2)(0, 0) < 0.0);
}

TEST_CASE("the interior cell lengthens where the path prefers it",
          "[neb][solid_state]") {
  Band solid(true, 0.4);
  auto neb = solid.neb();
  const double reactantL = neb->path.front()->getCell()(0, 0);
  const double productL = neb->path.back()->getCell()(0, 0);
  neb->updateForces();
  REQUIRE(neb->cellForce(2)(0, 0) > 0.0);
  const auto status = neb->compute();
  REQUIRE(status == eonc::NudgedElasticBand::NEBStatus::GOOD);
  REQUIRE(neb->path[2]->getCell()(0, 0) > reactantL + 0.05);
  REQUIRE(neb->path.front()->getCell()(0, 0) == Catch::Approx(reactantL));
  REQUIRE(neb->path.back()->getCell()(0, 0) == Catch::Approx(productL));
  for (long i = 1; i <= neb->numImages; ++i) {
    const Matrix3d cell = neb->path[static_cast<size_t>(i)]->getCell();
    REQUIRE(cell(0, 1) == Catch::Approx(0.0).margin(1e-8));
    REQUIRE(cell(0, 2) == Catch::Approx(0.0).margin(1e-8));
    REQUIRE(cell(1, 2) == Catch::Approx(0.0).margin(1e-8));
    REQUIRE(cell(1, 1) == Catch::Approx(4.0).margin(1e-6));
    REQUIRE(cell(2, 2) == Catch::Approx(4.0).margin(1e-6));
  }
}

TEST_CASE("solid_state refuses a nonperiodic band and OCINEB",
          "[neb][solid_state]") {
  Band solid(true);
  solid.reactant->setPeriodic(false);
  solid.product->setPeriodic(false);
  REQUIRE_THROWS_AS(solid.neb(), std::invalid_argument);

  Band walk(true);
  ParametersLoadAccess::neb_options(walk.params).climbing_image.ocineb.use_mmf =
      true;
  REQUIRE_THROWS_AS(walk.neb(), std::invalid_argument);

  Band idpp(true);
  ParametersLoadAccess::neb_options(idpp.params).initialization.method =
      eonc::NEBInit::IDPP;
  REQUIRE_THROWS_AS(idpp.neb(), std::invalid_argument);
}

TEST_CASE("solid_state ini keys set the lattice band", "[neb][solid_state]") {
  Parameters params;
  REQUIRE(params.load_ini_text("[Nudged Elastic Band]\n"
                               "solid_state = true\n"
                               "solid_state_weight = 2.5\n"
                               "solid_state_pressure = 0.01\n") == 0);
  REQUIRE(params.neb_options().solid_state.enabled);
  REQUIRE(params.neb_options().solid_state.weight == Catch::Approx(2.5));
  REQUIRE(params.neb_options().solid_state.pressure == Catch::Approx(0.01));
}

} // namespace
} // namespace tests

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
// A solid-state band of shifted Lennard-Jones argon from rgpot's LJ kernel.
//
// rgpot's pair kernels sum one minimum image per pair. In a cell narrower
// than twice the cutoff the energy jumps and the force is not its gradient,
// so FIRE sees P < 0 on every step, halves its time step below 1e-6 and
// stops with "[FIRE] m_dt is too small". L-BFGS on the same band runs to
// max_iterations. The band now refuses such a cell up front.
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/NudgedElasticBand.h"
#include "eon/SolidStateNEB.h"
#include "eon/potentials/RgpotAdapter/RgpotAdapter.h"
#include "rgpot/LennardJones/LJPot.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

namespace tests {
namespace {

static eonc::helpers::test::QuillTestLogger quillSetup;

// n x n x n conventional fcc cells of argon in the cell `cell`, placed at
// fixed fractional coordinates.
void fillFcc(eonc::Matter &matter, int n, const Matrix3d &cell) {
  const double basis[4][3] = {
      {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};
  const long count = 4L * n * n * n;
  matter.resize(count);
  AtomMatrix fractional(count, 3);
  long atom = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        for (const auto &site : basis) {
          fractional(atom, 0) = (site[0] + i) / n;
          fractional(atom, 1) = (site[1] + j) / n;
          fractional(atom, 2) = (site[2] + k) / n;
          matter.setAtomicNr(atom, 18);
          matter.setMass(atom, 39.948);
          ++atom;
        }
      }
    }
  }
  matter.setCell(cell);
  matter.setPeriodic(true);
  matter.setPositions(fractional * cell);
}

struct Band {
  Parameters params;
  std::shared_ptr<RgpotAdapter<rgpot::LJPot>> pot;
  std::shared_ptr<eonc::Matter> reactant;
  std::shared_ptr<eonc::Matter> product;

  // Reactant: cubic cell of side n a. Product: the cell stretched by 4 % in
  // x and y and compressed to match in z, at fixed fractional coordinates.
  Band(eonc::OptType method, int n, double lattice, double cutoff) {
    ParametersLoadAccess::main_options(params).parallel = false;
    auto &neb = ParametersLoadAccess::neb_options(params);
    neb.image_count = 5;
    neb.max_iterations = 2000;
    neb.opt_method = method;
    neb.force_tolerance = 1e-3;
    neb.solid_state.enabled = true;
    neb.climbing_image.enabled = true;
    neb.climbing_image.trigger_force = 0.5;
    auto &opt = ParametersLoadAccess::optimizer_options(params);
    opt.method = method;
    opt.max_move = 0.05;
    opt.max_iterations = 2000;
    opt.converged_force = neb.force_tolerance;

    pot = std::make_shared<RgpotAdapter<rgpot::LJPot>>(
        eonc::PotType::LJ, params, rgpot::LJConfig{.cutoff = cutoff});
    reactant = std::make_shared<eonc::Matter>(pot, params);
    product = std::make_shared<eonc::Matter>(pot, params);
    const double side = n * lattice;
    const double stretch = 1.04;
    Matrix3d h0 = Matrix3d::Zero();
    Matrix3d h1 = Matrix3d::Zero();
    h0.diagonal() << side, side, side;
    h1.diagonal() << side * stretch, side * stretch, side / (stretch * stretch);
    fillFcc(*reactant, n, h0);
    fillFcc(*product, n, h1);
  }

  std::unique_ptr<eonc::NudgedElasticBand> neb() const {
    return std::make_unique<eonc::NudgedElasticBand>(reactant, product, params,
                                                     pot);
  }
};

TEST_CASE("the perpendicular width of a cell is its closest face spacing",
          "[neb][solid_state]") {
  Matrix3d box = Matrix3d::Zero();
  box.diagonal() << 3.0, 2.0, 5.0;
  REQUIRE(eonc::neb::minimumPerpendicularWidth(box) ==
          Catch::Approx(2.0).margin(1e-12));

  // Rows are the lattice vectors. Shearing b along a leaves the spacing of
  // the faces of b at b_y.
  Matrix3d sheared = box;
  sheared(1, 0) = 1.5;
  REQUIRE(eonc::neb::minimumPerpendicularWidth(sheared) ==
          Catch::Approx(2.0).margin(1e-12));

  // c leaning over a with |b| = 3: the faces spanned by b and c are
  // V / |b x c| = 45 / |(15, 0, -12)| apart, closer than the other pairs.
  Matrix3d leaning = box;
  leaning(1, 1) = 3.0;
  leaning(2, 0) = 4.0;
  REQUIRE(eonc::neb::minimumPerpendicularWidth(leaning) ==
          Catch::Approx(45.0 / std::sqrt(369.0)).margin(1e-12));

  // Each width is the reciprocal of a column norm of the inverse cell.
  Matrix3d general;
  general << 3.1, 0.0, 0.0, 0.7, 2.9, 0.0, -0.4, 0.6, 3.3;
  const Matrix3d inverse = general.inverse();
  const double expected = 1.0 / inverse.colwise().norm().maxCoeff();
  REQUIRE(eonc::neb::minimumPerpendicularWidth(general) ==
          Catch::Approx(expected).margin(1e-12));
}

TEST_CASE("a cutoff of half the width or less is accepted",
          "[neb][solid_state]") {
  Matrix3d box = Matrix3d::Zero();
  box.diagonal() << 4.0, 4.0, 4.0;
  REQUIRE_NOTHROW(eonc::neb::requireMinimumImageCell(box, 2.0, "cell"));
  REQUIRE_NOTHROW(eonc::neb::requireMinimumImageCell(box, 0.0, "cell"));
  REQUIRE_THROWS_AS(eonc::neb::requireMinimumImageCell(box, 2.001, "cell"),
                    std::invalid_argument);
}

TEST_CASE("the rgpot Lennard-Jones adapter reports its minimum-image cutoff",
          "[neb][solid_state]") {
  Parameters params;
  RgpotAdapter<rgpot::LJPot> pot(eonc::PotType::LJ, params,
                                 rgpot::LJConfig{.cutoff = 3.25});
  REQUIRE(pot.minimumImageCutoff() == 3.25);
}

// The cell of four atoms is 1.6 Angstrom wide and the cutoff is 1.2. Every
// nearest-neighbor pair sits on a minimum-image tie, the force on an atom
// is not the derivative of the energy, and the FIRE time step collapses
// within a few dozen iterations. The refusal replaces that abort.
TEST_CASE("a solid-state band refuses a cell narrower than twice the cutoff",
          "[neb][solid_state][fire]") {
  for (const auto method : {eonc::OptType::FIRE, eonc::OptType::LBFGS}) {
    Band band(method, 1, 1.6, 1.2);
    REQUIRE_THROWS_WITH(band.neb(),
                        Catch::Matchers::ContainsSubstring(
                            "must be at least twice the cutoff wide"));
  }
}

TEST_CASE("the solid-state band converges with FIRE and with L-BFGS in a "
          "cell wider than twice the cutoff",
          "[neb][solid_state][fire]") {
  for (const auto method : {eonc::OptType::FIRE, eonc::OptType::LBFGS}) {
    Band setup(method, 2, 1.6, 1.4);
    auto band = setup.neb();
    const auto status = band->compute();
    const double force = band->convergenceForce();
    const std::string name = method == eonc::OptType::FIRE ? "FIRE" : "LBFGS";
    std::printf("%s solid-state band: status %d iterations %ld max force "
                "%.3e\n",
                name.c_str(), static_cast<int>(status), band->lastIteration(),
                force);
    REQUIRE(status == eonc::NudgedElasticBand::NEBStatus::GOOD);
    REQUIRE(band->lastIteration() > 0);
    REQUIRE(force < 1e-3);
  }
}

} // namespace
} // namespace tests

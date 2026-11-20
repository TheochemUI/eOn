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
// rgpot's periodic pair kernels (Lennard-Jones, Morse) sum one minimum
// image per pair. Unless twice the cutoff is strictly inside the cell, a
// half-cell tie is an energy cusp: the force jumps and is not the energy
// gradient, so FIRE sees P < 0 on every step, halves its time step below
// 1e-6 and stops with "[FIRE] m_dt is too small". L-BFGS on the same band
// runs to max_iterations. The band now refuses such a cell up front.
// Nonperiodic ZBL reports no minimum-image cutoff and is not checked.
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/NudgedElasticBand.h"
#include "eon/SolidStateNEB.h"
#include "eon/potentials/RgpotAdapter/RgpotAdapter.h"
#include "rgpot/LennardJones/LJPot.hpp"
#include "rgpot/Morse/MorsePot.hpp"
#include "rgpot/ZBL/ZBLPot.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
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

TEST_CASE("a zero or singular cell has perpendicular width zero",
          "[neb][solid_state]") {
  const double zero = eonc::neb::minimumPerpendicularWidth(Matrix3d::Zero());
  REQUIRE(std::isfinite(zero));
  REQUIRE(zero == 0.0);

  // Rank 2: the third row is the sum of the first two, so the volume is
  // zero and every face product that would divide it is zero as well.
  Matrix3d singular;
  singular << 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 2.0, 0.0;
  REQUIRE(singular.determinant() == 0.0);
  const double collapsed = eonc::neb::minimumPerpendicularWidth(singular);
  REQUIRE(std::isfinite(collapsed));
  REQUIRE(collapsed == 0.0);

  REQUIRE_THROWS_AS(
      eonc::neb::requireMinimumImageCell(Matrix3d::Zero(), 1.2, "cell"),
      std::invalid_argument);
  REQUIRE_THROWS_AS(eonc::neb::requireMinimumImageCell(singular, 1.2, "cell"),
                    std::invalid_argument);
}

TEST_CASE("only a cutoff strictly under half the width is accepted",
          "[neb][solid_state]") {
  Matrix3d box = Matrix3d::Zero();
  box.diagonal() << 4.0, 4.0, 4.0;
  // 2 and 4 are dyadic, so 2 * 2.0 == 4.0 with no rounding. Equality is
  // a tie and is refused; the next cutoff below 2 is accepted.
  REQUIRE_THROWS_AS(eonc::neb::requireMinimumImageCell(box, 2.0, "cell"),
                    std::invalid_argument);
  REQUIRE_NOTHROW(eonc::neb::requireMinimumImageCell(box, 0.0, "cell"));
  REQUIRE_NOTHROW(eonc::neb::requireMinimumImageCell(
      box, std::nextafter(2.0, 0.0), "cell"));
  REQUIRE_THROWS_AS(eonc::neb::requireMinimumImageCell(box, 2.001, "cell"),
                    std::invalid_argument);
}

// LJ cutoff 1.2 in a cube of side 2.4. 2 * 1.2 is exactly 2.4 in double,
// and the pair list includes r^2 <= cutoff^2, so this boundary is a tie.
TEST_CASE("a Lennard-Jones cell exactly twice the cutoff wide is refused",
          "[neb][solid_state]") {
  const double cutoff = 1.2;
  REQUIRE(2.0 * cutoff == 2.4);
  Matrix3d tie = Matrix3d::Zero();
  tie.diagonal().setConstant(2.4);
  REQUIRE_THROWS_AS(eonc::neb::requireMinimumImageCell(tie, cutoff, "cell"),
                    std::invalid_argument);

  const double justAbove =
      std::nextafter(2.0 * cutoff, std::numeric_limits<double>::infinity());
  Matrix3d wider = Matrix3d::Zero();
  wider.diagonal().setConstant(justAbove);
  REQUIRE(2.0 * cutoff < justAbove);
  REQUIRE_NOTHROW(eonc::neb::requireMinimumImageCell(wider, cutoff, "cell"));
}

TEST_CASE("the rgpot Lennard-Jones adapter reports its minimum-image cutoff",
          "[neb][solid_state]") {
  Parameters params;
  RgpotAdapter<rgpot::LJPot> pot(eonc::PotType::LJ, params,
                                 rgpot::LJConfig{.cutoff = 3.25});
  REQUIRE(pot.minimumImageCutoff() == 3.25);
}

TEST_CASE("the rgpot Morse adapter reports its configured cutoff",
          "[neb][solid_state]") {
  Parameters params;
  constexpr double cutoff = 5.0;
  RgpotAdapter<rgpot::MorsePot> pot(eonc::PotType::MORSE_PT, params,
                                    rgpot::MorseConfig{.cutoff = cutoff});
  REQUIRE(pot.minimumImageCutoff() == cutoff);
}

TEST_CASE("the rgpot ZBL adapter reports no minimum-image cutoff",
          "[neb][solid_state]") {
  Parameters params;
  // cut_global is the pair truncation, and it is not a minimum-image
  // cutoff: ZBL is nonperiodic and ignores the cell.
  RgpotAdapter<rgpot::ZBLPot> pot(
      eonc::PotType::ZBL, params,
      rgpot::ZBLConfig{.cut_inner = 1.0, .cut_global = 2.5});
  REQUIRE(pot.hessianRange() == 2.5);
  REQUIRE(pot.minimumImageCutoff() == 0.0);
}

// The cell of four atoms is 1.6 Angstrom wide and the cutoff is 1.2. Every
// nearest-neighbor pair sits on a minimum-image tie, the force on an atom
// is not the derivative of the energy, and the FIRE time step collapses
// within a few dozen iterations. The refusal replaces that abort.
TEST_CASE("a solid-state band refuses a cell narrower than twice the cutoff",
          "[neb][solid_state][fire]") {
  for (const auto method : {eonc::OptType::FIRE, eonc::OptType::LBFGS}) {
    Band band(method, 1, 1.6, 1.2);
    REQUIRE_THROWS_WITH(band.neb(), Catch::Matchers::ContainsSubstring(
                                        "must be wider than twice the cutoff"));
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

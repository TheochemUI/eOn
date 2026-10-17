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
// Bain path of a Lennard-Jones crystal on the solid-state band.
//
// LJ with u0 1 eV and psi 1 Angstrom, cut at 2.2 Angstrom with the C^2
// switch over the last 0.2 Angstrom. bcc and fcc are the two-atom
// body-centred tetragonal cell diag(a, a, c) repeated 5 x 5 x 4 (200
// atoms); every width stays above twice the cutoff, as rgpot's
// minimum-image pair list needs. Every atom sits on an inversion centre,
// so the transformation is the cell strain alone.
//
// Reference: validation/bain_lj_reference.py in rgsaddle, a relaxed scan
// of the energy per atom over c/a with spacing 2.5e-4:
//   bcc     a = 1.228441371,                  E = -7.4148256256 eV/atom
//   saddle  c/a = 1.09331985,                 E = -7.3557899883 eV/atom
//   fcc     a = 1.100691058, c = 1.556612222, E = -7.6698410342 eV/atom
// The climbing image must reach the scan's c/a within the spacing, and
// its barrier within kBarrierTol of 0.0590356373 eV/atom.
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/NudgedElasticBand.h"
#include "eon/Potential.h"

#include <cmath>
#include <memory>

namespace tests {
namespace {

static eonc::helpers::test::QuillTestLogger quillSetup;

constexpr int kNx = 5, kNy = 5, kNz = 4;
constexpr long kAtoms = 2L * kNx * kNy * kNz;
constexpr double kBccA = 1.228441371;
constexpr double kFccA = 1.100691058;
constexpr double kFccC = 1.556612222;
constexpr double kBccEnergy = -7.4148256256;
constexpr double kSaddleQ = 1.09331985;
constexpr double kBarrier = 0.0590356373;
constexpr double kScanSpacing = 2.5e-4;
// The energy change over one scan step at the saddle is 1.2e-6 eV/atom.
constexpr double kBarrierTol = 1.2e-6;

void place(eonc::Matter &matter, double a, double c) {
  matter.resize(kAtoms);
  Matrix3d cell = Matrix3d::Zero();
  cell(0, 0) = kNx * a;
  cell(1, 1) = kNy * a;
  cell(2, 2) = kNz * c;
  matter.setCell(cell);
  matter.setPeriodic(true);
  AtomMatrix positions(kAtoms, 3);
  long n = 0;
  for (int i = 0; i < kNx; ++i) {
    for (int j = 0; j < kNy; ++j) {
      for (int k = 0; k < kNz; ++k) {
        for (int b = 0; b < 2; ++b) {
          matter.setAtomicNr(n, 18);
          positions(n, 0) = (i + 0.5 * b) * a;
          positions(n, 1) = (j + 0.5 * b) * a;
          positions(n, 2) = (k + 0.5 * b) * c;
          ++n;
        }
      }
    }
  }
  matter.setPositions(positions);
}

TEST_CASE("the solid-state band climbs the Bain path of an LJ crystal",
          "[neb][solid_state][bain]") {
  Parameters params;
  REQUIRE(params.load_ini_text("[Potential]\n"
                               "potential = lj\n"
                               "pair_cutoff = 2.2\n"
                               "pair_switch_width = 0.2\n") == 0);
  ParametersLoadAccess::main_options(params).parallel = false;
  auto &neb = ParametersLoadAccess::neb_options(params);
  neb.image_count = 7;
  neb.max_iterations = 5000;
  neb.opt_method = eonc::OptType::FIRE;
  neb.force_tolerance = 1e-3;
  neb.solid_state.enabled = true;
  neb.climbing_image.enabled = true;
  // Arm the climbing image once the band force is under 0.5 eV/Angstrom
  // (ci_after), rgsaddle's default.
  neb.climbing_image.trigger_force = 0.5;
  auto &opt = ParametersLoadAccess::optimizer_options(params);
  opt.method = eonc::OptType::FIRE;
  opt.max_move = 0.1;
  opt.max_iterations = 5000;
  opt.converged_force = 1e-3;

  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = std::make_shared<eonc::Matter>(pot, params);
  auto product = std::make_shared<eonc::Matter>(pot, params);
  place(*reactant, kBccA, kBccA);
  place(*product, kFccA, kFccC);
  REQUIRE_THAT(reactant->getPotentialEnergy() / kAtoms,
               Catch::Matchers::WithinAbs(kBccEnergy, 1e-9));

  eonc::NudgedElasticBand band(reactant, product, params, pot);
  const auto status = band.compute();
  REQUIRE(status == eonc::NudgedElasticBand::NEBStatus::GOOD);

  const long ci = band.climbingImage;
  REQUIRE(ci >= 1);
  REQUIRE(ci <= band.numImages);
  const double e0 = band.path.front()->getPotentialEnergy();
  for (long i = 0; i <= band.numImages + 1; ++i) {
    const Matrix3d h = band.path[static_cast<size_t>(i)]->getCell();
    const double q = (h(2, 2) / kNz) / (h(0, 0) / kNx);
    const double de =
        (band.path[static_cast<size_t>(i)]->getPotentialEnergy() - e0) / kAtoms;
    UNSCOPED_INFO("image " << i << " c/a " << q << " E - E_bcc " << de);
    // No rotation, no shear drift, a square base.
    REQUIRE(h(0, 1) == 0.0);
    REQUIRE(h(0, 2) == 0.0);
    REQUIRE(h(1, 2) == 0.0);
    REQUIRE(std::abs(h(1, 0)) < 1e-10);
    REQUIRE(std::abs(h(2, 0)) < 1e-10);
    REQUIRE(std::abs(h(2, 1)) < 1e-10);
    REQUIRE(std::abs(h(0, 0) / kNx - h(1, 1) / kNy) < 1e-10);
  }
  const Matrix3d hci = band.path[static_cast<size_t>(ci)]->getCell();
  const double q = (hci(2, 2) / kNz) / (hci(0, 0) / kNx);
  const double barrier =
      (band.path[static_cast<size_t>(ci)]->getPotentialEnergy() - e0) /
      kAtoms;
  INFO("climbing image " << ci << " c/a " << q << " barrier " << barrier
                         << " eV/atom after " << band.lastIteration()
                         << " iterations, " << pot->forceCallCounter.load()
                         << " force calls");
  REQUIRE(std::abs(q - kSaddleQ) <= kScanSpacing);
  REQUIRE(std::abs(barrier - kBarrier) <= kBarrierTol);

  // No atom leaves its site: up to the lattice translation that wrapping
  // into the cell adds, the fractional coordinates of every image are
  // those of the reactant.
  const AtomMatrix s0 =
      reactant->getPositions() * reactant->getCell().inverse();
  for (long i = 1; i <= band.numImages; ++i) {
    const auto &image = *band.path[static_cast<size_t>(i)];
    const AtomMatrix ds =
        image.getPositions() * image.getCell().inverse() - s0;
    REQUIRE((ds.array() - ds.array().round()).abs().maxCoeff() < 1e-6);
  }
}

} // namespace
} // namespace tests

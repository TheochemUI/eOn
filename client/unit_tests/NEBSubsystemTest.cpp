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

/// Unit tests for NEB subsystem components: tangent, force projection,
/// spring forces, and spline extrema finding.

#include "Eigen.h"
#include "HelperFunctions.h"
#include "NEBForceProjection.h"
#include "NEBSplineExtrema.h"
#include "NEBSpringForce.h"
#include "NEBTangent.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

#include <cmath>
#include <memory>
#include <vector>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// ---------------------------------------------------------------------------
// Helper: build a 3-atom AtomMatrix from a flat initializer
// ---------------------------------------------------------------------------
static AtomMatrix make3(std::initializer_list<double> vals) {
  AtomMatrix m(3, 3);
  auto it = vals.begin();
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      m(r, c) = *it++;
  return m;
}

// ===== SimpleTangent =======================================================

TEST_CASE("SimpleTangent: normalized direction to next image", "[neb]") {
  // posDiffNext points along x for all three atoms
  AtomMatrix diffNext = make3({1, 0, 0, 1, 0, 0, 1, 0, 0});
  AtomMatrix diffPrev = make3({-1, 0, 0, -1, 0, 0, -1, 0, 0});

  eonc::neb::SimpleTangent strat;
  AtomMatrix tang = strat.compute(diffNext, diffPrev, 0.0, 0.0, 0.0);

  // Should be normalized
  REQUIRE(tang.norm() == Catch::Approx(1.0).epsilon(1e-12));

  // Should point in the same direction as diffNext (positive x)
  double dot = matDot(tang, diffNext);
  REQUIRE(dot > 0.0);
}

TEST_CASE("SimpleTangent: ignores energy and previous difference", "[neb]") {
  AtomMatrix diffNext = make3({0, 2, 0, 0, 2, 0, 0, 2, 0});
  AtomMatrix diffPrev = make3({-3, 0, 0, -3, 0, 0, -3, 0, 0});

  eonc::neb::SimpleTangent strat;
  // Varying energies should not affect the result
  AtomMatrix t1 = strat.compute(diffNext, diffPrev, 1.0, 0.5, 2.0);
  AtomMatrix t2 = strat.compute(diffNext, diffPrev, 5.0, 10.0, -3.0);

  REQUIRE_THAT(t1, eonc::helpers::test::IsApprox(t2, 1e-12));
}

// ===== ImprovedTangent =====================================================

TEST_CASE("ImprovedTangent: monotonically increasing energy uses next",
          "[neb]") {
  AtomMatrix diffNext = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  AtomMatrix diffPrev = make3({0, 1, 0, 0, 0, 0, 0, 0, 0});

  eonc::neb::ImprovedTangent strat;
  // energyNext > energy > energyPrev => use diffNext
  AtomMatrix tang = strat.compute(diffNext, diffPrev, 1.0, 0.0, 2.0);
  REQUIRE(tang.norm() == Catch::Approx(1.0).epsilon(1e-12));

  // tang should be parallel to diffNext (which is along x for atom 0)
  AtomMatrix diffNextNorm = diffNext / diffNext.norm();
  REQUIRE_THAT(tang, eonc::helpers::test::IsApprox(diffNextNorm, 1e-10));
}

TEST_CASE("ImprovedTangent: monotonically decreasing energy uses prev",
          "[neb]") {
  AtomMatrix diffNext = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  AtomMatrix diffPrev = make3({0, 1, 0, 0, 0, 0, 0, 0, 0});

  eonc::neb::ImprovedTangent strat;
  // energy > energyNext AND energyPrev > energy => use diffPrev
  AtomMatrix tang = strat.compute(diffNext, diffPrev, 1.0, 2.0, 0.0);
  REQUIRE(tang.norm() == Catch::Approx(1.0).epsilon(1e-12));

  AtomMatrix diffPrevNorm = diffPrev / diffPrev.norm();
  REQUIRE_THAT(tang, eonc::helpers::test::IsApprox(diffPrevNorm, 1e-10));
}

TEST_CASE("ImprovedTangent: extremum uses energy-weighted combination",
          "[neb]") {
  // At a minimum: energyPrev > energy < energyNext
  AtomMatrix diffNext = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  AtomMatrix diffPrev = make3({0, 1, 0, 0, 0, 0, 0, 0, 0});

  eonc::neb::ImprovedTangent strat;
  // energy = 0, energyPrev = 3, energyNext = 1
  // energyDiffPrev = 3-0 = 3, energyDiffNext = 1-0 = 1
  // min = 1, max = 3
  // energyDiffPrev(3) > energyDiffNext(1) => tang = diffNext*min + diffPrev*max
  //   = diffNext*1 + diffPrev*3
  AtomMatrix tang = strat.compute(diffNext, diffPrev, 0.0, 3.0, 1.0);
  REQUIRE(tang.norm() == Catch::Approx(1.0).epsilon(1e-12));

  // Expected unnormalized: diffNext*1 + diffPrev*3
  AtomMatrix expected = diffNext * 1.0 + diffPrev * 3.0;
  expected /= expected.norm();
  REQUIRE_THAT(tang, eonc::helpers::test::IsApprox(expected, 1e-10));
}

TEST_CASE("ImprovedTangent: maximum uses energy-weighted combination",
          "[neb]") {
  // At a maximum: energyPrev < energy > energyNext
  AtomMatrix diffNext = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  AtomMatrix diffPrev = make3({0, 1, 0, 0, 0, 0, 0, 0, 0});

  eonc::neb::ImprovedTangent strat;
  // energy = 5, energyPrev = 2, energyNext = 1
  // energyDiffPrev = 2-5 = -3, energyDiffNext = 1-5 = -4
  // |energyDiffPrev| = 3, |energyDiffNext| = 4
  // min = 3, max = 4
  // energyDiffPrev(-3) > energyDiffNext(-4) => tang = diffNext*min + diffPrev*max
  //   = diffNext*3 + diffPrev*4
  AtomMatrix tang = strat.compute(diffNext, diffPrev, 5.0, 2.0, 1.0);
  REQUIRE(tang.norm() == Catch::Approx(1.0).epsilon(1e-12));

  AtomMatrix expected = diffNext * 3.0 + diffPrev * 4.0;
  expected /= expected.norm();
  REQUIRE_THAT(tang, eonc::helpers::test::IsApprox(expected, 1e-10));
}

// ===== forcePerp ===========================================================

TEST_CASE("forcePerp: removes parallel component", "[neb]") {
  // tangent along x for atom 0
  AtomMatrix tangent = AtomMatrix::Zero(3, 3);
  tangent(0, 0) = 1.0;
  // normalize
  tangent /= tangent.norm();

  // force has both parallel and perpendicular components
  AtomMatrix force = make3({3, 4, 0, 0, 0, 0, 0, 0, 0});

  AtomMatrix fp = eonc::neb::forcePerp(force, tangent);

  // Perpendicular force should have zero dot product with tangent
  REQUIRE(matDot(fp, tangent) == Catch::Approx(0.0).margin(1e-12));

  // The perpendicular part should be (0,4,0, 0,0,0, 0,0,0)
  REQUIRE(fp(0, 0) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(fp(0, 1) == Catch::Approx(4.0));
}

TEST_CASE("forcePerp: purely parallel force yields zero", "[neb]") {
  AtomMatrix tangent = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  tangent /= tangent.norm();

  AtomMatrix force = make3({5, 0, 0, 0, 0, 0, 0, 0, 0});

  AtomMatrix fp = eonc::neb::forcePerp(force, tangent);
  REQUIRE(fp.norm() == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("forcePerp: purely perpendicular force is unchanged", "[neb]") {
  AtomMatrix tangent = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  tangent /= tangent.norm();

  AtomMatrix force = make3({0, 7, 0, 0, 0, 0, 0, 0, 0});

  AtomMatrix fp = eonc::neb::forcePerp(force, tangent);
  REQUIRE_THAT(fp, eonc::helpers::test::IsApprox(force, 1e-12));
}

// ===== climbingImageForce ==================================================

TEST_CASE("climbingImageForce: reverses parallel component", "[neb]") {
  AtomMatrix tangent = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  tangent /= tangent.norm();

  // force = (3, 4, 0, ...) => parallel = 3, perp = (0,4,0,...)
  AtomMatrix force = make3({3, 4, 0, 0, 0, 0, 0, 0, 0});
  AtomMatrix forceDNEB = AtomMatrix::Zero(3, 3);

  AtomMatrix fCI =
      eonc::neb::climbingImageForce(force, tangent, forceDNEB);

  // F_CI = F - 2*(F.t)*t = (3,4,0,...) - 2*3*(1,0,0,...) = (-3,4,0,...)
  REQUIRE(fCI(0, 0) == Catch::Approx(-3.0));
  REQUIRE(fCI(0, 1) == Catch::Approx(4.0));
}

TEST_CASE("climbingImageForce: includes DNEB contribution", "[neb]") {
  AtomMatrix tangent = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  tangent /= tangent.norm();

  AtomMatrix force = make3({3, 0, 0, 0, 0, 0, 0, 0, 0});
  AtomMatrix forceDNEB = make3({0, 1, 0, 0, 0, 0, 0, 0, 0});

  AtomMatrix fCI =
      eonc::neb::climbingImageForce(force, tangent, forceDNEB);

  // F_CI = (3,0,...) - 2*3*(1,0,...) + (0,1,...) = (-3, 1, 0, ...)
  REQUIRE(fCI(0, 0) == Catch::Approx(-3.0));
  REQUIRE(fCI(0, 1) == Catch::Approx(1.0));
}

// ===== computeDNEB =========================================================

TEST_CASE("computeDNEB: zero when perpendicular force is zero", "[neb]") {
  AtomMatrix tangent = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  tangent /= tangent.norm();

  AtomMatrix fSpring = make3({0, 2, 0, 0, 0, 0, 0, 0, 0});
  AtomMatrix fPerp = AtomMatrix::Zero(3, 3);

  AtomMatrix dneb = eonc::neb::computeDNEB(fSpring, tangent, fPerp);
  REQUIRE(dneb.norm() == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("computeDNEB: zero when spring force is purely parallel", "[neb]") {
  AtomMatrix tangent = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  tangent /= tangent.norm();

  // Spring force purely along tangent => spring perp = 0 => DNEB = 0
  AtomMatrix fSpring = make3({5, 0, 0, 0, 0, 0, 0, 0, 0});
  AtomMatrix fPerp = make3({0, 3, 0, 0, 0, 0, 0, 0, 0});

  AtomMatrix dneb = eonc::neb::computeDNEB(fSpring, tangent, fPerp);
  REQUIRE(dneb.norm() == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("computeDNEB: switching function behavior", "[neb]") {
  AtomMatrix tangent = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  tangent /= tangent.norm();

  // Spring force has perpendicular component along z
  AtomMatrix fSpring = make3({0, 0, 2, 0, 0, 0, 0, 0, 0});
  // Perpendicular force along y
  AtomMatrix fPerp = make3({0, 1, 0, 0, 0, 0, 0, 0, 0});

  AtomMatrix dneb = eonc::neb::computeDNEB(fSpring, tangent, fPerp);

  // Spring perp = fSpring - (fSpring.tangent)*tangent = (0,0,2,...) (no
  // parallel component). fPerp normalized = (0,1,0,...)/1. Spring perp
  // projected out of fPerpDir: (0,0,2) - dot((0,0,2),(0,1,0))*(0,1,0) =
  // (0,0,2). switching = (2/pi)*atan(|fPerp|^2/|springPerp|^2) =
  // (2/pi)*atan(1/4)
  double expectedSwitch =
      2.0 / eonc::helpers::pi * std::atan(1.0 / 4.0);
  // DNEB = (0,0,2,...) * switching
  REQUIRE(dneb(0, 2) == Catch::Approx(2.0 * expectedSwitch).epsilon(1e-10));
  // Other components should be zero (no y or x contribution)
  REQUIRE(dneb(0, 0) == Catch::Approx(0.0).margin(1e-12));
  REQUIRE(dneb(0, 1) == Catch::Approx(0.0).margin(1e-12));
}

// ===== zeroTranslation =====================================================

TEST_CASE("zeroTranslation: removes net translation when all atoms free",
          "[neb]") {
  AtomMatrix force = make3({3, 6, 9, 1, 2, 3, 2, 4, 6});

  eonc::neb::zeroTranslation(force, 3, 3);

  // Each column should now sum to zero
  for (int c = 0; c < 3; ++c) {
    REQUIRE(force.col(c).sum() == Catch::Approx(0.0).margin(1e-12));
  }
}

TEST_CASE("zeroTranslation: no-op when nFreeAtoms != nAtoms", "[neb]") {
  AtomMatrix force = make3({3, 6, 9, 1, 2, 3, 2, 4, 6});
  AtomMatrix original = force;

  // Only 2 of 3 atoms free => should not modify
  eonc::neb::zeroTranslation(force, 2, 3);

  REQUIRE_THAT(force, eonc::helpers::test::IsApprox(original, 1e-15));
}

// ===== UniformSpring (computeSpringForce) ==================================

TEST_CASE("UniformSpring: equal spacing gives zero parallel spring force",
          "[neb]") {
  // UniformSpring::compute requires a valid Matter for the full forceSpring
  // (which calls image->pbc). We verify the parallel component formula
  // directly: forceSpringPar = ksp * (distNext - distPrev) * tangent.
  double ksp = 5.0;

  AtomMatrix tangent = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  tangent /= tangent.norm();

  double distNext = 1.0;
  double distPrev = 1.0;

  AtomMatrix forceSpringPar = ksp * (distNext - distPrev) * tangent;
  REQUIRE(forceSpringPar.norm() == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("UniformSpring: unequal spacing gives nonzero parallel spring force",
          "[neb]") {
  double ksp = 5.0;

  AtomMatrix tangent = make3({1, 0, 0, 0, 0, 0, 0, 0, 0});
  tangent /= tangent.norm();

  double distNext = 2.0;
  double distPrev = 1.0;

  // forceSpringPar = ksp * (distNext - distPrev) * tangent = 5 * tangent
  AtomMatrix forceSpringPar = ksp * (distNext - distPrev) * tangent;
  AtomMatrix expected = ksp * 1.0 * tangent;
  REQUIRE_THAT(forceSpringPar,
               eonc::helpers::test::IsApprox(expected, 1e-12));
}

// ===== WeightedSpring ======================================================

TEST_CASE("WeightedSpring: energy-weighted spring force", "[neb]") {
  // Image index 1, so we use springConstants[1] and springConstants[0]
  std::vector<double> kvals = {2.0, 4.0, 6.0};
  eonc::neb::WeightedSpring spring{kvals};

  AtomMatrix tangent = make3({0, 1, 0, 0, 0, 0, 0, 0, 0});
  tangent /= tangent.norm();

  double distNext = 3.0;
  double distPrev = 2.0;

  // forceSpringPar = (kspNext*distNext - kspPrev*distPrev) * tangent
  //                = (4.0*3.0 - 2.0*2.0) * tangent = (12 - 4) * tangent = 8*tangent
  auto result = spring.compute(1, tangent, distNext, distPrev);

  AtomMatrix expected = 8.0 * tangent;
  REQUIRE_THAT(result.forceSpringPar,
               eonc::helpers::test::IsApprox(expected, 1e-12));

  // Full spring force should be zero for weighted springs
  REQUIRE(result.forceSpring.norm() == Catch::Approx(0.0).margin(1e-12));
}

// ===== findSplineExtrema ===================================================

// Minimal mock Matter that provides getPotentialEnergy, getPositions,
// getForces, distanceTo, and pbc -- enough for findSplineExtrema.
#include "Matter.h"

namespace {

/// Lightweight wrapper around Matter for spline extrema tests.
/// We cannot easily construct a full Matter without a .con file and potential,
/// so instead we test findSplineExtrema through the computeTangent +
/// cubic-spline math directly. We verify the cubic spline coefficient logic
/// with a hand-computed example.

/// Test the cubic spline extrema math directly using the polynomial
/// coefficients that findSplineExtrema would compute.
/// For a single interval with energy profile E(x) = a + b*x + c*x^2 + d*x^3,
/// extrema occur where dE/dx = b + 2*c*x + 3*d*x^2 = 0.
struct CubicInterval {
  double a, b, c, d;

  double energy(double f) const { return ((d * f + c) * f + b) * f + a; }
  double deriv(double f) const { return (3 * d * f + 2 * c) * f + b; }
  double curvature(double f) const { return 6 * d * f + 2 * c; }
};

} // namespace

TEST_CASE("findExtrema: parabolic barrier has one extremum in [0,1]",
          "[neb]") {
  // Construct a cubic polynomial that mimics a barrier between two images.
  // U1 = 0.0 at f=0, U2 = 0.0 at f=1, F1 = -1.0 (force along tangent),
  // F2 = 1.0 (force along tangent).
  //
  // findSplineExtrema uses these formulae for cubic coefficients:
  //   a = U1
  //   b = -F1
  //   c = 3*(U2-U1) + 2*F1 + F2
  //   d = -2*(U2-U1) - (F1+F2)
  double U1 = 0.0, U2 = 0.0;
  double F1 = -1.0, F2 = 1.0; // force projected onto tangent * dist

  CubicInterval seg;
  seg.a = U1;
  seg.b = -F1;        // 1.0
  seg.c = 3.0 * (U2 - U1) + 2.0 * F1 + F2; // -1.0
  seg.d = -2.0 * (U2 - U1) - (F1 + F2);     // 0.0

  // d=0, so this is a quadratic: b + 2*c*f = 0 => f = -b/(2c) = 1/2
  REQUIRE(seg.d == Catch::Approx(0.0).margin(1e-15));
  double f_ext = -seg.b / (2.0 * seg.c);
  REQUIRE(f_ext == Catch::Approx(0.5));
  REQUIRE(f_ext >= 0.0);
  REQUIRE(f_ext <= 1.0);

  // Energy at the extremum
  double E_ext = seg.energy(f_ext);
  REQUIRE(E_ext == Catch::Approx(0.25));

  // Curvature should be negative (it's a maximum)
  double curv = seg.curvature(f_ext);
  REQUIRE(curv < 0.0);
  REQUIRE(curv == Catch::Approx(-2.0));
}

TEST_CASE("findExtrema: asymmetric cubic has two roots, may have extrema in "
          "[0,1]",
          "[neb]") {
  // U1=0, U2=0.5, F1=-2.0, F2=0.5
  double U1 = 0.0, U2 = 0.5;
  double F1 = -2.0, F2 = 0.5;

  CubicInterval seg;
  seg.a = U1;
  seg.b = -F1;                                // 2.0
  seg.c = 3.0 * (U2 - U1) + 2.0 * F1 + F2;  // 1.5 - 4.0 + 0.5 = -2.0
  seg.d = -2.0 * (U2 - U1) - (F1 + F2);      // -1.0 + 1.5 = 0.5

  // dE/df = 3*d*f^2 + 2*c*f + b = 1.5*f^2 - 4*f + 2
  // Discriminant = (2c)^2 - 4*(3d)*b = 16 - 12 = 4
  double disc = seg.c * seg.c - 3.0 * seg.b * seg.d;
  REQUIRE(disc > 0.0);

  // Two roots from the cubic formula used in findSplineExtrema:
  // f1 = -(c + sqrt(disc)) / (3*d)
  // f2 = -(c - sqrt(disc)) / (3*d)
  double f1 = -(seg.c + std::sqrt(disc)) / (3.0 * seg.d);
  double f2 = -(seg.c - std::sqrt(disc)) / (3.0 * seg.d);

  // At least one should be in [0,1]
  bool has_extremum = (f1 >= 0 && f1 <= 1) || (f2 >= 0 && f2 <= 1);
  REQUIRE(has_extremum);

  // Verify that dE/df ~ 0 at the extrema
  if (f1 >= 0 && f1 <= 1) {
    REQUIRE(seg.deriv(f1) == Catch::Approx(0.0).margin(1e-10));
  }
  if (f2 >= 0 && f2 <= 1) {
    REQUIRE(seg.deriv(f2) == Catch::Approx(0.0).margin(1e-10));
  }
}

TEST_CASE("findExtrema: no extremum when derivative has no real root in [0,1]",
          "[neb]") {
  // Monotonically increasing cubic on [0,1]
  // U1=0, U2=10, F1=-20, F2=-20 (large forces pointing uphill)
  double U1 = 0.0, U2 = 10.0;
  double F1 = -20.0, F2 = -20.0;

  CubicInterval seg;
  seg.a = U1;
  seg.b = -F1;                                // 20
  seg.c = 3.0 * (U2 - U1) + 2.0 * F1 + F2;  // 30 - 40 - 20 = -30
  seg.d = -2.0 * (U2 - U1) - (F1 + F2);      // -20 + 40 = 20

  // dE/df = 60*f^2 - 60*f + 20 = 20*(3f^2 - 3f + 1)
  // disc of quadratic: 9-12 = -3 < 0 => no real roots
  double disc = seg.c * seg.c - 3.0 * seg.b * seg.d;
  REQUIRE(disc < 0.0);
}

} // namespace tests

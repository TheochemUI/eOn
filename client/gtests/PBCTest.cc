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
#include "catch2/catch_amalgamated.hpp"
#include "client/matter/PeriodicBoundaryConditions.hpp"

using namespace eonc;
using namespace eonc::mat;

Matrix3S createDefaultCell() {
  Matrix3S cell;
  cell << 10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0;
  return cell;
}

TEST_CASE("PBC_A1 wraps coordinates correctly", "[PBC_A1]") {
  Matrix3S cell = createDefaultCell();
  PBC_A1 pbc(cell);

  AtomMatrix diff(1, 3);
  diff << 15.0, -15.0, 25.0;

  AtomMatrix result = pbc(diff);

  REQUIRE(result(0, 0) == Catch::Approx(5.0));
  REQUIRE(result(0, 1) == Catch::Approx(-5.0));
  REQUIRE(result(0, 2) == Catch::Approx(5.0));
}

// TEST_CASE("PBC_B1 wraps coordinates correctly", "[PBC_B1]") {
//     Matrix3S cell = createDefaultCell();
//     PBC_B1 pbc(cell);

//     AtomMatrix diff(1, 3);
//     diff << 15.0, -15.0, 25.0;

//     AtomMatrix result = pbc(diff);

//     REQUIRE(result(0, 0) == Catch::Approx(5.0));
//     REQUIRE(result(0, 1) == Catch::Approx(-5.0));
//     REQUIRE(result(0, 2) == Catch::Approx(5.0));
// }

// TEST_CASE("PBC_B2 wraps coordinates correctly", "[PBC_B2]") {
//     Matrix3S cell = createDefaultCell();
//     PBC_B2 pbc(cell);

//     AtomMatrix diff(1, 3);
//     diff << 15.0, -15.0, 25.0;

//     AtomMatrix result = pbc(diff);

//     REQUIRE(result(0, 0) == Catch::Approx(5.0));
//     REQUIRE(result(0, 1) == Catch::Approx(-5.0));
//     REQUIRE(result(0, 2) == Catch::Approx(5.0));
// }

// TEST_CASE("PBC_B3 wraps coordinates correctly", "[PBC_B3]") {
//     Matrix3S cell = createDefaultCell();
//     PBC_B3 pbc(cell);

//     AtomMatrix diff(1, 3);
//     diff << 15.0, -15.0, 25.0;

//     AtomMatrix result = pbc(diff);

//     REQUIRE(result(0, 0) == Catch::Approx(5.0));
//     REQUIRE(result(0, 1) == Catch::Approx(-5.0));
//     REQUIRE(result(0, 2) == Catch::Approx(5.0));
// }

// TEST_CASE("PBC_B4 wraps coordinates correctly", "[PBC_B4]") {
//     Matrix3S cell = createDefaultCell();
//     PBC_B4 pbc(cell);

//     AtomMatrix diff(1, 3);
//     diff << 15.0, -15.0, 25.0;

//     AtomMatrix result = pbc(diff);

//     REQUIRE(result(0, 0) == Catch::Approx(5.0));
//     REQUIRE(result(0, 1) == Catch::Approx(-5.0));
//     REQUIRE(result(0, 2) == Catch::Approx(5.0));
// }

// TEST_CASE("PBC_C1 wraps coordinates correctly", "[PBC_C1]") {
//     Matrix3S cell = createDefaultCell();
//     PBC_C1 pbc(cell);

//     AtomMatrix diff(1, 3);
//     diff << 15.0, -15.0, 25.0;

//     AtomMatrix result = pbc(diff);

//     REQUIRE(result(0, 0) == Catch::Approx(5.0));
//     REQUIRE(result(0, 1) == Catch::Approx(-5.0));
//     REQUIRE(result(0, 2) == Catch::Approx(5.0));
// }

// TEST_CASE("PBC_C5 wraps coordinates correctly", "[PBC_C5]") {
//     Matrix3S cell = createDefaultCell();
//     PBC_C5 pbc(cell);

//     AtomMatrix diff(1, 3);
//     diff << 15.0, -15.0, 25.0;

//     AtomMatrix result = pbc(diff);

//     REQUIRE(result(0, 0) == Catch::Approx(5.0));
//     REQUIRE(result(0, 1) == Catch::Approx(-5.0));
//     REQUIRE(result(0, 2) == Catch::Approx(5.0));
// }

// TEST_CASE("PBC_C6 wraps coordinates correctly", "[PBC_C6]") {
//     Matrix3S cell = createDefaultCell();
//     PBC_C6 pbc(cell);

//     AtomMatrix diff(1, 3);
//     diff << 15.0, -15.0, 25.0;

//     AtomMatrix result = pbc(diff);

//     REQUIRE(result(0, 0) == Catch::Approx(5.0));
//     REQUIRE(result(0, 1) == Catch::Approx(-5.0));
//     REQUIRE(result(0, 2) == Catch::Approx(5.0));
// }

// TEST_CASE("PBC_ESVN wraps coordinates correctly", "[PBC_ESVN]") {
//     Matrix3S cell = createDefaultCell();
//     PBC_ESVN pbc(cell);

//     AtomMatrix diff(1, 3);
//     diff << 15.0, -15.0, 25.0;

//     AtomMatrix result = pbc(diff);

//     REQUIRE(result(0, 0) == Catch::Approx(5.0));
//     REQUIRE(result(0, 1) == Catch::Approx(-5.0));
//     REQUIRE(result(0, 2) == Catch::Approx(5.0));
// }

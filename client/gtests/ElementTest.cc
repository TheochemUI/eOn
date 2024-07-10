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
#include "Element.hpp"
#include "catch2/catch_amalgamated.hpp"

using namespace eonc;

TEST_CASE("mass2atom function test", "[mass2atom]") {
  SECTION("Valid atomic masses") {
    REQUIRE(mass2atom(1.008) == "H");
    REQUIRE(mass2atom(4.0026) == "He");
    REQUIRE(mass2atom(12.011) == "C");
    REQUIRE(mass2atom(238.03) == "U");
  }

  SECTION("Close to valid atomic masses") {
    REQUIRE(mass2atom(1.3) == "H");
    REQUIRE(mass2atom(4.2) == "He");
    REQUIRE(mass2atom(12.0) == "C");
    REQUIRE(mass2atom(238.0) == "U");
  }

  SECTION("Invalid atomic masses") {
    REQUIRE(mass2atom(0.0) == "Unknown");
    REQUIRE(mass2atom(1000.0) == "Unknown");
  }
}

TEST_CASE("symbol2atomicNumber function test", "[symbol2atomicNumber]") {
  SECTION("Valid symbols") {
    REQUIRE(symbol2atomicNumber("H") == 1);
    REQUIRE(symbol2atomicNumber("He") == 2);
    REQUIRE(symbol2atomicNumber("C") == 6);
    REQUIRE(symbol2atomicNumber("U") == 92);
  }

  SECTION("Invalid symbols") {
    REQUIRE(symbol2atomicNumber("X") == static_cast<size_t>(-1));
    REQUIRE(symbol2atomicNumber("Hydrogen") == static_cast<size_t>(-1));
    REQUIRE(symbol2atomicNumber("") == static_cast<size_t>(-1));
  }
}

TEST_CASE("atomicNumber2symbol function test", "[atomicNumber2symbol]") {
  SECTION("Valid atomic numbers") {
    REQUIRE(atomicNumber2symbol(1) == "H");
    REQUIRE(atomicNumber2symbol(2) == "He");
    REQUIRE(atomicNumber2symbol(6) == "C");
    REQUIRE(atomicNumber2symbol(92) == "U");
  }

  SECTION("Invalid atomic numbers") {
    REQUIRE(atomicNumber2symbol(0) == "Unknown");
    REQUIRE(atomicNumber2symbol(93) == "Unknown");
    REQUIRE(atomicNumber2symbol(100) == "Unknown");
  }
}

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
#include "eon/OHTSTJob.h"

TEST_CASE("pmf scan steps s uniformly from reactant to product", "[oh_tst]") {
  const double guideLen = 10.0;
  const long nScan = 40;
  REQUIRE(eonc::pmfScanS(0, nScan, guideLen) == 0.0);
  REQUIRE(eonc::pmfScanS(nScan - 1, nScan, guideLen) ==
          Catch::Approx(guideLen));
  double previous = -1.0;
  for (long plane = 0; plane < nScan; ++plane) {
    const double s = eonc::pmfScanS(plane, nScan, guideLen);
    REQUIRE(s == Catch::Approx(static_cast<double>(plane) * guideLen /
                               static_cast<double>(nScan - 1)));
    REQUIRE(s > previous);
    previous = s;
  }
  // s_init = 0.05 sits past the second uniform station, so a scan that
  // starts there and then takes (plane+1)*L/(n-1) steps backward.
  const double sInit = 0.05 * guideLen;
  const double oldSecond = guideLen / static_cast<double>(nScan - 1);
  REQUIRE(sInit > oldSecond);
  REQUIRE(eonc::pmfScanS(0, nScan, guideLen) <
          eonc::pmfScanS(1, nScan, guideLen));
}

TEST_CASE("pmf scan with fewer than two planes still spans the guideline",
          "[oh_tst]") {
  REQUIRE(eonc::pmfScanS(0, 1, 4.0) == 0.0);
  REQUIRE(eonc::pmfScanS(1, 1, 4.0) == Catch::Approx(4.0));
}

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
#include "eon/potentials/Rgpot/GroupSchedule.h"

#include "catch2/catch_amalgamated.hpp"

#include <cstdint>
#include <vector>

using eonc::GroupSchedule;

namespace {
std::vector<long> loads(const std::vector<int> &route, int groups) {
  std::vector<long> n(static_cast<size_t>(groups), 0);
  for (int g : route)
    ++n[static_cast<size_t>(g)];
  return n;
}
} // namespace

TEST_CASE("A full band keeps every image on its group", "[RGPOT][schedule]") {
  GroupSchedule s(4);
  const std::vector<std::int64_t> images{0, 1, 2, 3, 4, 5, 6, 7};
  const auto first = s.assign(images);
  REQUIRE(first == std::vector<int>{0, 1, 2, 3, 0, 1, 2, 3});
  REQUIRE(s.assign(images) == first);
}

TEST_CASE("A ring of more beads than groups runs in even rounds",
          "[RGPOT][schedule]") {
  GroupSchedule s(3);
  std::vector<std::int64_t> beads(16);
  for (int j = 0; j < 16; ++j)
    beads[static_cast<size_t>(j)] = j;
  const auto route = s.assign(beads);
  // ceil(16 / 3) = 6: no group takes more than one round beyond the rest.
  const auto n = loads(route, 3);
  REQUIRE(*std::max_element(n.begin(), n.end()) == 6);
  for (int step = 0; step < 5; ++step)
    REQUIRE(s.assign(beads) == route);
}

TEST_CASE("A partial band that would stack on one group is spread",
          "[RGPOT][schedule]") {
  GroupSchedule s(2);
  const std::vector<std::int64_t> band{0, 1, 2, 3, 4, 5, 6};
  static_cast<void>(s.assign(band));
  // Images 0, 2, 4 and 6 live on group 0. A batch of those four alone
  // would put all four there; the cap is two.
  const std::vector<std::int64_t> dirty{0, 2, 4, 6};
  const auto route = s.assign(dirty);
  REQUIRE(loads(route, 2) == std::vector<long>{2, 2});
  REQUIRE(route[0] == 0);
  REQUIRE(route[1] == 0);
  // The moved images take the new group as their home: the same batch
  // again lands where it did.
  REQUIRE(s.assign(dirty) == route);
}

TEST_CASE("A single request follows the nearest recorded geometry",
          "[RGPOT][schedule]") {
  GroupSchedule s(3);
  const std::vector<std::int64_t> band{0, 1, 2};
  const auto route = s.assign(band);
  const double a[3] = {0.0, 0.0, 0.0};
  const double b[3] = {1.0, 0.0, 0.0};
  const double c[3] = {2.0, 0.0, 0.0};
  s.record(0, a, 3);
  s.record(1, b, 3);
  s.record(2, c, 3);
  const double near_c[3] = {1.9, 0.1, 0.0};
  REQUIRE(s.nearest(near_c, 3) == route[2]);
  const double near_a[3] = {-0.2, 0.0, 0.0};
  REQUIRE(s.nearest(near_a, 3) == route[0]);
  // Another atom count matches nothing recorded: group 0.
  const double other[6] = {};
  REQUIRE(s.nearest(other, 6) == 0);
}

TEST_CASE("Keys without an owner fall back to the batch position",
          "[RGPOT][schedule]") {
  GroupSchedule s(2);
  REQUIRE(s.defaultHome(-2, 0) == 0);
  REQUIRE(s.defaultHome(-3, 1) == 1);
  REQUIRE(s.defaultHome(5, 0) == 1);
  REQUIRE(s.assign(std::vector<std::int64_t>{-2, -3}) ==
          std::vector<int>{0, 1});
}

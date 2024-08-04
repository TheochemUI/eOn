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
#include "client/IOHelpers.hpp"
#include "catch2/catch_amalgamated.hpp"

using namespace eonc;
TEST_CASE("getUniqueElements - Basic functionality") {
  VectorType vec(10);
  vec << 1, 2, 3, 2, 1, 4, 5, 4, 6, 5;
  VectorType result = io::getUniqueElements(vec);

  VectorType expected(6);
  expected << 1, 2, 3, 4, 5, 6;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

TEST_CASE("getUniqueElements - All unique elements") {
  VectorType vec(5);
  vec << 1, 2, 3, 4, 5;
  VectorType result = io::getUniqueElements(vec);

  VectorType expected(5);
  expected << 1, 2, 3, 4, 5;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

TEST_CASE("getUniqueElements - All elements are the same") {
  VectorType vec(5);
  vec << 1, 1, 1, 1, 1;
  VectorType result = io::getUniqueElements(vec);

  VectorType expected(1);
  expected << 1;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

TEST_CASE("getUniqueElements - Empty vector") {
  VectorType vec(0);
  VectorType result = io::getUniqueElements(vec);

  VectorType expected(0);

  REQUIRE(result.size() == expected.size());
}

TEST_CASE("getUniqueElements - Large vector with mixed elements") {
  VectorType vec(20);
  vec << 5, 3, 5, 7, 2, 2, 9, 7, 1, 3, 4, 5, 6, 9, 8, 0, 7, 6, 4, 3;
  VectorType result = io::getUniqueElements(vec);

  VectorType expected(10);
  expected << 5, 3, 7, 2, 9, 1, 4, 6, 8, 0;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

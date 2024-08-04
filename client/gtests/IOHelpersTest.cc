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
TEST_CASE("getUniqueValues - Basic functionality") {
  VectorType vec(10);
  vec << 1, 2, 3, 2, 1, 4, 5, 4, 6, 5;
  VectorType result = io::getUniqueValues(vec);

  VectorType expected(6);
  expected << 1, 2, 3, 4, 5, 6;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

TEST_CASE("getUniqueValues - All unique elements") {
  VectorType vec(5);
  vec << 1, 2, 3, 4, 5;
  VectorType result = io::getUniqueValues(vec);

  VectorType expected(5);
  expected << 1, 2, 3, 4, 5;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

TEST_CASE("getUniqueValues - All elements are the same") {
  VectorType vec(5);
  vec << 1, 1, 1, 1, 1;
  VectorType result = io::getUniqueValues(vec);

  VectorType expected(1);
  expected << 1;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

TEST_CASE("getUniqueValues - Empty vector") {
  VectorType vec(0);
  VectorType result = io::getUniqueValues(vec);

  VectorType expected(0);

  REQUIRE(result.size() == expected.size());
}

TEST_CASE("getUniqueValues - Large vector with mixed elements") {
  VectorType vec(20);
  vec << 5, 3, 5, 7, 2, 2, 9, 7, 1, 3, 4, 5, 6, 9, 8, 0, 7, 6, 4, 3;
  VectorType result = io::getUniqueValues(vec);

  VectorType expected(10);
  expected << 5, 3, 7, 2, 9, 1, 4, 6, 8, 0;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

TEST_CASE("getUniqueCounts - Basic functionality") {
  Vector<size_t> vec(10);
  vec << 1, 2, 3, 2, 1, 4, 5, 4, 6, 5;
  Vector<size_t> result = io::getUniqueCounts(vec);

  Vector<size_t> expected(6);
  expected << 2, 2, 1, 2, 2, 1;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

TEST_CASE("getUniqueCounts - All unique elements") {
  Vector<size_t> vec(5);
  vec << 1, 2, 3, 4, 5;
  Vector<size_t> result = io::getUniqueCounts(vec);

  Vector<size_t> expected(5);
  expected << 1, 1, 1, 1, 1;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

TEST_CASE("getUniqueCounts - All elements are the same") {
  Vector<size_t> vec(5);
  vec << 1, 1, 1, 1, 1;
  Vector<size_t> result = io::getUniqueCounts(vec);

  Vector<size_t> expected(1);
  expected << 5;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

TEST_CASE("getUniqueCounts - Empty vector") {
  Vector<size_t> vec(0);
  Vector<size_t> result = io::getUniqueCounts(vec);

  Vector<size_t> expected(0);

  REQUIRE(result.size() == expected.size());
}

TEST_CASE("getUniqueCounts - Large vector with mixed elements") {
  Vector<size_t> vec(20);
  vec << 5, 3, 5, 7, 2, 2, 9, 7, 1, 3, 4, 5, 6, 9, 8, 0, 7, 6, 4, 3;
  Vector<size_t> result = io::getUniqueCounts(vec);

  Vector<size_t> expected(10);
  expected << 3, 3, 3, 2, 2, 1, 2, 2, 1, 1;

  REQUIRE(result.size() == expected.size());
  for (int i = 0; i < result.size(); ++i) {
    REQUIRE(result(i) == Catch::Approx(expected(i)).epsilon(1e-6));
  }
}

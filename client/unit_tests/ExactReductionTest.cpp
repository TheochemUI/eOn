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
#include "eon/ExactReduction.h"
#include "catch2/catch_amalgamated.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Case {
  bool isDot = false;
  std::vector<double> x;
  std::vector<double> y;
  double expected = 0.0;
};

std::vector<double> readRow(std::istream &in, std::size_t n) {
  std::string line;
  std::getline(in, line);
  std::istringstream row(line);
  std::vector<double> v;
  v.reserve(n);
  std::string tok;
  while (row >> tok) {
    v.push_back(std::strtod(tok.c_str(), nullptr));
  }
  REQUIRE(v.size() == n);
  return v;
}

// Written by client/validation/exact_reduction_cases.py: the exact rational
// value rounded once to nearest, ties to even.
std::vector<Case> loadCases() {
  std::ifstream in("cases.txt");
  REQUIRE(in.good());
  std::vector<Case> cases;
  std::string head;
  while (std::getline(in, head)) {
    if (head.empty()) {
      continue;
    }
    std::istringstream h(head);
    std::string kind;
    std::size_t n = 0;
    h >> kind >> n;
    Case c;
    c.isDot = kind == "dot";
    c.x = readRow(in, n);
    if (c.isDot) {
      c.y = readRow(in, n);
    }
    c.expected = readRow(in, 1)[0];
    cases.push_back(std::move(c));
  }
  return cases;
}

uint64_t bitsOf(double v) { return std::bit_cast<uint64_t>(v); }

double evaluate(const Case &c) {
  return c.isDot ? eonc::exact::dot(c.x.data(), c.y.data(), c.x.size())
                 : eonc::exact::sum(c.x.data(), c.x.size());
}

} // namespace

TEST_CASE("exact reductions round the exact value once to nearest even",
          "[exact]") {
  const auto cases = loadCases();
  REQUIRE(cases.size() > 40);
  for (std::size_t k = 0; k < cases.size(); ++k) {
    INFO("case " << k << " n=" << cases[k].x.size());
    const double got = evaluate(cases[k]);
    // Bitwise: an exact zero is +0 by contract, the reference agrees.
    CHECK(bitsOf(got) == bitsOf(cases[k].expected));
  }
}

TEST_CASE("exact reductions do not depend on the summation order", "[exact]") {
  const auto cases = loadCases();
  std::mt19937_64 gen(42);
  for (std::size_t k = 0; k < cases.size(); ++k) {
    Case c = cases[k];
    const double first = evaluate(c);
    std::vector<std::size_t> perm(c.x.size());
    std::iota(perm.begin(), perm.end(), 0);
    for (int trial = 0; trial < 16; ++trial) {
      std::shuffle(perm.begin(), perm.end(), gen);
      Case p = c;
      for (std::size_t i = 0; i < perm.size(); ++i) {
        p.x[i] = c.x[perm[i]];
        if (c.isDot) {
          p.y[i] = c.y[perm[i]];
        }
      }
      INFO("case " << k << " trial " << trial);
      CHECK(bitsOf(evaluate(p)) == bitsOf(first));
    }
  }
}

TEST_CASE("exact reductions pass non-finite inputs through as IEEE sums",
          "[exact]") {
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<double> a{1.0, inf, 2.0};
  const std::vector<double> b{1.0, 1.0, 1.0};
  CHECK(eonc::exact::dot(a.data(), b.data(), 3) == inf);
  CHECK(eonc::exact::sum(a.data(), 3) == inf);
  const std::vector<double> c{1.0, nan};
  CHECK(std::isnan(eonc::exact::sum(c.data(), 2)));
  const std::vector<double> d{inf, -inf};
  CHECK(std::isnan(eonc::exact::sum(d.data(), 2)));
  CHECK(eonc::exact::sum(nullptr, 0) == 0.0);
  CHECK(!std::signbit(eonc::exact::sum(nullptr, 0)));
}

TEST_CASE("exact norm and normalized agree with their definitions", "[exact]") {
  Eigen::VectorXd v(4);
  v << 3.0, 4.0, 12.0, 84.0;
  CHECK(eonc::exact::squaredNorm(v) == 7225.0);
  CHECK(eonc::exact::norm(v) == 85.0);
  const Eigen::VectorXd u = eonc::exact::normalized(v);
  CHECK(u[0] == 3.0 / 85.0);
  CHECK(eonc::exact::normalized(Eigen::VectorXd::Zero(3)).isZero(0.0));
  // A segment and a strided expression go through a contiguous copy.
  CHECK(eonc::exact::dot(v.head(2), v.tail(2)) == 3.0 * 12.0 + 4.0 * 84.0);
}

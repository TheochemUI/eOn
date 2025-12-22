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
#pragma once

#include "catch2/catch_amalgamated.hpp"

namespace helper_functions::test {

template <typename T>
class EigenMatcher : public Catch::Matchers::MatcherGenericBase {
  const T &m_expected;
  double m_threshold;

public:
  EigenMatcher(const T &expected, double threshold)
      : m_expected(expected),
        m_threshold(threshold) {}

  bool match(const T &actual) const {
    if (actual.rows() != m_expected.rows() ||
        actual.cols() != m_expected.cols()) {
      return false;
    }
    return actual.isApprox(m_expected, m_threshold);
  }

  std::string describe() const override {
    std::stringstream ss;
    ss << "\nis approximately equal to:\n"
       << m_expected << "\nwith threshold: " << m_threshold;
    return ss.str();
  }
};

// Catch2-style helper function to keep syntax clean
template <typename T>
auto IsApprox(const T &expected, double threshold = 1e-4) {
  return EigenMatcher<T>(expected, threshold);
}

} // namespace helper_functions::test

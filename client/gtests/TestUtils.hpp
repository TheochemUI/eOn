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

#include "EonLogger.h"
#include "catch2/catch_amalgamated.hpp"
#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/sinks/NullSink.h"

namespace eonc::helpers::test {

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

/// \brief RAII helper to initialize quill logging for tests.
///
/// Many eOn internals (CIniFile::GetValue, Matter, Potential) call
/// quill::Frontend::get_logger("combi"). If the backend isn't started
/// or the logger isn't registered, get_logger() returns nullptr → segfault.
///
/// This struct ensures:
/// 1. quill::Backend is started before tests run
/// 2. A "combi" logger exists (with NullSink to discard output)
/// 3. Backend is stopped cleanly on teardown
///
/// Usage: Just instantiate as a static global in your test file:
///   static helpers::test::QuillTestLogger _quill_setup;
struct QuillTestLogger {
  QuillTestLogger() {
    if (!eonc::log::get()) {
      quill::Backend::start();
      auto null_sink = quill::Frontend::create_or_get_sink<quill::NullSink>(
          "null_test_sink");
      quill::Frontend::create_or_get_logger(
          "combi", std::move(null_sink),
          quill::PatternFormatterOptions{"%(message)"},
          quill::ClockSourceType::System);
    }
  }
  ~QuillTestLogger() = default;
};

} // namespace eonc::helpers::test

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
#include "eon/EonLogger.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eonc_test_aliases.hpp"
#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/sinks/NullSink.h"
#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <set>
#include <thread>

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
/// or the logger isn't registered, get_logger() returns nullptr and the
/// next log crashes.
///
/// This struct ensures:
/// 1. quill::Backend is started before tests run
/// 2. A "combi" logger exists (with NullSink to discard output)
/// 3. Backend is stopped cleanly on teardown
///
/// Do not probe the logger with eonc::log::get(). That fallback creates
/// a console logger named "combi" when none exists, and this constructor
/// then leaves it in place. Dimer and Lanczos info lines go to that
/// stdout sink. Meson captures the pipe; on Windows the pipe buffer is
/// a few kilobytes and fwrite blocks the logging thread until the
/// runner reads. test_bench_parity fills the pipe and misses its 30s
/// limit. Register the null sink against Frontend directly.
///
/// Usage: Just instantiate as a static global in your test file:
///   static helpers::test::QuillTestLogger _quill_setup;
struct QuillTestLogger {
  QuillTestLogger() {
    if (quill::Frontend::get_logger("combi")) {
      return;
    }
    try {
      quill::Backend::start();
    } catch (...) {
    }
    try {
      auto null_sink = quill::Frontend::create_or_get_sink<quill::NullSink>(
          "null_test_sink");
      quill::Frontend::create_or_get_logger(
          "combi", std::move(null_sink),
          quill::PatternFormatterOptions{"%(message)"},
          quill::ClockSourceType::System);
    } catch (...) {
    }
  }
  ~QuillTestLogger() = default;
};

/// LJ forces through a potential that asks for one instance per image, as
/// rgpot's neighbour-table kernels do. Every instance records the threads
/// that called it into one shared set.
class PerImageLJ final : public eonc::Potential {
public:
  struct Threads {
    std::mutex mu;
    std::set<std::thread::id> ids;
  };

  explicit PerImageLJ(std::shared_ptr<Threads> threads,
                      std::function<void()> on_force = {})
      : Potential(PotType::LJ),
        threads_{std::move(threads)},
        on_force_{std::move(on_force)},
        lj_{makeLJ()} {}
  PerImageLJ(std::shared_ptr<Threads> threads, eonc::IPotRegistry &registry,
             std::function<void()> on_force = {})
      : Potential(PotType::LJ, registry),
        threads_{std::move(threads)},
        on_force_{std::move(on_force)},
        lj_{makeLJ()} {}

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return true;
  }
  [[nodiscard]] std::shared_ptr<eonc::Potential>
  clonePotential() const override {
    return std::make_shared<PerImageLJ>(threads_, on_force_);
  }

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    {
      std::lock_guard<std::mutex> lock(threads_->mu);
      threads_->ids.insert(std::this_thread::get_id());
    }
    if (on_force_) {
      on_force_();
    }
    lj_->force(nAtoms, positions, atomicNrs, forces, energy, variance, box);
  }

private:
  static std::unique_ptr<eonc::Potential> makeLJ() {
    Parameters params;
    return eonc::helpers::makePotential(PotType::LJ, params);
  }

  std::shared_ptr<Threads> threads_;
  std::function<void()> on_force_;
  std::unique_ptr<eonc::Potential> lj_;
};

/// Counts the kernel evaluations behind eOn's own forceCallCounter, so a
/// test can check that a path through a library still counts every one.
class CountingPotential final : public eonc::Potential {
public:
  explicit CountingPotential(std::shared_ptr<eonc::Potential> inner)
      : Potential(PotType::LJ),
        inner_{std::move(inner)} {}

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    ++kernelCalls;
    inner_->force(nAtoms, positions, atomicNrs, forces, energy, variance, box);
  }
  [[nodiscard]] double finiteCutoff() const noexcept override {
    return inner_->finiteCutoff();
  }

  std::atomic<size_t> kernelCalls{0};

private:
  std::shared_ptr<eonc::Potential> inner_;
};

} // namespace eonc::helpers::test

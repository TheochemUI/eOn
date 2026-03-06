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

#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/Logger.h"
#include "quill/sinks/FileSink.h"
#include <source_location>
#include <string_view>

namespace eonc::log {

/// \brief Get or create the default "combi" logger.
///
/// This is the primary logger used throughout eOn. Returns a raw pointer
/// (quill's API) that is valid for the lifetime of the program.
///
/// Usage:
///   auto* log = eonc::log::get();
///   LOG_INFO(log, "message");
[[nodiscard]] inline quill::Logger *get() noexcept {
  return quill::Frontend::get_logger("combi");
}

/// \brief Get or create a named logger with a file sink.
///
/// Creates a logger that writes to the specified file with the given pattern.
/// If the logger already exists, returns the existing instance.
///
/// \param name Logger name (must be unique)
/// \param filename Output file path
/// \param pattern Quill format pattern (default: "%(message)")
/// \return Raw pointer to the logger (valid for program lifetime)
[[nodiscard]] inline quill::Logger *
get_file(std::string_view name, std::string_view filename,
         std::string_view pattern = "%(message)") {
  quill::FileSinkConfig cfg;
  cfg.set_open_mode('w');

  auto file_sink = quill::Frontend::create_or_get_sink<quill::FileSink>(
      std::string(filename), cfg);

  return quill::Frontend::create_or_get_logger(
      std::string(name), std::move(file_sink),
      quill::PatternFormatterOptions{std::string(pattern)},
      quill::ClockSourceType::System);
}

/// \brief RAII helper for class-scoped logging.
///
/// Automatically retrieves the "combi" logger on construction. Designed for
/// member variables: just declare `eonc::log::Scoped m_log;` and use it.
///
/// Example:
/// \code
///   class MyClass {
///     eonc::log::Scoped m_log;  // No manual initialization needed!
///   public:
///     void work() {
///       LOG_INFO(m_log, "doing work");
///     }
///   };
/// \endcode
struct Scoped {
  quill::Logger *logger{nullptr};

  Scoped() noexcept
      : logger(get()) {}

  // Implicit conversion to Logger* for seamless use with LOG_* macros
  operator quill::Logger *() const noexcept { return logger; }

  // Support LOG macro pointer dereference
  quill::Logger *operator->() const noexcept { return logger; }

  // Explicit access if needed
  [[nodiscard]] quill::Logger *ptr() const noexcept { return logger; }
};

/// \brief RAII helper for file-based logging in a class.
///
/// Similar to Scoped, but writes to a dedicated file. Useful for optimizer
/// traces, detailed diagnostics, etc.
///
/// Example:
/// \code
///   class Optimizer {
///     eonc::log::FileScoped m_log{"opt_trace", "optimizer.log"};
///   public:
///     void step() {
///       LOG_INFO(m_log, "iteration {}", iter);
///     }
///   };
/// \endcode
struct FileScoped {
  quill::Logger *logger{nullptr};

  FileScoped(std::string_view name, std::string_view filename,
             std::string_view pattern = "%(message)") noexcept
      : logger(get_file(name, filename, pattern)) {}

  operator quill::Logger *() const noexcept { return logger; }
  [[nodiscard]] quill::Logger *ptr() const noexcept { return logger; }
};

} // namespace eonc::log

// Force MSVC to evaluate __VA_ARGS__ properly
#define EONC_EXPAND(x) x

/// \brief Convenience macro for one-shot logging without storing a logger.
///
/// Use when you only need to log once or twice in a function and don't want
/// to declare a logger variable.
///
/// Example:
///   EONC_LOG_INFO("quick message");
///   EONC_LOG_DEBUG("value: {}", x);
#define EONC_LOG_TRACE(...)                                                    \
  EONC_EXPAND(LOG_TRACE_L1(eonc::log::get(), __VA_ARGS__))
#define EONC_LOG_DEBUG(...)                                                    \
  EONC_EXPAND(LOG_DEBUG(eonc::log::get(), __VA_ARGS__))
#define EONC_LOG_INFO(...) EONC_EXPAND(LOG_INFO(eonc::log::get(), __VA_ARGS__))
#define EONC_LOG_WARNING(...)                                                  \
  EONC_EXPAND(LOG_WARNING(eonc::log::get(), __VA_ARGS__))
#define EONC_LOG_ERROR(...)                                                    \
  EONC_EXPAND(LOG_ERROR(eonc::log::get(), __VA_ARGS__))
#define EONC_LOG_CRITICAL(...)                                                 \
  EONC_EXPAND(LOG_CRITICAL(eonc::log::get(), __VA_ARGS__))

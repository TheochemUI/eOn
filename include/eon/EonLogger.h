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
#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include "quill/sinks/ConsoleSink.h"
#include "quill/sinks/FileSink.h"
#include <atomic>
#include <chrono>
#include <filesystem>
#include <mutex>
#include <source_location>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>

namespace eonc::log {

// Logger name constants for internal use
namespace detail {
inline constexpr std::string_view kDefaultLoggerName{"combi"};
inline constexpr std::string_view kTracebackLoggerName{"_traceback"};
} // namespace detail


/// How queued log records reach the sinks.
///
/// Thread runs Quill's own backend thread. Manual runs no thread at all:
/// the thread that logs polls the queue itself through poll(), flush()
/// and shutdown(), so a client can be strictly single threaded.
enum class BackendMode { Thread, Manual };

namespace detail {

inline std::atomic<quill::ManualBackendWorker *> g_manual_worker{nullptr};
inline std::atomic<std::thread::id> g_manual_owner{};
inline std::atomic<std::chrono::steady_clock::rep> g_last_poll{0};

/// Poll no more often than this in Manual mode. Matches the flush
/// interval of the Thread backend, so output is as timely.
inline constexpr std::chrono::milliseconds kManualPollInterval{100};

[[nodiscard]] inline quill::ManualBackendWorker *manual_worker_for_caller() {
  auto *worker = g_manual_worker.load(std::memory_order_acquire);
  if (worker && g_manual_owner.load() == std::this_thread::get_id()) {
    return worker;
  }
  return nullptr;
}

} // namespace detail

inline BackendMode start_client_backend(BackendMode mode) noexcept;

/// \brief True once the client queue is drained by its own thread.
[[nodiscard]] inline bool manual_backend() noexcept {
  return detail::g_manual_worker.load(std::memory_order_acquire) != nullptr;
}

/// \brief Write every queued record to the sinks and flush the sinks.
///
/// In Manual mode only the thread that started the backend can do this;
/// a call from another thread leaves the records queued for that thread.
/// In Thread mode the backend thread does the work and this returns once
/// it has.
inline void flush(quill::Logger *logger) noexcept {
  if (manual_backend()) {
    if (auto *worker = detail::manual_worker_for_caller()) {
      worker->poll();
      // An idle poll flushes the sinks (sink_min_flush_interval is zero).
      worker->poll_one();
    }
    return;
  }
  // The flush record waits for a backend. One that was never started would
  // block here for good, and one that was stopped cannot answer.
  if (!quill::Backend::is_running()) {
    start_client_backend(BackendMode::Thread);
    if (!quill::Backend::is_running()) {
      return;
    }
  }
  if (logger) {
    logger->flush_log();
  }
}

/// \brief Drain the queue now and then, from the polling thread.
///
/// Cheap when the backend runs its own thread or when called from a
/// thread that does not own the queue. Place it on a path the compute
/// thread revisits often.
inline void poll_if_due() noexcept {
  if (!manual_backend()) {
    return;
  }
  auto *worker = detail::manual_worker_for_caller();
  if (!worker) {
    return;
  }
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto last = detail::g_last_poll.load(std::memory_order_relaxed);
  if (now - last <
      std::chrono::duration_cast<std::chrono::steady_clock::duration>(
          detail::kManualPollInterval)
          .count()) {
    return;
  }
  detail::g_last_poll.store(now, std::memory_order_relaxed);
  worker->poll();
  worker->poll_one();
}

/// \brief Drain the queue and stop the backend. Safe to call twice.
inline void shutdown() noexcept {
  if (manual_backend()) {
    flush(nullptr);
    return;
  }
  // Records logged before any backend existed are written by this one.
  if (!quill::Backend::is_running()) {
    start_client_backend(BackendMode::Thread);
  }
  quill::Backend::stop();
}

/// \brief Create the combi and traceback loggers once, without a backend.
///
/// Records logged before start_client_backend() wait in the queue, so the
/// parameter file can be read first and still report its own errors.
/// Not called from process load. Help, version, and feature flags return
/// before this. Tests that already registered "combi" are left alone.
[[nodiscard]] inline quill::Logger *init_client_loggers() noexcept {
  static std::once_flag once;
  std::call_once(once, [] {
    if (quill::Frontend::get_logger(std::string(detail::kDefaultLoggerName))) {
      return;
    }
    try {
      auto console_sink =
          quill::Frontend::create_or_get_sink<quill::ConsoleSink>("console");
      auto file_sink = quill::Frontend::create_or_get_sink<quill::FileSink>(
          "client_quill.log",
          []() {
            quill::FileSinkConfig cfg;
            cfg.set_open_mode('w');
            return cfg;
          }(),
          quill::FileEventNotifier{});
      auto *logger = quill::Frontend::create_or_get_logger(
          "combi", {std::move(console_sink), std::move(file_sink)},
          quill::PatternFormatterOptions{"%(message)"},
          quill::ClockSourceType::System);
      logger->set_log_level(quill::LogLevel::TraceL3);
      auto trace_csink =
          quill::Frontend::create_or_get_sink<quill::ConsoleSink>(
              "trace_console");
      auto trace_fsink = quill::Frontend::create_or_get_sink<quill::FileSink>(
          "client_traceback.log",
          []() {
            quill::FileSinkConfig cfg;
            cfg.set_open_mode('w');
            return cfg;
          }(),
          quill::FileEventNotifier{});
      quill::Frontend::create_or_get_logger(
          "_traceback", {std::move(trace_csink), std::move(trace_fsink)},
          quill::PatternFormatterOptions{
              " [%(log_level)] [%(source_location)] [%(caller_function)] \n "
              "%(message)\n[end %(log_level)]"},
          quill::ClockSourceType::System);
      // Build the logger and thread context singletons before the backend
      // one. Statics are destroyed in reverse order, so the manual worker
      // then drains its queue at exit while both are still alive.
      quill::Frontend::preallocate();
    } catch (...) {
    }
  });
  return quill::Frontend::get_logger(std::string(detail::kDefaultLoggerName));
}

/// \brief Start draining the client queue, once. Later calls do nothing.
///
/// Manual installs an exit handler that drains the queue, since Quill's
/// own exit handler belongs to its backend thread. A backend that is
/// already running (a test fixture started one) is left in charge.
/// \return the mode now in effect.
inline BackendMode start_client_backend(BackendMode mode) noexcept {
  static std::once_flag once;
  static BackendMode effective{BackendMode::Thread};
  std::call_once(once, [mode] {
    if (quill::Backend::is_running() || manual_backend()) {
      return;
    }
    try {
      quill::BackendOptions backend_options;
      backend_options.transit_event_buffer_initial_capacity = 2048;
      backend_options.log_timestamp_ordering_grace_period =
          std::chrono::microseconds{0};
      backend_options.check_printable_char = {};
      if (mode == BackendMode::Manual) {
        backend_options.sink_min_flush_interval = std::chrono::milliseconds{0};
        auto *worker = quill::Backend::acquire_manual_backend_worker();
        worker->init(backend_options);
        detail::g_manual_owner.store(std::this_thread::get_id());
        detail::g_manual_worker.store(worker, std::memory_order_release);
        std::atexit([] { shutdown(); });
        effective = BackendMode::Manual;
        return;
      }
      backend_options.sleep_duration = std::chrono::microseconds{10};
      backend_options.sink_min_flush_interval = std::chrono::milliseconds{100};
      quill::Backend::start(backend_options);
    } catch (...) {
    }
  });
  return effective;
}

/// \brief Start the client backend thread and the combi / traceback
/// loggers once.
[[nodiscard]] inline quill::Logger *init_client() noexcept {
  auto *logger = init_client_loggers();
  start_client_backend(BackendMode::Thread);
  return logger;
}

/// \brief Get or create the default "combi" logger.
///
/// This is the primary logger used throughout eOn. Returns a raw pointer
/// (quill's API) that is valid for the lifetime of the program.
///
/// Usage:
///   auto* log = eonc::log::get();
///   QUILL_LOG_INFO(log, "message");
[[nodiscard]] inline quill::Logger *get() noexcept {
  static std::atomic<quill::Logger *> cached_logger{nullptr};
  if (quill::Logger *l = cached_logger.load(std::memory_order_relaxed)) {
    return l;
  }

  if (auto *logger = quill::Frontend::get_logger(
          std::string(detail::kDefaultLoggerName))) {
    cached_logger.store(logger, std::memory_order_relaxed);
    return logger;
  }

  // Fallback for Python bindings/tests that didn't run ClientEON's main()
  static quill::Logger *fallback = []() -> quill::Logger * {
    try {
      quill::Backend::start();
    } catch (...) {
    } // Ignore if backend already started

    try {
      auto console_sink =
          quill::Frontend::create_or_get_sink<quill::ConsoleSink>(
              "fallback_console");
      return quill::Frontend::create_or_get_logger(
          "combi", std::move(console_sink),
          quill::PatternFormatterOptions{"%(message)"},
          quill::ClockSourceType::System);
    } catch (...) {
      return nullptr;
    }
  }();
  if (fallback)
    cached_logger.store(fallback, std::memory_order_relaxed);
  return fallback;
}

/// \brief Get or create the "_traceback" logger for traceback logging.
///
/// This logger is used for detailed traceback and debugging information.
/// Returns a raw pointer (quill's API) valid for the lifetime of the program.
///
/// Usage:
///   auto* log = eonc::log::traceback();
///   QUILL_LOG_INFO(log, "traceback message");
[[nodiscard]] inline quill::Logger *traceback() noexcept {
  static std::atomic<quill::Logger *> cached{nullptr};
  if (quill::Logger *l = cached.load(std::memory_order_relaxed)) {
    return l;
  }
  if (auto *logger = quill::Frontend::get_logger(
          std::string(detail::kTracebackLoggerName))) {
    cached.store(logger, std::memory_order_relaxed);
    return logger;
  }
  return nullptr;
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
  // If the logger already exists, return it immediately.
  // Quill strictly checks configuration matches on recreation and throws if
  // PatternFormatterOptions differ, which crashes instances (like
  // Dimer/Optimizers) that are created multiple times in loops (e.g., NEB).
  if (auto *existing = quill::Frontend::get_logger(std::string(name))) {
    return existing;
  }
  static std::mutex file_logger_mutex;
  std::lock_guard<std::mutex> lock(file_logger_mutex);
  if (auto *existing = quill::Frontend::get_logger(std::string(name))) {
    return existing;
  }

  try {
    // Durable absolute path: optimizers open "_lbfgs.log" etc. under CWD.
    // Job integration fixtures chdir into /tmp/eon_test_job_* then delete that
    // tree when the case ends, while Quill's process-lifetime FileSink still
    // flushes — "fopen failed after 3 attempts". Keep logs under a stable
    // temp subdir instead of the disposable workdir.
    namespace fs = std::filesystem;
    fs::path path{std::string(filename)};
    if (!path.is_absolute()) {
      static const fs::path log_dir = [] {
        fs::path d = fs::temp_directory_path() / "eon_logs";
        std::error_code ec;
        fs::create_directories(d, ec);
        return d;
      }();
      path = log_dir / path.filename();
    }
    quill::FileSinkConfig cfg;
    cfg.set_open_mode('w');
    auto file_sink = quill::Frontend::create_or_get_sink<quill::FileSink>(
        path.string(), cfg);
    return quill::Frontend::create_or_get_logger(
        std::string(name), std::move(file_sink),
        quill::PatternFormatterOptions{std::string(pattern)},
        quill::ClockSourceType::System);
  } catch (...) {
    return get(); // Fallback to default combi logger if creation fails
  }
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
///       QUILL_LOG_INFO(m_log, "doing work");
///     }
///   };
/// \endcode
struct Scoped {
  [[no_unique_address]] quill::Logger *logger{nullptr};

  Scoped() noexcept : logger(get()) {}

  // Implicit conversion to Logger* for seamless use with QUILL_LOG_* macros
  operator quill::Logger *() const noexcept { return logger; }

  // Support QUILL_LOG_ macro pointer dereference
  quill::Logger *operator->() const noexcept { return logger; }

  // Allow dynamic re-assignment (e.g. switching to _traceback logger)
  Scoped &operator=(quill::Logger *l) noexcept {
    logger = l;
    return *this;
  }

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
///       QUILL_LOG_INFO(m_log, "iteration {}", iter);
///     }
///   };
/// \endcode
struct FileScoped {
  [[no_unique_address]] quill::Logger *logger{nullptr};

  FileScoped(std::string_view name, std::string_view filename,
             std::string_view pattern = "%(message)") noexcept
      : logger(get_file(name, filename, pattern)) {}

  operator quill::Logger *() const noexcept { return logger; }
  // Support QUILL_LOG_ macro pointer dereference
  quill::Logger *operator->() const noexcept { return logger; }
  FileScoped &operator=(quill::Logger *l) noexcept {
    logger = l;
    return *this;
  }
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
  do {                                                                         \
    if (auto *l = eonc::log::get()) {                                          \
      EONC_EXPAND(QUILL_LOG_TRACE_L1(l, __VA_ARGS__));                         \
    }                                                                          \
  } while (0)
#define EONC_LOG_DEBUG(...)                                                    \
  do {                                                                         \
    if (auto *l = eonc::log::get()) {                                          \
      EONC_EXPAND(QUILL_LOG_DEBUG(l, __VA_ARGS__));                            \
    }                                                                          \
  } while (0)
#define EONC_LOG_INFO(...)                                                     \
  do {                                                                         \
    if (auto *l = eonc::log::get()) {                                          \
      EONC_EXPAND(QUILL_LOG_INFO(l, __VA_ARGS__));                             \
    }                                                                          \
  } while (0)
#define EONC_LOG_WARNING(...)                                                  \
  do {                                                                         \
    if (auto *l = eonc::log::get()) {                                          \
      EONC_EXPAND(QUILL_LOG_WARNING(l, __VA_ARGS__));                          \
    }                                                                          \
  } while (0)
#define EONC_LOG_ERROR(...)                                                    \
  do {                                                                         \
    if (auto *l = eonc::log::get()) {                                          \
      EONC_EXPAND(QUILL_LOG_ERROR(l, __VA_ARGS__));                            \
    }                                                                          \
  } while (0)
#define EONC_LOG_CRITICAL(...)                                                 \
  do {                                                                         \
    if (auto *l = eonc::log::get()) {                                          \
      EONC_EXPAND(QUILL_LOG_CRITICAL(l, __VA_ARGS__));                         \
    }                                                                          \
  } while (0)

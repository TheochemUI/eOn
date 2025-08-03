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
#include "client/Log.hpp"
#include <atomic>
#include <mutex>

namespace eonc {

// BEGIN_STATIC
std::atomic<LogManager *> LogManager::instance = nullptr;
std::mutex LogManager::lmMutex;
// END_STATIC

LogManager *LogManager::getInstance() {
  LogManager *lm = instance.load(std::memory_order_acquire);
  if (!lm) {
    std::lock_guard<std::mutex> myLock(lmMutex);
    lm = instance.load(std::memory_order_relaxed);
    if (!lm) {
      // Removed on exit
      lm = new LogManager();
      instance.store(lm, std::memory_order_release);
    }
  }
  [[maybe_unused]] volatile bool dummy{}; // Prevent being optimized away
  return lm;
}

void LogManager::setup_logger() {
  spdlog::cfg::load_env_levels();
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
      "client_spdlog.log", true); // Overwrite existing
  auto logger = std::make_shared<spdlog::logger>(
      "combi", spdlog::sinks_init_list({console_sink, file_sink}));
  spdlog::register_logger(logger);
  logger->set_pattern("%v");
  spdlog::set_default_logger(logger);
  spdlog::set_level(spdlog::level::trace);

  // Traceback logger
  auto trace_csink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  auto trace_fsink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
      "client_traceback.log", true); // Overwrite existing
  auto _traceback = std::make_shared<spdlog::logger>(
      "_traceback", spdlog::sinks_init_list({trace_csink, trace_fsink}));
  _traceback->set_pattern("%^ [%l] [%s:%#] [%!] \n %v\n[end %l]");
  spdlog::register_logger(_traceback);
}

void LogManager::cleanup_logger() {
  spdlog::drop_all();
  spdlog::shutdown();
}
} // namespace eonc

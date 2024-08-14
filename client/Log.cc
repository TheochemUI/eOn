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

namespace eonc {

void LogManager::setup_logger() {
  spdlog::cfg::load_env_levels();
  // Sinks
  spdlog::flush_every(std::chrono::seconds(3));
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
      "client_spdlog.log", true); // Overwrite existing
  auto logger = std::make_shared<spdlog::logger>(
      "combi", spdlog::sinks_init_list({console_sink, file_sink}));
  spdlog::register_logger(logger);
  logger->set_pattern("%v");
  spdlog::set_default_logger(logger);
  // Traceback logger
  spdlog::set_level(spdlog::level::trace);
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

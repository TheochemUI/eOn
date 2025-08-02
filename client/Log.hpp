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
#include <spdlog/cfg/env.h> // support for loading levels from the environment variable
#include <spdlog/fmt/ostr.h> // support for user defined types
#include <spdlog/logger.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace eonc {

// As per
// https://www.modernescpp.com/index.php/thread-safe-initialization-of-a-singleton/
class LogManager {
public:
  static LogManager *getInstance();

private:
  void setup_logger();
  void cleanup_logger();
  LogManager() { setup_logger(); }
  ~LogManager() { cleanup_logger(); }
  LogManager(const LogManager &) = delete;
  LogManager &operator=(const LogManager &) = delete;

  static std::atomic<LogManager *> instance;
  static std::mutex lmMutex;
};

} // namespace eonc

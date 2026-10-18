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

// The client logger as eonclient starts it: init_client, not the null
// test logger, so the backend runs with the client's options.

#include "catch2/catch_amalgamated.hpp"
#include "eon/EonLogger.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

namespace {

namespace fs = std::filesystem;

#if defined(__linux__)
// Voluntary context switches of the backend thread, read from /proc.
long backendSwitches() {
  for (const auto &task : fs::directory_iterator("/proc/self/task")) {
    std::ifstream comm(task.path() / "comm");
    std::string name;
    std::getline(comm, name);
    if (name.rfind("Quill", 0) != 0) {
      continue;
    }
    std::ifstream status(task.path() / "status");
    std::string line;
    while (std::getline(status, line)) {
      if (line.rfind("voluntary_ctxt_switches:", 0) == 0) {
        return std::stol(line.substr(line.find(':') + 1));
      }
    }
  }
  return -1;
}
#endif

} // namespace

TEST_CASE("client logger keeps order, flushes, and sleeps when idle",
          "[logger]") {
  const fs::path dir =
      fs::temp_directory_path() /
      ("eon-logger-" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()));
  fs::create_directories(dir);
  const fs::path cwd = fs::current_path();
  fs::current_path(dir);
  quill::Logger *log = eonc::log::init_client();
  fs::current_path(cwd);
  REQUIRE(log != nullptr);

  constexpr int lines = 500;
  for (int i = 0; i < lines; ++i) {
    QUILL_LOG_INFO(log, "line {}", i);
  }
  const auto start = std::chrono::steady_clock::now();
  log->flush_log();
  const auto flushed = std::chrono::steady_clock::now() - start;
  // One idle sleep at most before the backend sees the flush.
  REQUIRE(flushed < std::chrono::milliseconds(500));

  std::ifstream in(dir / "client_quill.log");
  std::vector<std::string> got;
  for (std::string line; std::getline(in, line);) {
    if (line.rfind("line ", 0) == 0) {
      got.push_back(line);
    }
  }
  REQUIRE(got.size() == static_cast<size_t>(lines));
  for (int i = 0; i < lines; ++i) {
    REQUIRE(got[static_cast<size_t>(i)] == "line " + std::to_string(i));
  }

#if defined(__linux__)
  // Idle, the backend wakes about once per kBackendSleep: some 200 times
  // in 200 ms at 1 ms. A 10 us sleep woke it some 2,600 times.
  const long before = backendSwitches();
  REQUIRE(before >= 0);
  std::this_thread::sleep_for(std::chrono::milliseconds(200));
  const long wakeups = backendSwitches() - before;
  CAPTURE(wakeups);
  REQUIRE(wakeups < 600);
#endif
  quill::Backend::stop();
  std::error_code ec;
  fs::remove_all(dir, ec);
}

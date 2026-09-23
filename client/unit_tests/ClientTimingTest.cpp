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

#include "catch2/catch_amalgamated.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#ifndef _WIN32
#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace {

double readTimeSeconds(const fs::path &path) {
  std::ifstream in(path);
  if (!in) {
    return -1.0;
  }
  std::string line;
  double found = -1.0;
  while (std::getline(in, line)) {
    std::istringstream row(line);
    double value = 0.0;
    std::string key;
    if (row >> value >> key && key == "time_seconds") {
      found = value;
    }
  }
  return found;
}

void writeText(const fs::path &path, const std::string &text) {
  std::ofstream out(path);
  REQUIRE(out.is_open());
  out << text;
}

struct CwdGuard {
  fs::path previous;
  explicit CwdGuard(const fs::path &next)
      : previous(fs::current_path()) {
    fs::current_path(next);
  }
  ~CwdGuard() {
    std::error_code ec;
    fs::current_path(previous, ec);
  }
};

} // namespace

TEST_CASE("bundled jobs record their own time_seconds", "[client][timing]") {
#ifndef _WIN32
  const char *client = std::getenv("EONCLIENT");
  const char *systems = std::getenv("EON_TEST_SYSTEMS_DIR");
  if (client == nullptr || systems == nullptr) {
    SKIP("EONCLIENT and EON_TEST_SYSTEMS_DIR are set by meson test");
  }

  const fs::path reactant = fs::path(systems) / "neb_morse" / "reactant.con";
  REQUIRE(fs::is_regular_file(reactant));

  const fs::path dir =
      fs::temp_directory_path() /
      ("eon_job_wall_" + std::to_string(static_cast<long long>(::getpid())));
  std::error_code ec;
  fs::remove_all(dir, ec);
  fs::create_directories(dir);
  fs::copy_file(reactant, dir / "pos_0.con");
  fs::copy_file(reactant, dir / "pos_1.con");

  // Slot 0 is a real dynamics run. Slot 1 is one force evaluation. A clock
  // started before the bundle loop makes the second time_seconds larger.
  writeText(dir / "config_0.ini", "[Main]\n"
                                  "job = dynamics\n"
                                  "temperature = 300\n"
                                  "quiet = true\n"
                                  "random_seed = 42\n"
                                  "[Potential]\n"
                                  "potential = lj\n"
                                  "[Dynamics]\n"
                                  "time_step = 1.0\n"
                                  "time = 10000\n");
  writeText(dir / "config_1.ini", "[Main]\n"
                                  "job = point\n"
                                  "quiet = true\n"
                                  "random_seed = 42\n"
                                  "[Potential]\n"
                                  "potential = lj\n");

  const fs::path log = dir / "spawn.log";
  const std::string log_name = log.string();
  extern char **environ;
  char *argv[] = {const_cast<char *>(client), nullptr};
  posix_spawn_file_actions_t actions;
  posix_spawn_file_actions_init(&actions);
  posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, log_name.c_str(),
                                   O_WRONLY | O_CREAT | O_TRUNC, 0644);
  posix_spawn_file_actions_adddup2(&actions, STDOUT_FILENO, STDERR_FILENO);

  pid_t pid = 0;
  int spawned = 0;
  int status = 1;
  {
    const CwdGuard cwd(dir);
    spawned = posix_spawn(&pid, client, &actions, nullptr, argv, environ);
    if (spawned == 0) {
      if (waitpid(pid, &status, 0) < 0) {
        status = 1;
      }
    }
  }
  posix_spawn_file_actions_destroy(&actions);

  if (spawned != 0 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    std::ifstream in(log);
    std::stringstream buffer;
    buffer << in.rdbuf();
    INFO(buffer.str());
    REQUIRE(spawned == 0);
    REQUIRE(WIFEXITED(status));
    REQUIRE(WEXITSTATUS(status) == 0);
  }

  const double first = readTimeSeconds(dir / "results_0.dat");
  const double second = readTimeSeconds(dir / "results_1.dat");
  INFO("time_seconds first=" << first << " second=" << second);
  REQUIRE(first > 0.0);
  REQUIRE(second > 0.0);
  REQUIRE(second < first);
  fs::remove_all(dir, ec);
#else
  SKIP("eonclient spawn covers the POSIX client");
#endif
}

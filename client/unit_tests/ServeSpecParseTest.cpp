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
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/ServeMode.h"

#include <chrono>
#include <cstdlib>
#include <string>
#include <vector>

#ifndef _WIN32
#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

static eonc::helpers::test::QuillTestLogger _quill_setup;

#ifndef _WIN32
namespace {

struct SpawnResult {
  bool running = false;
  int status = -1;
  std::string err;
};

// argv is built before fork. The child only dup2/execs, so a live logger
// thread in the parent cannot deadlock malloc in the child.
SpawnResult spawnClient(const std::vector<std::string> &args, int wait_ms) {
  SpawnResult out;
  const char *bin = std::getenv("EONCLIENT");
  if (bin == nullptr || bin[0] == '\0') {
    return out;
  }
  int fds[2];
  if (pipe(fds) != 0) {
    return out;
  }
  std::vector<std::string> storage;
  storage.reserve(args.size() + 1);
  storage.emplace_back(bin);
  storage.insert(storage.end(), args.begin(), args.end());
  std::vector<char *> av;
  av.reserve(storage.size() + 1);
  for (auto &s : storage) {
    av.push_back(s.data());
  }
  av.push_back(nullptr);

  pid_t pid = fork();
  if (pid < 0) {
    close(fds[0]);
    close(fds[1]);
    return out;
  }
  if (pid == 0) {
    dup2(fds[1], STDERR_FILENO);
    close(fds[0]);
    close(fds[1]);
    execv(bin, av.data());
    _exit(127);
  }
  close(fds[1]);
  int flags = fcntl(fds[0], F_GETFL, 0);
  if (flags >= 0) {
    fcntl(fds[0], F_SETFL, flags | O_NONBLOCK);
  }

  bool exited = false;
  const auto start = std::chrono::steady_clock::now();
  while (true) {
    char buf[512];
    ssize_t n = read(fds[0], buf, sizeof(buf));
    if (n > 0) {
      out.err.append(buf, static_cast<size_t>(n));
    }
    int status = 0;
    pid_t r = waitpid(pid, &status, WNOHANG);
    if (r == pid) {
      exited = true;
      out.status = status;
      break;
    }
    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                             std::chrono::steady_clock::now() - start)
                             .count();
    if (elapsed >= wait_ms) {
      break;
    }
    pollfd pfd{};
    pfd.fd = fds[0];
    pfd.events = POLLIN;
    auto remain = static_cast<int>(wait_ms - elapsed);
    if (remain > 50) {
      remain = 50;
    }
    poll(&pfd, 1, remain);
  }
  if (!exited) {
    kill(pid, SIGKILL);
    int status = 0;
    waitpid(pid, &status, 0);
    out.running = true;
    out.status = status;
  } else {
    char buf[512];
    ssize_t n = 0;
    while ((n = read(fds[0], buf, sizeof(buf))) > 0) {
      out.err.append(buf, static_cast<size_t>(n));
    }
  }
  close(fds[0]);
  return out;
}

} // namespace

TEST_CASE("-p without serve flags does not serve", "[serve]") {
  const char *bin = std::getenv("EONCLIENT");
  if (bin == nullptr || bin[0] == '\0') {
    SKIP("EONCLIENT is not set");
  }

  auto missing = spawnClient({"-p", "lj"}, 1500);
  INFO(missing.err);
  REQUIRE_FALSE(missing.running);
  REQUIRE(WIFEXITED(missing.status));
  CHECK(WEXITSTATUS(missing.status) == EXIT_FAILURE);
  CHECK(missing.err.find("con file") != std::string::npos);

  auto with_file = spawnClient({"-p", "lj", "no-such-structure.con"}, 1500);
  INFO(with_file.err);
  REQUIRE_FALSE(with_file.running);
  REQUIRE(WIFEXITED(with_file.status));
  CHECK(WEXITSTATUS(with_file.status) == EXIT_FAILURE);
  CHECK(with_file.err.find("Failed to load") != std::string::npos);

  auto serving = spawnClient({"-p", "lj", "--serve-port", "45921"}, 800);
  INFO(serving.err);
  CHECK(serving.running);
}
#endif

TEST_CASE("parseServeSpec single endpoint", "[serve]") {
  auto eps = parseServeSpec("lj:12345");
  REQUIRE(eps.size() == 1);
  CHECK(eps[0].potential == PotType::LJ);
  CHECK(eps[0].host == "localhost");
  CHECK(eps[0].port == 12345);
}

TEST_CASE("parseServeSpec multiple endpoints", "[serve]") {
  auto eps = parseServeSpec("lj:12345,eam_al:12346");
  REQUIRE(eps.size() == 2);
  CHECK(eps[0].potential == PotType::LJ);
  CHECK(eps[0].port == 12345);
  CHECK(eps[1].potential == PotType::EAM_AL);
  CHECK(eps[1].port == 12346);
}

TEST_CASE("parseServeSpec with host", "[serve]") {
  auto eps = parseServeSpec("lj:0.0.0.0:9999");
  REQUIRE(eps.size() == 1);
  CHECK(eps[0].potential == PotType::LJ);
  CHECK(eps[0].host == "0.0.0.0");
  CHECK(eps[0].port == 9999);
}

TEST_CASE("parseServeSpec mixed format", "[serve]") {
  auto eps = parseServeSpec("lj:12345, eam_al:0.0.0.0:12346");
  REQUIRE(eps.size() == 2);
  CHECK(eps[0].potential == PotType::LJ);
  CHECK(eps[0].host == "localhost");
  CHECK(eps[0].port == 12345);
  CHECK(eps[1].potential == PotType::EAM_AL);
  CHECK(eps[1].host == "0.0.0.0");
  CHECK(eps[1].port == 12346);
}

TEST_CASE("parseServeSpec unknown potential is skipped", "[serve]") {
  auto eps = parseServeSpec("nonexistent:12345");
  CHECK(eps.empty());
}

TEST_CASE("parseServeSpec empty spec", "[serve]") {
  auto eps = parseServeSpec("");
  CHECK(eps.empty());
}

TEST_CASE("parseServeSpec whitespace handling", "[serve]") {
  auto eps = parseServeSpec("  lj : 12345 , eam_al : 12346  ");
  // "lj " will be trimmed, ": 12345" -- the port parsing needs the colon
  // The actual parsing trims the token but the colon position is found first.
  // "  lj : 12345 " -> trimmed -> "lj : 12345"
  // first_colon at 2, pot_str="lj", rest=" 12345"
  // Since port is " 12345", stoi skips leading whitespace.
  REQUIRE(eps.size() == 2);
  CHECK(eps[0].potential == PotType::LJ);
  CHECK(eps[0].port == 12345);
}

TEST_CASE("parseServeSpec case insensitive", "[serve]") {
  auto eps = parseServeSpec("LJ:12345");
  REQUIRE(eps.size() == 1);
  CHECK(eps[0].potential == PotType::LJ);
}

TEST_CASE("ServeEndpoint struct members", "[serve]") {
  ServeEndpoint ep;
  ep.potential = PotType::LJ;
  ep.host = "localhost";
  ep.port = 12345;
  CHECK(ep.potential == PotType::LJ);
  CHECK(ep.host == "localhost");
  CHECK(ep.port == 12345);
}

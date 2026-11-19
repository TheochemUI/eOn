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

// Thread counts of a real eonclient process search, sampled from outside
// through /proc, plus the manual Quill backend that makes one thread
// possible. Linux only: the counts come from /proc/<pid>/task.

#include "catch2/catch_amalgamated.hpp"

#include "eon/EonLogger.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iterator>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#ifdef __linux__
#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

#ifdef __linux__

extern char **environ;

namespace {

constexpr const char *kQuillThreadName = "QuillBackend";

std::string slurp(const fs::path &path) {
  std::ifstream in(path);
  std::stringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

std::vector<std::string> lines(const std::string &text) {
  std::vector<std::string> out;
  std::istringstream in(text);
  for (std::string line; std::getline(in, line);) {
    out.push_back(line);
  }
  return out;
}

void writeText(const fs::path &path, const std::string &text) {
  std::ofstream out(path);
  REQUIRE(out.is_open());
  out << text;
}

size_t ownThreadCount() {
  return static_cast<size_t>(std::distance(
      fs::directory_iterator("/proc/self/task"), fs::directory_iterator{}));
}

struct CwdGuard {
  fs::path previous;
  explicit CwdGuard(const fs::path &next) : previous(fs::current_path()) {
    fs::current_path(next);
  }
  ~CwdGuard() {
    std::error_code ec;
    fs::current_path(previous, ec);
  }
};

struct ClientRun {
  int exitCode{-1};
  long samples{0};
  pid_t pid{0};
  // Every thread id seen, with every name it carried while it was seen.
  std::map<long, std::set<std::string>> threads;
  fs::path dir;

  // Threads that are neither the compute thread nor Quill's backend.
  [[nodiscard]] std::vector<long> workers() const {
    std::vector<long> out;
    for (const auto &[tid, names] : threads) {
      if (tid != pid && names.count(kQuillThreadName) == 0) {
        out.push_back(tid);
      }
    }
    return out;
  }
  [[nodiscard]] size_t quillThreads() const {
    size_t n = 0;
    for (const auto &entry : threads) {
      n += entry.second.count(kQuillThreadName);
    }
    return n;
  }
  [[nodiscard]] std::string describe() const {
    std::ostringstream out;
    out << "samples=" << samples << " threads=" << threads.size() << " [";
    for (const auto &[tid, names] : threads) {
      out << " " << (tid == pid ? "main" : std::to_string(tid)) << ":";
      for (const auto &name : names) {
        out << name << ",";
      }
    }
    out << " ]";
    return out.str();
  }
};

// Runs eonclient in `dir` and samples /proc/<pid>/task until it exits.
ClientRun runClient(const fs::path &dir) {
  const char *client = std::getenv("EONCLIENT");
  REQUIRE(client != nullptr);

  ClientRun run;
  run.dir = dir;
  const std::string logName = (dir / "stdout.log").string();
  char *argv[] = {const_cast<char *>(client), nullptr};
  posix_spawn_file_actions_t actions;
  posix_spawn_file_actions_init(&actions);
  posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, logName.c_str(),
                                   O_WRONLY | O_CREAT | O_TRUNC, 0644);
  posix_spawn_file_actions_adddup2(&actions, STDOUT_FILENO, STDERR_FILENO);

  pid_t pid = 0;
  int spawned = 0;
  {
    const CwdGuard cwd(dir);
    spawned = posix_spawn(&pid, client, &actions, nullptr, argv, environ);
  }
  posix_spawn_file_actions_destroy(&actions);
  REQUIRE(spawned == 0);
  run.pid = pid;

  const fs::path taskDir = fs::path("/proc") / std::to_string(pid) / "task";
  const timespec pause{0, 100000};
  int status = 0;
  while (true) {
    const pid_t done = waitpid(pid, &status, WNOHANG);
    if (done == pid) {
      break;
    }
    if (done < 0) {
      FAIL("waitpid failed while sampling the client");
    }
    std::error_code ec;
    for (fs::directory_iterator it(taskDir, ec), end; !ec && it != end;
         it.increment(ec)) {
      const long tid = std::atol(it->path().filename().c_str());
      std::ifstream comm(it->path() / "comm");
      std::string name;
      if (tid > 0 && std::getline(comm, name)) {
        run.threads[tid].insert(name);
      }
    }
    ++run.samples;
    nanosleep(&pause, nullptr);
  }
  run.exitCode = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
  return run;
}

fs::path makeJobDir(const char *tag, bool withDisplacement) {
  const char *systems = std::getenv("EON_TEST_SYSTEMS_DIR");
  REQUIRE(systems != nullptr);
  const fs::path source = fs::path(systems) / "saddle_search";
  REQUIRE(fs::is_regular_file(source / "pos.con"));

  const fs::path dir =
      fs::temp_directory_path() /
      (std::string("eon_threads_") + tag + "_" + std::to_string(::getpid()));
  std::error_code ec;
  fs::remove_all(dir, ec);
  fs::create_directories(dir);
  fs::copy_file(source / "pos.con", dir / "pos.con");
  if (withDisplacement) {
    fs::copy_file(source / "displacement.con", dir / "displacement.con");
    fs::copy_file(source / "direction.dat", dir / "direction.dat");
  }
  return dir;
}

void writeProcessSearchConfig(const fs::path &dir, const std::string &main) {
  writeText(dir / "config.ini", "[Main]\n"
                                "job = process_search\n"
                                "temperature = 300\n"
                                "random_seed = 706253457\n" +
                                    main +
                                    "\n[Potential]\n"
                                    "potential = morse_pt\n"
                                    "\n[Optimizer]\n"
                                    "opt_method = lbfgs\n"
                                    "converged_force = 0.001\n"
                                    "max_iterations = 1000\n"
                                    "max_move = 0.2\n"
                                    "\n[Saddle Search]\n"
                                    "displace_least_coordinated_weight = 1.0\n"
                                    "displace_radius = 3.3\n"
                                    "displace_magnitude = 0.01\n"
                                    "min_mode_method = dimer\n"
                                    "max_energy = 10.0\n");
}

bool hasResult(const fs::path &dir, const std::string &line) {
  for (const auto &row : lines(slurp(dir / "results.dat"))) {
    if (row == line) {
      return true;
    }
  }
  return false;
}

// Positions of `needles` in `text`, which must appear in the given order.
void requireInOrder(const std::string &text,
                    const std::vector<std::string> &needles) {
  size_t from = 0;
  for (const auto &needle : needles) {
    const size_t at = text.find(needle, from);
    INFO("missing or out of order: " << needle);
    REQUIRE(at != std::string::npos);
    from = at + needle.size();
  }
}

const std::vector<std::string> kRunOrder = {
    "eOn Client",   "Loading parameter file config.ini",
    "Timing Information:", "Real time:", "User time:", "System time:"};

} // namespace

TEST_CASE("a single threaded client runs the process search on one thread",
          "[client][threads][process_search]") {
  const fs::path dir = makeJobDir("single", true);
  writeProcessSearchConfig(dir, "parallel = false\n"
                                "single_threaded_client = true");
  const ClientRun run = runClient(dir);
  INFO(run.describe());
  INFO(slurp(dir / "stdout.log"));

  REQUIRE(run.exitCode == 0);
  REQUIRE(hasResult(dir, "0 termination_reason"));
  // The job ran for long enough to be sampled many times.
  REQUIRE(run.samples > 100);
  REQUIRE(run.threads.size() == 1);
  REQUIRE(run.threads.count(run.pid) == 1);
  REQUIRE(run.threads.at(run.pid) == std::set<std::string>{"eonclient"});

  // Nothing is lost on the way out: the last record of the run is in the
  // file and the console copy.
  const std::string log = slurp(dir / "client_quill.log");
  requireInOrder(log, kRunOrder);
  requireInOrder(slurp(dir / "stdout.log"), kRunOrder);

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("parallel = false leaves the logger thread and no worker",
          "[client][threads][process_search]") {
  const fs::path dir = makeJobDir("serial", true);
  writeProcessSearchConfig(dir, "parallel = false");
  const ClientRun run = runClient(dir);
  INFO(run.describe());

  REQUIRE(run.exitCode == 0);
  REQUIRE(hasResult(dir, "0 termination_reason"));
  REQUIRE(run.samples > 100);
  REQUIRE(run.workers().empty());
  REQUIRE(run.quillThreads() == 1);
  REQUIRE(run.threads.size() == 2);

  const std::string log = slurp(dir / "client_quill.log");
  requireInOrder(log, kRunOrder);

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("both logging modes write the same records",
          "[client][threads][process_search]") {
  const fs::path single = makeJobDir("same_single", true);
  const fs::path threaded = makeJobDir("same_thread", true);
  writeProcessSearchConfig(single, "parallel = false\n"
                                   "single_threaded_client = true");
  writeProcessSearchConfig(threaded, "parallel = false");
  REQUIRE(runClient(single).exitCode == 0);
  REQUIRE(runClient(threaded).exitCode == 0);

  const auto a = lines(slurp(single / "client_quill.log"));
  const auto b = lines(slurp(threaded / "client_quill.log"));
  REQUIRE(a.size() > 10);
  REQUIRE(a.size() == b.size());
  // Timings, the pid, and the directory differ between runs. The text of
  // every other record is the same, in the same order.
  size_t same = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    same += a[i] == b[i] ? 1 : 0;
  }
  REQUIRE(same + 8 >= a.size());

  std::error_code ec;
  fs::remove_all(single, ec);
  fs::remove_all(threaded, ec);
}

TEST_CASE("a single threaded client still reports a failed job",
          "[client][threads][process_search]") {
  // No displacement.con and no direction.dat: the search cannot start.
  const fs::path dir = makeJobDir("failed", false);
  writeProcessSearchConfig(dir, "parallel = false\n"
                                "single_threaded_client = true\n"
                                "\n[Saddle Search]\n"
                                "displace_type = load");
  const ClientRun run = runClient(dir);
  INFO(run.describe());
  const std::string out = slurp(dir / "stdout.log");
  INFO(out);

  REQUIRE(run.exitCode != 0);
  REQUIRE(run.threads.size() == 1);
  const std::string cause = "failed to load displacement.con";
  REQUIRE(out.find(cause) != std::string::npos);
  const std::string log = slurp(dir / "client_quill.log");
  REQUIRE(log.find("[ERROR] unhandled exception: " + cause) !=
          std::string::npos);
  requireInOrder(log, {"Loading parameter file config.ini",
                       "[ERROR] unhandled exception"});

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("single_threaded_client with parallel left on is refused",
          "[client][threads]") {
  const fs::path dir = makeJobDir("refused", true);
  writeProcessSearchConfig(dir, "single_threaded_client = true");
  const ClientRun run = runClient(dir);
  INFO(run.describe());
  const std::string out = slurp(dir / "stdout.log");
  INFO(out);

  REQUIRE(run.exitCode != 0);
  REQUIRE(out.find("single_threaded_client = true needs parallel = false") !=
          std::string::npos);
  REQUIRE_FALSE(fs::exists(dir / "results.dat"));

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("the manual Quill backend drains on the logging thread",
          "[logger][threads]") {
  const fs::path file =
      fs::temp_directory_path() /
      ("eon_manual_backend_" + std::to_string(::getpid()) + ".log");
  std::error_code ec;
  fs::remove(file, ec);

  quill::FileSinkConfig sinkConfig;
  sinkConfig.set_open_mode('w');
  auto sink = quill::Frontend::create_or_get_sink<quill::FileSink>(
      file.string(), sinkConfig);
  auto *logger = quill::Frontend::create_or_get_logger(
      "manual_probe", std::move(sink),
      quill::PatternFormatterOptions{"%(message)"},
      quill::ClockSourceType::System);
  quill::Frontend::preallocate();

  const size_t before = ownThreadCount();
  REQUIRE(eonc::log::start_client_backend(eonc::log::BackendMode::Manual) ==
          eonc::log::BackendMode::Manual);
  REQUIRE(eonc::log::manual_backend());
  REQUIRE(ownThreadCount() == before);

  constexpr int kLines = 500;
  for (int i = 0; i < kLines; ++i) {
    QUILL_LOG_INFO(logger, "line {}", i);
  }
  eonc::log::flush(logger);
  {
    const auto written = lines(slurp(file));
    REQUIRE(written.size() == kLines);
    for (int i = 0; i < kLines; ++i) {
      REQUIRE(written[static_cast<size_t>(i)] == "line " + std::to_string(i));
    }
  }
  REQUIRE(ownThreadCount() == before);

  // Polling is timed. Before the interval is up a record stays queued.
  QUILL_LOG_INFO(logger, "timed");
  std::this_thread::sleep_for(eonc::log::detail::kManualPollInterval +
                              std::chrono::milliseconds(50));
  eonc::log::poll_if_due();
  REQUIRE(lines(slurp(file)).back() == "timed");

  // Another thread cannot poll. It returns instead of waiting for a backend
  // that does not exist, and the record stays queued for the owner.
  QUILL_LOG_INFO(logger, "queued");
  std::thread other([&] { eonc::log::flush(logger); });
  other.join();
  REQUIRE(lines(slurp(file)).back() == "timed");
  eonc::log::flush(logger);
  REQUIRE(lines(slurp(file)).back() == "queued");

  // The last records written before exit are not dropped.
  QUILL_LOG_INFO(logger, "last");
  eonc::log::shutdown();
  REQUIRE(lines(slurp(file)).back() == "last");
  REQUIRE(ownThreadCount() == before);

  fs::remove(file, ec);
}

#else

TEST_CASE("process search thread counts need /proc", "[client][threads]") {
  SKIP("thread counts are read from /proc/<pid>/task");
}

#endif

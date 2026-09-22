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

#include "eon/Bundling.h"

#include "catch2/catch_amalgamated.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace {

struct CwdGuard {
  std::filesystem::path previous;
  explicit CwdGuard(const std::filesystem::path &next)
      : previous(std::filesystem::current_path()) {
    std::filesystem::current_path(next);
  }
  ~CwdGuard() {
    std::error_code ec;
    std::filesystem::current_path(previous, ec);
  }
};

void writeBytes(const std::filesystem::path &path, const std::string &body) {
  std::ofstream out(path, std::ios::binary);
  REQUIRE(out);
  out << body;
}

std::string readBytes(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  REQUIRE(in);
  return std::string(std::istreambuf_iterator<char>(in), {});
}

bool containsName(const std::vector<std::string> &names,
                  const std::string &want) {
  return std::find(names.begin(), names.end(), want) != names.end();
}

} // namespace

TEST_CASE("unbundle restores compressed names that bundle writes",
          "[bundling]") {
  namespace fs = std::filesystem;
  const auto dir = fs::temp_directory_path() / "eon_bundle_con_gz";
  fs::remove_all(dir);
  fs::create_directories(dir);
  writeBytes(dir / "results.con", "con-bytes");
  writeBytes(dir / "results.con.gz", "gz-bytes");
  writeBytes(dir / "results.dat", "dat-bytes");
  writeBytes(dir / "pos_final.con", "leave-me");

  {
    CwdGuard guard(dir);
    std::vector<std::string> bundled;
    eonc::bundle(3, {"results.con", "results.con.gz", "results.dat"}, &bundled);
    REQUIRE(fs::is_regular_file("results_3.con"));
    REQUIRE(fs::is_regular_file("results_3.con.gz"));
    REQUIRE(fs::is_regular_file("results_3.dat"));
    REQUIRE_FALSE(fs::exists("results.con.gz"));
    REQUIRE(fs::is_regular_file("pos_final.con"));

    const auto restored = eonc::unbundle(3);
    REQUIRE(restored.size() == 3);
    REQUIRE(containsName(restored, "results.con"));
    REQUIRE(containsName(restored, "results.con.gz"));
    REQUIRE(containsName(restored, "results.dat"));
    REQUIRE_FALSE(containsName(restored, "pos_final.con"));
    REQUIRE(fs::is_regular_file("results.con.gz"));
    REQUIRE(readBytes("results.con.gz") == "gz-bytes");
    REQUIRE(readBytes("results.con") == "con-bytes");
    REQUIRE(readBytes("results.dat") == "dat-bytes");
    REQUIRE(readBytes("pos_final.con") == "leave-me");
    REQUIRE_FALSE(fs::exists("pos.con"));

    REQUIRE(eonc::unbundle(2).empty());
  }

  fs::remove_all(dir);
}

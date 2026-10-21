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

#include "ReadconDbMirror.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/ConFileIO.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"

#include <cstdlib>
#include <string>

namespace tests {
namespace {

void set_library_override(const char *path) {
#if defined(_WIN32)
  _putenv_s("EON_READCON_DB_LIBRARY", path == nullptr ? "" : path);
#else
  if (path == nullptr) {
    unsetenv("EON_READCON_DB_LIBRARY");
  } else {
    setenv("EON_READCON_DB_LIBRARY", path, 1);
  }
#endif
}

} // namespace

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("loaded rkrdb_open matches the python package version",
          "[readcon-db]") {
  set_library_override(nullptr);
  eonc::io::readcon_db_mirror_reset();
  REQUIRE(eonc::io::readcon_db_mirror_ok());
  REQUIRE(std::string(eonc::io::readcon_db_loaded_version()) == "0.1.6");
}

TEST_CASE("missing readcon-db leaves IoStatus unchanged", "[readcon-db]") {
#if defined(_WIN32)
  set_library_override("C:/no/such/readcon_db.dll");
#else
  set_library_override("/no/such/libreadcon_db.so");
#endif
  eonc::io::readcon_db_mirror_reset();

  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  Matter matter(pot, params);
  const eonc::io::IoStatus status =
      matter.con2matter(std::string("reactant.con"));
  REQUIRE(status == eonc::io::IoStatus::Ok);
  REQUIRE_FALSE(eonc::io::readcon_db_mirror_ok());
}

} // namespace tests

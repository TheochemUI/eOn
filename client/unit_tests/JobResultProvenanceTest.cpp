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
#include "eon/JobResult.h"

TEST_CASE("legacy optimizer provenance omits xtsci and eindir stamps",
          "[jobresult][provenance]") {
  eonc::JobResultEnvelope env;
  env.job_type = "minimization";
  env.provenance.backend = "cg";
  env.provenance.engine_version = "3.3.1";
  env.provenance.engine_build_identity = "abc123";
  env.provenance.rgpot_name = "lj";
  const std::string text = env.toString();
  REQUIRE(text.find("cg optimizer_backend\n") != std::string::npos);
  REQUIRE(text.find("eon.optimizer.v1 optimizer_provenance_schema\n") !=
          std::string::npos);
  REQUIRE(text.find("eon engine_id\n") != std::string::npos);
  REQUIRE(text.find("lj rgpot_name\n") != std::string::npos);
  REQUIRE(text.find("optimizer_xts_abi_major") == std::string::npos);
  REQUIRE(text.find("optimizer_eindir_abi_major") == std::string::npos);
  REQUIRE(text.find("compatibility_engine_protocol_family") ==
          std::string::npos);
}

TEST_CASE("xtsci optimizer provenance records xts and eindir ABI",
          "[jobresult][provenance]") {
  eonc::JobResultProvenance provenance;
  provenance.backend = "xtsci";
  provenance.xtsci = true;
  provenance.xts_abi_major = 1;
  provenance.xts_abi_minor = 10;
  provenance.xts_abi_layout = 2;
  provenance.engine_version = "3.3.1";
  provenance.engine_build_identity = "abc123";
  provenance.rgpot_name = "morse_pt";
  provenance.rgpot_version = std::string(eonc::JobResultProvenance::rgpot_pin);
  provenance.has_eindir = true;
  provenance.eindir_abi_major = 1;
  provenance.eindir_abi_minor = 0;
  provenance.eindir_objective_layout = 1;
  provenance.eindir_objective_size = 64;
  provenance.eindir_objective_align = 8;
  provenance.eindir_dlpack_major = 1;
  provenance.eindir_dlpack_minor = 0;
  provenance.eindir_features = 3;
  const std::string text = provenance.text();
  REQUIRE(text.find("xtsci optimizer_backend\n") != std::string::npos);
  REQUIRE(text.find("1 optimizer_xts_abi_major\n") != std::string::npos);
  REQUIRE(text.find("10 optimizer_xts_abi_minor\n") != std::string::npos);
  REQUIRE(text.find("2 optimizer_xts_abi_layout\n") != std::string::npos);
  REQUIRE(text.find("eon.objective compatibility_engine_protocol_family\n") !=
          std::string::npos);
  REQUIRE(text.find("1 optimizer_eindir_abi_major\n") != std::string::npos);
  REQUIRE(text.find("64 optimizer_eindir_objective_size\n") !=
          std::string::npos);
  REQUIRE(text.find("3 optimizer_eindir_features\n") != std::string::npos);
  REQUIRE(text.find("eon.rgpot.v1 rgpot_schema\n") != std::string::npos);
  REQUIRE(text.find("3.2.0 rgpot_version\n") != std::string::npos);
  REQUIRE(text.find("0.14.9 compatibility_readcon_min_version\n") !=
          std::string::npos);
}

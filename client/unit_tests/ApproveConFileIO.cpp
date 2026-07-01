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
#include "ConFileIO.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <string>

namespace {

namespace fs = std::filesystem;

static eonc::helpers::test::QuillTestLogger _quill_setup;

/// Unique temp .con path (avoids parallel/manual temp-dir name clashes).
fs::path make_tmp_con(std::string_view tag) {
  const auto stamp =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();
  return fs::temp_directory_path() /
         std::format("approve_confileio_{}_{}.con", tag, stamp);
}

std::string dump_geometry(const Matter &m) {
  std::ostringstream os;
  os << std::format("n_atoms={}\n", m.numberOfAtoms());
  const auto cell = m.getCell();
  os << std::format(
      "cell={:.8f},{:.8f},{:.8f};{:.8f},{:.8f},{:.8f};{:.8f},{:.8f},{:.8f}\n",
      cell(0, 0), cell(0, 1), cell(0, 2), cell(1, 0), cell(1, 1), cell(1, 2),
      cell(2, 0), cell(2, 1), cell(2, 2));
  const auto pos = m.getPositions();
  for (long i = 0; i < m.numberOfAtoms(); ++i) {
    os << std::format(
        "a {} Z={} mass={:.8f} fixed={} id={} pos={:.8f},{:.8f},{:.8f}\n", i,
        m.getAtomicNr(i), m.getMass(i), m.getFixed(i), m.getAtomIndex(i),
        pos(i, 0), pos(i, 1), pos(i, 2));
  }
  return os.str();
}

std::string read_master(const std::string &name) {
  // Meson workdir is neb_morse; masters live in approval_tests/ there.
  const fs::path p = fs::path("approval_tests") / name;
  std::ifstream in(p);
  REQUIRE(in.good());
  return std::string((std::istreambuf_iterator<char>(in)),
                     std::istreambuf_iterator<char>());
}

std::shared_ptr<Matter> load_reactant() {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(m->con2matter(std::string("reactant.con"))));
  return m;
}

size_t count_substr(const std::string &hay, const std::string &needle) {
  size_t n = 0;
  size_t pos = 0;
  while ((pos = hay.find(needle, pos)) != std::string::npos) {
    ++n;
    pos += needle.size();
  }
  return n;
}

} // namespace

TEST_CASE("VerifyClassicFixtureRead", "[approval][confileio][compat]") {
  auto m = load_reactant();
  REQUIRE(
      dump_geometry(*m) ==
      read_master("ApproveConFileIO.VerifyClassicFixtureRead.approved.txt"));
}

TEST_CASE("VerifyGeometryOnlyRoundTrip", "[approval][confileio][compat]") {
  auto m = load_reactant();
  REQUIRE(m->needsForceUpdate());

  const auto tmp = make_tmp_con("geom_roundtrip");
  REQUIRE(eonc::io::io_ok(m->matter2con(tmp.string(), false)));

  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m2 = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(m2->con2matter(tmp.string())));
  REQUIRE(
      dump_geometry(*m2) ==
      read_master("ApproveConFileIO.VerifyGeometryOnlyRoundTrip.approved.txt"));
  fs::remove(tmp);
}

TEST_CASE("VerifyForceBearingWrite", "[approval][confileio][modern]") {
  auto m = load_reactant();
  const double E = m->getPotentialEnergy();
  REQUIRE_FALSE(m->needsForceUpdate());

  const auto tmp = make_tmp_con("forces_energy");
  REQUIRE(eonc::io::io_ok(m->matter2con(tmp.string(), false)));

  std::ifstream in(tmp);
  std::string body((std::istreambuf_iterator<char>(in)),
                   std::istreambuf_iterator<char>());

  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m2 = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(m2->con2matter(tmp.string())));
  fs::remove(tmp);

  std::ostringstream summary;
  summary << "energy_finite=" << (std::isfinite(E) ? "true" : "false") << "\n";
  summary << "has_forces_token="
          << (body.find("forces") != std::string::npos ? "true" : "false")
          << "\n";
  summary << "has_energy_token="
          << (body.find("\"energy\"") != std::string::npos ? "true" : "false")
          << "\n";
  summary << "needs_force_update_after_reload="
          << (m2->needsForceUpdate() ? "true" : "false") << "\n";
  summary << "geometry_after_force_write:\n" << dump_geometry(*m2);
  REQUIRE(summary.str() ==
          read_master("ApproveConFileIO.VerifyForceBearingWrite.approved.txt"));
}

TEST_CASE("VerifyNebPathCloneWrite", "[approval][confileio][neb]") {
  auto a = load_reactant();
  auto b = load_reactant();
  auto c = load_reactant();
  std::vector<std::shared_ptr<Matter>> path{a, b, c};
  std::vector<eonc::io::ConFrameMetadata> metas(3);
  for (size_t i = 0; i < 3; ++i) {
    metas[i].frame_index = i;
    metas[i].neb_bead = i;
    metas[i].neb_band = 7;
    metas[i].energy = -1.0 * static_cast<double>(i);
  }
  const auto tmp = make_tmp_con("neb_band");
  REQUIRE(eonc::io::io_ok(eonc::io::writeNebPath(tmp.string(), path, metas)));

  std::ifstream in(tmp);
  std::string body((std::istreambuf_iterator<char>(in)),
                   std::istreambuf_iterator<char>());
  fs::remove(tmp);

  std::ostringstream summary;
  summary << "neb_band_markers\n";
  summary << "count_neb_band_7=" << count_substr(body, "\"neb_band\":7")
          << "\n";
  summary << "count_neb_bead=" << count_substr(body, "\"neb_bead\":") << "\n";
  summary << "has_con_spec="
          << (body.find("con_spec_version") != std::string::npos ? "true"
                                                                 : "false")
          << "\n";
  REQUIRE(summary.str() ==
          read_master("ApproveConFileIO.VerifyNebPathCloneWrite.approved.txt"));
}

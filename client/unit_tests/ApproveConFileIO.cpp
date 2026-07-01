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
#define APPROVALS_CATCH2_V3
#include "ApprovalTests.hpp"

// Masters live under workdir approval_tests/ (meson workdir = neb_morse).
auto directoryDisposer =
    ApprovalTests::Approvals::useApprovalsSubdirectory("approval_tests");

#include "ConFileIO.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <string>

namespace {

namespace fs = std::filesystem;

static eonc::helpers::test::QuillTestLogger _quill_setup;

/// Deterministic geometry dump for golden masters (main-era shared fields).
/// Omits forces/energy so classic fixtures compare with pre-0.13 dumps.
std::string dump_geometry(const Matter &m) {
  std::ostringstream os;
  os << std::format("n_atoms={}
", m.numberOfAtoms());
  const auto cell = m.getCell();
  os << std::format(
      "cell={:.8f},{:.8f},{:.8f};{:.8f},{:.8f},{:.8f};{:.8f},{:.8f},{:.8f}
",
      cell(0, 0), cell(0, 1), cell(0, 2), cell(1, 0), cell(1, 1), cell(1, 2),
      cell(2, 0), cell(2, 1), cell(2, 2));
  const auto pos = m.getPositions();
  for (long i = 0; i < m.numberOfAtoms(); ++i) {
    os << std::format(
        "a {} Z={} mass={:.8f} fixed={} id={} pos={:.8f},{:.8f},{:.8f}
        ", i,
        m.getAtomicNr(i),
        m.getMass(i), m.getFixed(i), m.getAtomIndex(i), pos(i, 0), pos(i, 1),
        pos(i, 2));
  }
  return os.str();
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
  // Back-compat: reading a main-era classic .con must yield stable geometry.
  auto m = load_reactant();
  ApprovalTests::Approvals::verify(dump_geometry(*m));
}

TEST_CASE("VerifyGeometryOnlyRoundTrip", "[approval][confileio][compat]") {
  // Dirty pot cache → no forces/energy on write (classic movie layout).
  auto m = load_reactant();
  REQUIRE(m->needsForceUpdate());

  const auto tmp =
      fs::temp_directory_path() / "approve_confileio_geom_roundtrip.con";
  REQUIRE(eonc::io::io_ok(m->matter2con(tmp.string(), false)));

  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m2 = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(m2->con2matter(tmp.string())));
  ApprovalTests::Approvals::verify(dump_geometry(*m2));
  fs::remove(tmp);
}

TEST_CASE("VerifyForceBearingWrite", "[approval][confileio][modern]") {
  // Intentional 0.13 layout (forces + energy); own golden master.
  auto m = load_reactant();
  const double E = m->getPotentialEnergy();
  REQUIRE_FALSE(m->needsForceUpdate());

  const auto tmp =
      fs::temp_directory_path() / "approve_confileio_forces_energy.con";
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

  // Exact potential energy is pot-implementation detail; lock I/O shape +
  // geometry (compat) and that energy was restored so forces stay trusted.
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
  ApprovalTests::Approvals::verify(summary.str());
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
  const auto tmp = fs::temp_directory_path() / "approve_confileio_neb_band.con";
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
  ApprovalTests::Approvals::verify(summary.str());
}

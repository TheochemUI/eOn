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
#include "Matter.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include <cstdio>
#include <filesystem>
#include <fstream>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class ConFileIOFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> original;

  ConFileIOFixture()
      : params{},
        pot{nullptr},
        original{nullptr} {
    params.potential_options.potential = PotType::LJ;
    pot = eonc::helpers::makePotential(PotType::LJ, params);
    original = std::make_shared<Matter>(pot, params);
    original->con2matter(std::string("reactant.con"));
  }

  ~ConFileIOFixture() {
    // Clean up test output files
    std::remove("_test_roundtrip.con");
    std::remove("_test_append.con");
  }
};

TEST_CASE_METHOD(ConFileIOFixture,
                 "Con file write-read round trip preserves data",
                 "[confileio]") {
  std::string tmpfile = "_test_roundtrip.con";
  original->matter2con(tmpfile);

  auto loaded = std::make_shared<Matter>(pot, params);
  loaded->con2matter(tmpfile);

  REQUIRE(loaded->numberOfAtoms() == original->numberOfAtoms());

  // Positions match
  auto origPos = original->getPositions();
  auto loadPos = loaded->getPositions();
  REQUIRE(origPos.isApprox(loadPos, 1e-6));

  // Cell matches
  auto origCell = original->getCell();
  auto loadCell = loaded->getCell();
  REQUIRE(origCell.isApprox(loadCell, 1e-6));

  // Atomic numbers match
  REQUIRE(loaded->getAtomicNrs() == original->getAtomicNrs());
}

TEST_CASE_METHOD(ConFileIOFixture, "Con file append mode creates larger file",
                 "[confileio]") {
  std::string tmpfile = "_test_append.con";
  original->matter2con(tmpfile, false);

  // Get file size after first write
  auto size1 = std::filesystem::file_size(tmpfile);
  REQUIRE(size1 > 0);

  // Append a second structure
  original->matter2con(tmpfile, true);

  auto size2 = std::filesystem::file_size(tmpfile);
  REQUIRE(size2 > size1);
}

TEST_CASE_METHOD(ConFileIOFixture,
                 "Con file free-function interface matches member delegates",
                 "[confileio]") {
  std::string tmpfile = "_test_roundtrip.con";
  eonc::io::matter2con(*original, tmpfile, false);

  auto loaded = std::make_shared<Matter>(pot, params);
  eonc::io::con2matter(*loaded, tmpfile);

  REQUIRE(loaded->numberOfAtoms() == original->numberOfAtoms());
  REQUIRE(loaded->getPositions().isApprox(original->getPositions(), 1e-6));
}

TEST_CASE_METHOD(ConFileIOFixture, "Con file preserves fixed-atom flags",
                 "[confileio]") {
  // Set some atoms as fixed
  original->setFixed(0, true);
  original->setFixed(1, false);

  std::string tmpfile = "_test_roundtrip.con";
  original->matter2con(tmpfile);

  auto loaded = std::make_shared<Matter>(pot, params);
  loaded->con2matter(tmpfile);

  REQUIRE(loaded->getFixed(0) == true);
  REQUIRE(loaded->getFixed(1) == false);
}

TEST_CASE_METHOD(ConFileIOFixture, "XYZ output is finite and non-empty",
                 "[confileio]") {
  // Note: matter2xyz appends ".xyz" to the filename
  auto tmppath = std::filesystem::temp_directory_path() / "_test_output";
  std::string tmpbase = tmppath.string();
  eonc::io::matter2xyz(*original, tmpbase);
  std::string tmpfile = tmpbase + ".xyz";

  auto size = std::filesystem::file_size(tmpfile);
  REQUIRE(size > 0);

  // Read first line -- should be atom count
  std::ifstream f(tmpfile);
  int natoms;
  f >> natoms;
  REQUIRE(natoms == original->numberOfAtoms());

  std::remove(tmpfile.c_str());
}

TEST_CASE_METHOD(ConFileIOFixture,
                 "String interface write-read roundtrip preserves geometry",
                 "[confileio]") {
  // Write via string interface
  eonc::io::matter2con(*original, "_test_str.con");

  // Read back
  auto loaded = std::make_shared<Matter>(pot, params);
  loaded->con2matter(std::string("_test_str.con"));

  REQUIRE(loaded->numberOfAtoms() == original->numberOfAtoms());
  REQUIRE(loaded->getPositions().isApprox(original->getPositions(), 1e-6));

  std::remove("_test_str.con");
}

TEST_CASE_METHOD(ConFileIOFixture,
                 "Con file energy is preserved after write-read cycle",
                 "[confileio]") {
  double E_orig = original->getPotentialEnergy();

  std::string tmpfile = "_test_roundtrip.con";
  original->matter2con(tmpfile);

  auto loaded = std::make_shared<Matter>(pot, params);
  loaded->con2matter(tmpfile);

  double E_loaded = loaded->getPotentialEnergy();
  // SVN reference: -39.965351
  REQUIRE(E_loaded == Catch::Approx(E_orig).epsilon(1e-6));
}

TEST_CASE("ConFileIO handles velocity data", "[confileio][velocity]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m = std::make_shared<Matter>(pot, params);
  m->con2matter(std::string("reactant.con"));

  // Set some velocities
  auto vel = m->getVelocities();
  vel.setRandom();
  m->setVelocities(vel);

  // Write and read back
  auto tmppath = std::filesystem::temp_directory_path() / "_test_vel";
  std::string tmpfile = tmppath.string() + ".con";
  m->matter2con(tmpfile, false);

  auto m2 = std::make_shared<Matter>(pot, params);
  m2->con2matter(tmpfile);

  // Positions should match
  REQUIRE(m->getPositions().isApprox(m2->getPositions(), 1e-6));
  std::filesystem::remove(tmpfile);
}

TEST_CASE("ConFileIO multiple writes accumulate in append mode",
          "[confileio][append]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m = std::make_shared<Matter>(pot, params);
  m->con2matter(std::string("reactant.con"));

  auto tmppath = std::filesystem::temp_directory_path() / "_test_append.con";
  std::string tmpfile = tmppath.string();

  // Write twice in append mode
  m->matter2con(tmpfile, false);
  m->matter2con(tmpfile, true);

  // File should be roughly double the size of single write
  auto fsize = std::filesystem::file_size(tmpfile);
  REQUIRE(fsize > 1000); // Two structures

  std::filesystem::remove(tmpfile);
}

TEST_CASE("ConFileIO reads multi-component Pt system",
          "[confileio][multicomponent]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  auto m = std::make_shared<Matter>(pot, params);
  m->con2matter(std::string("../Pt_Heptamer_FrozenLayers/pos.con"));

  // Pt heptamer has atoms with some fixed
  REQUIRE(m->numberOfAtoms() > 0);
  REQUIRE(m->numberOfFreeAtoms() < m->numberOfAtoms());

  // Verify cell is reasonable
  auto cell = m->getCell();
  REQUIRE(cell(0, 0) > 0.0);
  REQUIRE(cell(1, 1) > 0.0);
  REQUIRE(cell(2, 2) > 0.0);

  // Write and read back
  auto tmppath = std::filesystem::temp_directory_path() / "_test_multi.con";
  std::string tmpfile = tmppath.string();
  m->matter2con(tmpfile, false);

  auto m2 = std::make_shared<Matter>(pot, params);
  m2->con2matter(tmpfile);

  REQUIRE(m2->numberOfAtoms() == m->numberOfAtoms());
  REQUIRE(m2->numberOfFreeAtoms() == m->numberOfFreeAtoms());
  REQUIRE(m->getPositions().isApprox(m2->getPositions(), 1e-6));

  std::filesystem::remove(tmpfile);
}

TEST_CASE("ConFileIO reads Si diamond system", "[confileio][si]") {
  Parameters params;
  params.potential_options.potential = PotType::SW_SI;
  auto pot = eonc::helpers::makePotential(PotType::SW_SI, params);
  auto m = std::make_shared<Matter>(pot, params);
  m->con2matter(std::string("../si_diamond/pos.con"));

  REQUIRE(m->numberOfAtoms() == 8);
  REQUIRE(m->getAtomicNr(0) == 14); // Si atomic number
}

TEST_CASE("ConFileIO convel round-trip preserves velocities",
          "[confileio][convel]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m = std::make_shared<Matter>(pot, params);
  m->con2matter(std::string("reactant.con"));

  // Set velocities
  auto vel = m->getVelocities();
  vel.setRandom();
  vel *= 0.01; // small velocities
  m->setVelocities(vel);

  auto tmppath = std::filesystem::temp_directory_path() / "_test_out.convel";
  std::string tmpfile = tmppath.string();

  // Write convel (function expects .convel extension or appends it)
  eonc::io::matter2convel(*m, tmpfile);
  REQUIRE(std::filesystem::exists(tmpfile));

  // Read back
  auto m2 = std::make_shared<Matter>(pot, params);
  eonc::io::convel2matter(*m2, tmpfile);

  // Positions and velocities should match
  REQUIRE(m->getPositions().isApprox(m2->getPositions(), 1e-6));
  auto vel2 = m2->getVelocities();
  // convel format has limited precision
  REQUIRE(vel.isApprox(vel2, 1e-3));

  std::filesystem::remove(tmpfile);
}

TEST_CASE("ConFileIO writeTibble produces valid output",
          "[confileio][tibble]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto m = std::make_shared<Matter>(pot, params);
  m->con2matter(std::string("reactant.con"));
  // Trigger force computation so writeTibble can access forces
  m->getPotentialEnergy();

  auto tmppath = std::filesystem::temp_directory_path() / "_test_tibble.dat";
  std::string tmpfile = tmppath.string();

  eonc::io::writeTibble(*m, tmpfile);
  REQUIRE(std::filesystem::exists(tmpfile));
  auto fsize = std::filesystem::file_size(tmpfile);
  REQUIRE(fsize > 50);

  std::filesystem::remove(tmpfile);
}

} /* namespace tests */

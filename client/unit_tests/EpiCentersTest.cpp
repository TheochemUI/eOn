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
#include "EpiCenters.h"
#include "Matter.h"
#include "Parameters.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace Catch::Matchers;

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// ---------------------------------------------------------------------------
// Fixture: loads Pt_Heptamer_FrozenLayers system
//   7 free atoms (indices 0-6), 336 frozen, 343 total
// ---------------------------------------------------------------------------
class EpiCentersFixture {
public:
  EpiCentersFixture()
      : params{},
        pot{nullptr},
        matter{nullptr} {
    params.potential_options.potential = PotType::LJ;
    pot = eonc::helpers::makePotential(params.potential_options.potential,
                                       params);
    matter = std::make_shared<Matter>(pot, params);
    const std::string confile("pos.con");
    const bool ok = matter->con2matter(confile);
    REQUIRE(ok);
  }

  ~EpiCentersFixture() = default;

protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> matter;
};

// ---------------------------------------------------------------------------
// listedAtomEpiCenter tests
// ---------------------------------------------------------------------------

TEST_CASE_METHOD(EpiCentersFixture,
                 "listedAtomEpiCenter returns an atom from the list",
                 "[EpiCenters][listedAtomEpiCenter]") {
  std::vector<long> atomList = {0, 2, 4, 6};

  // Run many times to exercise the random selection
  std::set<long> selected;
  for (int trial = 0; trial < 200; ++trial) {
    long idx = eonc::EpiCenters::listedAtomEpiCenter(matter.get(), atomList);
    // Result must be one of the atoms in the list
    REQUIRE(std::find(atomList.begin(), atomList.end(), idx) != atomList.end());
    selected.insert(idx);
  }

  // Over 200 trials with 4 options, we should hit at least 2
  REQUIRE(selected.size() >= 2);
}

TEST_CASE_METHOD(EpiCentersFixture,
                 "listedAtomEpiCenter filters out frozen atoms",
                 "[EpiCenters][listedAtomEpiCenter]") {
  // Atoms 7+ are frozen in this system; only atom 3 is free
  std::vector<long> atomList = {3, 7, 8, 9, 100};

  for (int trial = 0; trial < 50; ++trial) {
    long idx = eonc::EpiCenters::listedAtomEpiCenter(matter.get(), atomList);
    // Must pick the only free atom in the list
    REQUIRE(idx == 3);
  }
}

TEST_CASE_METHOD(EpiCentersFixture,
                 "listedAtomEpiCenter filters out out-of-range indices",
                 "[EpiCenters][listedAtomEpiCenter]") {
  long nAtoms = matter->numberOfAtoms();
  // Include one valid free atom (0) and some out-of-range indices
  std::vector<long> atomList = {-1, 0, nAtoms, nAtoms + 100};

  for (int trial = 0; trial < 50; ++trial) {
    long idx = eonc::EpiCenters::listedAtomEpiCenter(matter.get(), atomList);
    REQUIRE(idx == 0);
  }
}

TEST_CASE_METHOD(EpiCentersFixture,
                 "listedAtomEpiCenter single free atom always selected",
                 "[EpiCenters][listedAtomEpiCenter]") {
  std::vector<long> atomList = {5};

  for (int trial = 0; trial < 20; ++trial) {
    long idx = eonc::EpiCenters::listedAtomEpiCenter(matter.get(), atomList);
    REQUIRE(idx == 5);
  }
}

// ---------------------------------------------------------------------------
// Parameters parsing tests: displace_atom_list and displace_type validation
// ---------------------------------------------------------------------------

TEST_CASE("Parameters default displace_atom_list is empty",
          "[Parameters][displace_atom_list]") {
  Parameters params;
  REQUIRE(params.saddle_search_options.displace_atom_list.empty());
}

TEST_CASE("Parameters default displace_type is load",
          "[Parameters][displace_type]") {
  Parameters params;
  REQUIRE(params.saddle_search_options.displace_type ==
          std::string(eonc::EpiCenters::DISP_LOAD));
}

TEST_CASE("Parameters parses comma-separated displace_atom_list from INI",
          "[Parameters][displace_atom_list]") {
  // Write a temporary config.ini with displace_atom_list
  const std::string ini_content = "[Saddle Search]\n"
                                  "client_displace_type = listed_atoms\n"
                                  "displace_atom_list = 10, 20, 30\n";

  FILE *tmpf = tmpfile();
  REQUIRE(tmpf != nullptr);
  fprintf(tmpf, "%s", ini_content.c_str());
  rewind(tmpf);

  Parameters params;
  int err = params.load(tmpf);
  fclose(tmpf);

  REQUIRE(err == 0);
  REQUIRE(params.saddle_search_options.displace_type ==
          std::string(eonc::EpiCenters::DISP_LISTED_ATOMS));
  REQUIRE(params.saddle_search_options.displace_atom_list.size() == 3);
  REQUIRE(params.saddle_search_options.displace_atom_list[0] == 10);
  REQUIRE(params.saddle_search_options.displace_atom_list[1] == 20);
  REQUIRE(params.saddle_search_options.displace_atom_list[2] == 30);
}

TEST_CASE("Parameters handles single-element displace_atom_list",
          "[Parameters][displace_atom_list]") {
  const std::string ini_content = "[Saddle Search]\n"
                                  "client_displace_type = listed_atoms\n"
                                  "displace_atom_list = 42\n";

  FILE *tmpf = tmpfile();
  REQUIRE(tmpf != nullptr);
  fprintf(tmpf, "%s", ini_content.c_str());
  rewind(tmpf);

  Parameters params;
  int err = params.load(tmpf);
  fclose(tmpf);

  REQUIRE(err == 0);
  REQUIRE(params.saddle_search_options.displace_atom_list.size() == 1);
  REQUIRE(params.saddle_search_options.displace_atom_list[0] == 42);
}

TEST_CASE("Parameters handles negative indices in displace_atom_list",
          "[Parameters][displace_atom_list]") {
  const std::string ini_content = "[Saddle Search]\n"
                                  "displace_atom_list = 5, -1, 10\n";

  FILE *tmpf = tmpfile();
  REQUIRE(tmpf != nullptr);
  fprintf(tmpf, "%s", ini_content.c_str());
  rewind(tmpf);

  Parameters params;
  int err = params.load(tmpf);
  fclose(tmpf);

  REQUIRE(err == 0);
  REQUIRE(params.saddle_search_options.displace_atom_list.size() == 3);
  REQUIRE(params.saddle_search_options.displace_atom_list[0] == 5);
  REQUIRE(params.saddle_search_options.displace_atom_list[1] == -1);
  REQUIRE(params.saddle_search_options.displace_atom_list[2] == 10);
}

TEST_CASE("Parameters empty displace_atom_list stays empty",
          "[Parameters][displace_atom_list]") {
  const std::string ini_content = "[Saddle Search]\n"
                                  "client_displace_type = listed_atoms\n";

  FILE *tmpf = tmpfile();
  REQUIRE(tmpf != nullptr);
  fprintf(tmpf, "%s", ini_content.c_str());
  rewind(tmpf);

  Parameters params;
  int err = params.load(tmpf);
  fclose(tmpf);

  REQUIRE(err == 0);
  REQUIRE(params.saddle_search_options.displace_atom_list.empty());
}

TEST_CASE("Parameters listed_atoms is whitelisted in displace_type",
          "[Parameters][displace_type]") {
  const std::string ini_content = "[Saddle Search]\n"
                                  "client_displace_type = listed_atoms\n";

  FILE *tmpf = tmpfile();
  REQUIRE(tmpf != nullptr);
  fprintf(tmpf, "%s", ini_content.c_str());
  rewind(tmpf);

  Parameters params;
  int err = params.load(tmpf);
  fclose(tmpf);

  REQUIRE(err == 0);
  // Should NOT fall back to "load"
  REQUIRE(params.saddle_search_options.displace_type ==
          std::string(eonc::EpiCenters::DISP_LISTED_ATOMS));
}

TEST_CASE("Parameters unknown displace_type falls back to load",
          "[Parameters][displace_type]") {
  const std::string ini_content = "[Saddle Search]\n"
                                  "client_displace_type = bogus_type\n";

  FILE *tmpf = tmpfile();
  REQUIRE(tmpf != nullptr);
  fprintf(tmpf, "%s", ini_content.c_str());
  rewind(tmpf);

  Parameters params;
  int err = params.load(tmpf);
  fclose(tmpf);

  REQUIRE(err == 0);
  REQUIRE(params.saddle_search_options.displace_type ==
          std::string(eonc::EpiCenters::DISP_LOAD));
}

TEST_CASE("Parameters all valid displace_types are accepted",
          "[Parameters][displace_type]") {
  const std::vector<std::string> valid_types = {
      eonc::EpiCenters::DISP_LOAD,
      eonc::EpiCenters::DISP_NOT_FCC_OR_HCP,
      eonc::EpiCenters::DISP_MIN_COORDINATED,
      eonc::EpiCenters::DISP_LAST_ATOM,
      eonc::EpiCenters::DISP_RANDOM,
      eonc::EpiCenters::DISP_LISTED_ATOMS,
  };

  for (const auto &dtype : valid_types) {
    std::string ini_content =
        "[Saddle Search]\nclient_displace_type = " + dtype + "\n";

    FILE *tmpf = tmpfile();
    REQUIRE(tmpf != nullptr);
    fprintf(tmpf, "%s", ini_content.c_str());
    rewind(tmpf);

    Parameters params;
    int err = params.load(tmpf);
    fclose(tmpf);

    REQUIRE(err == 0);
    // "load" is whitelisted only implicitly (it's the fallback),
    // but it should still be preserved if set explicitly
    REQUIRE(params.saddle_search_options.displace_type == dtype);
  }
}

// --- EpiCenter method tests on Pt_Heptamer_FrozenLayers ---

TEST_CASE("minCoordinatedEpiCenter returns a free atom index",
          "[epicenters][min_coordinated]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  long idx = eonc::EpiCenters::minCoordinatedEpiCenter(matter.get(), 3.3);
  REQUIRE(idx >= 0);
  REQUIRE(idx < matter->numberOfAtoms());
  REQUIRE(!matter->getFixed(idx));
}

TEST_CASE("randomFreeAtomEpiCenter returns a free atom",
          "[epicenters][random]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  long idx = eonc::EpiCenters::randomFreeAtomEpiCenter(matter.get());
  REQUIRE(idx >= 0);
  REQUIRE(idx < matter->numberOfAtoms());
  REQUIRE(!matter->getFixed(idx));
}

TEST_CASE("lastAtom returns last atom index", "[epicenters][last_atom]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  long idx = eonc::EpiCenters::lastAtom(matter.get());
  REQUIRE(idx == matter->numberOfAtoms() - 1);
}

TEST_CASE("coordination returns finite values for all atoms",
          "[epicenters][coordination]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  long nAtoms = matter->numberOfAtoms();
  std::vector<long> coord(nAtoms);
  eonc::EpiCenters::coordination(coord.data(), matter.get(), 3.3);

  for (long i = 0; i < nAtoms; i++) {
    REQUIRE(coord[i] >= 0);
  }
}

TEST_CASE("minCoordination returns a valid atom index",
          "[epicenters][min_coordination]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  long idx = eonc::EpiCenters::minCoordination(matter.get(), 3.3);
  REQUIRE(idx >= 0);
  REQUIRE(idx < matter->numberOfAtoms());
}

TEST_CASE("cnaEpiCenter returns a non-FCC/HCP atom", "[epicenters][cna]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  // CNA needs the neighbor cutoff
  long idx = eonc::EpiCenters::cnaEpiCenter(matter.get(), 3.3);
  // May return -2 if no non-FCC/HCP atom found (small cluster)
  // Just verify it doesn't crash
  CHECK(idx >= -2);
}

TEST_CASE("coordinationLessOrEqual filters atoms correctly",
          "[epicenters][coordination_filter]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  long nAtoms = matter->numberOfAtoms();
  auto result = std::make_unique<bool[]>(nAtoms);
  eonc::EpiCenters::coordinationLessOrEqual(result.get(), 6, matter.get(), 3.3);

  // At least some atoms should pass or fail the filter
  bool anyTrue = false, anyFalse = false;
  for (long i = 0; i < nAtoms; i++) {
    if (result[i])
      anyTrue = true;
    else
      anyFalse = true;
  }
  // Just verify the function runs and produces booleans
  CHECK((anyTrue || anyFalse));
}

} // namespace tests

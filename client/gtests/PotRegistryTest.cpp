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

#include "PotRegistry.h"
#include "Potential.h"
#include "catch2/catch_amalgamated.hpp"
#include <fstream>
#include <memory>

using namespace eonc;

// Minimal concrete potential for testing
class DummyPotential : public Potential {
public:
  DummyPotential(PotType pt, const Parameters &p) : Potential(pt, p) {}

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    *energy = 0.0;
    for (long i = 0; i < nAtoms * 3; ++i)
      forces[i] = 0.0;
  }
};

TEST_CASE("PotRegistry tracks creation and destruction",
          "[PotRegistry][lifecycle]") {
  PotRegistry::get().reset();

  Parameters params;
  params.potential_options.potential = PotType::LJ;

  SECTION("on_created increments alive and created counts") {
    REQUIRE(PotRegistry::get().type_alive(PotType::LJ) == 0);
    {
      DummyPotential pot(PotType::LJ, params);
      REQUIRE(PotRegistry::get().type_alive(PotType::LJ) == 1);
    }
    REQUIRE(PotRegistry::get().type_alive(PotType::LJ) == 0);
  }

  SECTION("multiple instances tracked independently") {
    {
      DummyPotential pot1(PotType::LJ, params);
      DummyPotential pot2(PotType::LJ, params);
      DummyPotential pot3(PotType::EMT, params);
      REQUIRE(PotRegistry::get().type_alive(PotType::LJ) == 2);
      REQUIRE(PotRegistry::get().type_alive(PotType::EMT) == 1);
    }
    REQUIRE(PotRegistry::get().type_alive(PotType::LJ) == 0);
    REQUIRE(PotRegistry::get().type_alive(PotType::EMT) == 0);
  }
}

TEST_CASE("PotRegistry tracks force calls", "[PotRegistry][force_calls]") {
  PotRegistry::get().reset();

  Parameters params;
  params.potential_options.potential = PotType::LJ;

  SECTION("on_force_call increments per-type and total counters") {
    DummyPotential pot(PotType::LJ, params);
    REQUIRE(PotRegistry::get().type_force_calls(PotType::LJ) == 0);
    REQUIRE(PotRegistry::get().total_force_calls() == 0);

    PotRegistry::get().on_force_call(PotType::LJ);
    PotRegistry::get().on_force_call(PotType::LJ);
    PotRegistry::get().on_force_call(PotType::LJ);

    REQUIRE(PotRegistry::get().type_force_calls(PotType::LJ) == 3);
    REQUIRE(PotRegistry::get().total_force_calls() == 3);
  }

  SECTION("get_ef increments both forceCallCounter and registry") {
    auto pot = std::make_shared<DummyPotential>(PotType::LJ, params);

    // Two atoms, simple setup
    AtomMatrix pos(2, 3);
    pos << 0.0, 0.0, 0.0, 1.5, 0.0, 0.0;
    VectorXi atmnrs(2);
    atmnrs << 79, 79;
    Matrix3d box = Matrix3d::Identity() * 10.0;

    auto [e1, f1] = pot->get_ef(pos, atmnrs, box);
    auto [e2, f2] = pot->get_ef(pos, atmnrs, box);

    REQUIRE(pot->forceCallCounter == 2);
    REQUIRE(PotRegistry::get().type_force_calls(PotType::LJ) == 2);
    REQUIRE(PotRegistry::get().total_force_calls() == 2);
  }

  SECTION("delta tracking pattern works") {
    auto pot = std::make_shared<DummyPotential>(PotType::LJ, params);
    AtomMatrix pos(2, 3);
    pos << 0.0, 0.0, 0.0, 1.5, 0.0, 0.0;
    VectorXi atmnrs(2);
    atmnrs << 79, 79;
    Matrix3d box = Matrix3d::Identity() * 10.0;

    // Simulate the delta tracking pattern used by jobs
    size_t baseline = PotRegistry::get().total_force_calls();
    pot->get_ef(pos, atmnrs, box);
    pot->get_ef(pos, atmnrs, box);
    pot->get_ef(pos, atmnrs, box);
    size_t delta = PotRegistry::get().total_force_calls() - baseline;

    REQUIRE(delta == 3);

    // Second phase
    baseline = PotRegistry::get().total_force_calls();
    pot->get_ef(pos, atmnrs, box);
    pot->get_ef(pos, atmnrs, box);
    delta = PotRegistry::get().total_force_calls() - baseline;

    REQUIRE(delta == 2);
    REQUIRE(PotRegistry::get().total_force_calls() == 5);
  }
}

TEST_CASE("PotRegistry reset clears all state", "[PotRegistry][reset]") {
  PotRegistry::get().reset();

  Parameters params;
  params.potential_options.potential = PotType::EMT;

  // Create some state
  PotRegistry::get().on_force_call(PotType::EMT);
  PotRegistry::get().on_force_call(PotType::EMT);
  REQUIRE(PotRegistry::get().total_force_calls() == 2);

  // Create and destroy a potential to generate a record
  { DummyPotential pot(PotType::EMT, params); }

  PotRegistry::get().reset();

  REQUIRE(PotRegistry::get().total_force_calls() == 0);
  REQUIRE(PotRegistry::get().type_force_calls(PotType::EMT) == 0);
  REQUIRE(PotRegistry::get().type_alive(PotType::EMT) == 0);
}

TEST_CASE("PotRegistry write_summary produces valid JSON",
          "[PotRegistry][json]") {
  PotRegistry::get().reset();

  Parameters params;
  params.potential_options.potential = PotType::LJ;

  // Create and destroy a potential, calling get_ef to register force calls
  {
    DummyPotential pot(PotType::LJ, params);
    AtomMatrix pos(2, 3);
    pos << 0.0, 0.0, 0.0, 1.5, 0.0, 0.0;
    VectorXi atmnrs(2);
    atmnrs << 79, 79;
    Matrix3d box = Matrix3d::Identity() * 10.0;
    pot.get_ef(pos, atmnrs, box);
    pot.get_ef(pos, atmnrs, box);
  }

  std::string path = "_test_potcalls.json";
  PotRegistry::get().write_summary(path);

  // Read and verify basic structure
  std::ifstream ifs(path);
  REQUIRE(ifs.is_open());
  std::string content((std::istreambuf_iterator<char>(ifs)),
                      std::istreambuf_iterator<char>());
  ifs.close();
  std::remove(path.c_str());

  // Verify JSON array structure
  REQUIRE(content.front() == '[');
  REQUIRE(content.find("\"id\": 1") != std::string::npos);
  REQUIRE(content.find("\"type\": \"LJ\"") != std::string::npos);
  REQUIRE(content.find("\"force_calls\": 2") != std::string::npos);
  REQUIRE(content.find("\"created_at\"") != std::string::npos);
  REQUIRE(content.find("\"destroyed_at\"") != std::string::npos);
}

TEST_CASE("PotRegistry handles multiple types correctly",
          "[PotRegistry][multi_type]") {
  PotRegistry::get().reset();

  Parameters params;

  {
    DummyPotential lj(PotType::LJ, params);
    DummyPotential emt(PotType::EMT, params);

    PotRegistry::get().on_force_call(PotType::LJ);
    PotRegistry::get().on_force_call(PotType::LJ);
    PotRegistry::get().on_force_call(PotType::EMT);

    REQUIRE(PotRegistry::get().type_force_calls(PotType::LJ) == 2);
    REQUIRE(PotRegistry::get().type_force_calls(PotType::EMT) == 1);
    REQUIRE(PotRegistry::get().total_force_calls() == 3);
    REQUIRE(PotRegistry::get().type_alive(PotType::LJ) == 1);
    REQUIRE(PotRegistry::get().type_alive(PotType::EMT) == 1);
  }

  REQUIRE(PotRegistry::get().type_alive(PotType::LJ) == 0);
  REQUIRE(PotRegistry::get().type_alive(PotType::EMT) == 0);
  REQUIRE(PotRegistry::get().total_force_calls() == 3);
}

TEST_CASE("PotRegistry instance records capture correct force_calls",
          "[PotRegistry][records]") {
  PotRegistry::get().reset();

  Parameters params;
  params.potential_options.potential = PotType::LJ;

  // Instance with 5 force calls
  {
    DummyPotential pot(PotType::LJ, params);
    AtomMatrix pos(2, 3);
    pos << 0.0, 0.0, 0.0, 1.5, 0.0, 0.0;
    VectorXi atmnrs(2);
    atmnrs << 79, 79;
    Matrix3d box = Matrix3d::Identity() * 10.0;

    for (int i = 0; i < 5; ++i)
      pot.get_ef(pos, atmnrs, box);
  }

  // Instance with 0 force calls (diagnostic signal)
  { DummyPotential pot(PotType::LJ, params); }

  // Write and verify both records appear
  std::string path = "_test_potcalls_records.json";
  PotRegistry::get().write_summary(path);

  std::ifstream ifs(path);
  std::string content((std::istreambuf_iterator<char>(ifs)),
                      std::istreambuf_iterator<char>());
  ifs.close();
  std::remove(path.c_str());

  REQUIRE(content.find("\"force_calls\": 5") != std::string::npos);
  REQUIRE(content.find("\"force_calls\": 0") != std::string::npos);
}

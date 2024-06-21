#include "Matter.h"
#include "Parameters.h"
#include "catch2/catch_amalgamated.hpp"
#include <memory>

using namespace Catch::Matchers;

TEST_CASE("TestCell", "[MatterTest]") {
  auto params = std::make_shared<Parameters>();
  auto pot_default = helper_functions::makePotential(PotType::LJ, params);
  auto m1 = std::make_shared<Matter>(pot_default, params);
  std::string confile("pos.con");
  m1->con2matter(confile);

  Eigen::Matrix3d _cell;
  Eigen::Matrix3d _cellInverse;
  // clang-format off
    _cell << // Comma initialized
        25.0, 0.0, 0.0,
        0.0, 25.0, 0.0,
        0.0, 0.0, 25.0;
  // clang-format on

  REQUIRE_THAT(m1->getCell()(0, 0), WithinAbs(_cell(0, 0), 0.01));
  REQUIRE_THAT(m1->getCell()(1, 1), WithinAbs(_cell(1, 1), 0.01));
  REQUIRE_THAT(m1->getCell()(2, 2), WithinAbs(_cell(2, 2), 0.01));
}

TEST_CASE("SetGetAtomicNrs", "[MatterTest]") {
  auto params = std::make_shared<Parameters>();
  auto pot_default = helper_functions::makePotential(PotType::LJ, params);
  auto m1 = std::make_shared<Matter>(pot_default, params);
  std::string confile("pos.con");
  m1->con2matter(confile);

  Eigen::VectorXi _atmnrs(13);
  _atmnrs << 8, 8, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 16;
  Eigen::VectorXi _atmnrs2(13);
  _atmnrs2 << 16, 16, 12, 12, 12, 12, 2, 2, 2, 2, 2, 2, 32;

  REQUIRE(m1->getAtomicNrs() == _atmnrs);

  for (auto &atmnr : _atmnrs) {
    atmnr *= 2;
  }
  m1->setAtomicNrs(_atmnrs);

  REQUIRE(m1->getAtomicNrs() == _atmnrs2);
}

TEST_CASE("SetPotential", "[MatterTest]") {
  auto params = std::make_shared<Parameters>();
  auto pot_default = helper_functions::makePotential(PotType::LJ, params);
  auto m1 = std::make_shared<Matter>(pot_default, params);
  std::string confile("pos.con");
  m1->con2matter(confile);

  double m1_ipot = m1->getPotentialEnergy();
  params->potential = PotType::MORSE_PT;
  auto pot = helper_functions::makePotential(params->potential, params);

  REQUIRE(m1->getPotential() != pot);
  m1->setPotential(pot);

  double m1_fpot = m1->getPotentialEnergy();
  REQUIRE_THAT(m1_ipot, WithinAbs(-8.9245813315, 0.01));
  REQUIRE_THAT(m1_fpot, WithinAbs(1611.8672392832, 0.01));
  REQUIRE(m1_ipot != m1_fpot);
}

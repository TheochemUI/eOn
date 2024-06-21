#include "../MatrixHelpers.hpp"
#include "Matter.h"
#include "catch2/catch_amalgamated.hpp"

using namespace std::placeholders;
using namespace Catch::Matchers;

namespace tests {

class PotTest {
public:
  PotTest()
      : params{std::make_shared<Parameters>()},
        m1{nullptr},
        pot_default{nullptr},
        threshold{1e-2} {}

  ~PotTest() {}

  void SetUp() {
    pot_default = helper_functions::makePotential(PotType::LJ, params);
    m1 = std::make_shared<Matter>(pot_default, params);
    std::string confile("pos.con");
    m1->con2matter(confile);
  }

  void TearDown() {}

protected:
  std::shared_ptr<Parameters> params;
  std::shared_ptr<Matter> m1;
  std::shared_ptr<Potential> pot_default;
  double threshold;
};

TEST_CASE_METHOD(PotTest, "getType", "[PotTest]") {
  SetUp();
  auto params = std::make_shared<Parameters>();
  params->potential = PotType::LJ;
  std::shared_ptr<Potential> pot =
      helper_functions::makePotential(PotType::LJ, params);
  REQUIRE(pot->getType() == PotType::LJ);
  params->potential = PotType::MORSE_PT;
  std::shared_ptr<Potential> pot2 =
      helper_functions::makePotential(PotType::MORSE_PT, params);
  REQUIRE(pot2->getType() == PotType::MORSE_PT);
  TearDown();
}

TEST_CASE_METHOD(PotTest, "callForce", "[PotTest]") {
  SetUp();
  auto matEq =
      std::bind(helper_functions::eigenEquality<AtomMatrix>, _1, _2, threshold);

  double energy_lj{-8.9245539406};
  AtomMatrix forces_lj{Eigen::MatrixXd::Zero(m1->numberOfAtoms(), 3)};
  forces_lj << 0.463529, -0.24841, -1.30892, 1.08313, 0.0575693, 0.886135,
      -2.20116, -2.23402, -2.09064, 2.23797, -2.40278, 2.50396, -2.71153,
      0.972603, -2.41778, 2.14953, 3.62404, 2.3242, 0.823661, 1.02942, 1.04033,
      -0.488582, 1.86173, 0.036945, -0.234744, -0.460397, -1.21061, -0.0650386,
      -2.24842, 0.560231, 1.03832, -0.512915, 1.15749, -0.813122, 0.310969,
      -1.78941, -1.28196, 0.250611, 0.30807;

  double energy_morse{1611.8700112470};
  AtomMatrix forces_morse{Eigen::MatrixXd::Zero(m1->numberOfAtoms(), 3)};
  forces_morse << -64.9256, 38.3478, 194.53, -157.289, -7.61147, -134.612,
      -164.748, -274.908, -285.305, 76.0495, -338.668, 421.113, -226.069,
      166.355, -318.858, 106.552, 385.35, 298.415, 364.555, 448.297, 440.502,
      -201.26, 742.277, 84.4901, -106.268, -272.671, -696.314, -12.0385,
      -743.517, 239.385, 470, -229.572, 504.702, -293.633, 125.068, -697.975,
      209.073, -38.7473, -50.0734;

  params->potential = PotType::LJ;
  std::shared_ptr<Potential> pot =
      helper_functions::makePotential(params->potential, params);
  double e_lj{0};
  AtomMatrix f_lj = Eigen::MatrixXd::Ones(m1->numberOfAtoms(), 3);

  pot->force(m1->numberOfAtoms(), m1->getPositions().data(),
             m1->getAtomicNrs().data(), f_lj.data(), &e_lj, nullptr,
             m1->getCell().data());
  REQUIRE_THAT(e_lj, WithinAbs(energy_lj, threshold));
  REQUIRE(matEq(forces_lj, f_lj));

  double e_morse{0};
  AtomMatrix f_morse = Eigen::MatrixXd::Ones(m1->numberOfAtoms(), 3);
  params->potential = PotType::MORSE_PT;
  std::shared_ptr<Potential> pot2 =
      helper_functions::makePotential(params->potential, params);
  REQUIRE(pot2 != pot);

  pot2->force(m1->numberOfAtoms(), m1->getPositions().data(),
              m1->getAtomicNrs().data(), f_morse.data(), &e_morse, nullptr,
              m1->getCell().data());
  REQUIRE_THAT(e_morse, WithinAbs(energy_morse, threshold));
  REQUIRE(matEq(forces_morse, f_morse));

  TearDown();
}

} /* namespace tests */

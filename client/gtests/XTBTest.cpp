#include "../MatrixHelpers.hpp"
#include "Matter.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

using namespace std::placeholders;
using namespace Catch::Matchers;

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class PotTest {
public:
  PotTest()
      : params{},
        m1{nullptr},
        pot_default{nullptr},
        threshold{1e-2} {}

  ~PotTest() {}

  void SetUp() {
    pot_default = eonc::helpers::makePotential(PotType::LJ, params);
    m1 = std::make_shared<Matter>(pot_default, params);
    std::string confile("pos.con");
    m1->con2matter(confile);
  }

  void TearDown() {}

protected:
  Parameters params;
  std::shared_ptr<Matter> m1;
  std::shared_ptr<Potential> pot_default;
  double threshold;
};

TEST_CASE_METHOD(PotTest, "XTB", "[PotTest]") {
  SetUp();
  auto matEq =
      std::bind(eonc::helpers::eigenEquality<AtomMatrix>, _1, _2, threshold);
  double expected_energy = -623.58142247693342597;
  AtomMatrix expected_forces(m1->numberOfAtoms(), 3);
  // clang-format off
  expected_forces <<
     1.22362,   -0.111937,  -0.65068,
     1.41682,   -0.0333704, -0.0175255,
     0.191715,  -0.169379,   0.123635,
    -0.0569636,  0.58173,   -0.475738,
    -0.0185091,  0.273378,   0.105843,
    -0.41846,   -0.276634,  -0.174549,
    -0.396333,  -0.344539,  -0.21208,
     0.145317,  -0.274345,   0.0193664,
    -0.0648892, -0.0536216,  0.597702,
     0.134743,   0.290109,  -0.0652169,
    -0.449305,   0.144147,  -0.289649,
     0.121169,   0.112835,   0.485777,
    -1.82893,   -0.138373,   0.553116;
  // clang-format on

  double e_mta{0};
  AtomMatrix f_mta = MatrixXd::Ones(m1->numberOfAtoms(), 3);
  params.potential_options.potential = PotType::XTB;
  params.xtb_options.paramset = "GFN2xTB";
  params.xtb_options.acc = 1.0;
  params.xtb_options.elec_temperature = 300.0;
  params.xtb_options.maxiter = 250;
  params.xtb_options.charge = 0.0;
  params.xtb_options.uhf = 0;
  auto pot = eonc::helpers::makePotential(params.potential_options.potential,
                                             params);
  pot->force(m1->numberOfAtoms(), m1->getPositions().data(),
             m1->getAtomicNrs().data(), f_mta.data(), &e_mta, nullptr,
             m1->getCell().data());
  SECTION("Energy Check") {
    REQUIRE_THAT(e_mta, WithinAbs(expected_energy, threshold));
  }
  SECTION("Force Matrix Check") {
    REQUIRE_THAT(f_mta,
                 eonc::helpers::test::IsApprox(expected_forces, threshold));
  }
  // Call again to see that update works
  pot->force(m1->numberOfAtoms(), m1->getPositions().data(),
             m1->getAtomicNrs().data(), f_mta.data(), &e_mta, nullptr,
             m1->getCell().data());
  SECTION("Energy Check") {
    REQUIRE_THAT(e_mta, WithinAbs(expected_energy, threshold));
  }
  SECTION("Force Matrix Check") {
    REQUIRE_THAT(f_mta,
                 eonc::helpers::test::IsApprox(expected_forces, threshold));
  }
  TearDown();
}

} /* namespace tests */

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

TEST_CASE_METHOD(PotTest, "Metatomic", "[PotTest]") {
  SetUp();
  auto matEq =
      std::bind(helper_functions::eigenEquality<AtomMatrix>, _1, _2, threshold);
  double expected_energy = 98374.87753058573;
  AtomMatrix expected_forces(m1->numberOfAtoms(), 3);
  expected_forces << -90017.67874663, 14048.6323898, 42667.56998833,
      51715.62564964, 81020.09231365, -33288.72879168, 13893.5617694,
      89793.76258628, 33937.11367155, -71755.76235155, 52324.79462503,
      -32402.45955531, 46268.61994191, -89837.93451292, 11509.42733453,
      -67770.6981497, -6678.90441423, -33619.79643378, -17329.67691782,
      -68054.08864133, -59929.39232182, -14176.97088206, 30833.36510504,
      -85864.44933633, -75239.50566824, -57104.7492775, -3715.98865599,
      90169.5778592, 5214.90295971, 30786.57650636, 20281.84422414,
      21840.57389606, 85752.46910389, 104913.11853287, -11120.68664665,
      -29664.02315515, 9047.94473885, -62279.76038293, 73831.6816454;

  double e_mta{0};
  AtomMatrix f_mta = Eigen::MatrixXd::Ones(m1->numberOfAtoms(), 3);
  params->potential = PotType::METATOMIC;
  params->metatomic_options.model_path = "lennard-jones.pt";
  auto pot = helper_functions::makePotential(params->potential, params);
  pot->force(m1->numberOfAtoms(), m1->getPositions().data(),
             m1->getAtomicNrs().data(), f_mta.data(), &e_mta, nullptr,
             m1->getCell().data());
  REQUIRE_THAT(e_mta, WithinAbs(expected_energy, threshold));
  REQUIRE(matEq(f_mta, expected_forces));
  TearDown();
}

} /* namespace tests */

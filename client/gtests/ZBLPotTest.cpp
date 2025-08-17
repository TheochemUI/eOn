#include "../MatrixHelpers.hpp"
#include "Matter.h"
#include "catch2/catch_amalgamated.hpp"
#include <memory>
#include <string>

using namespace std::placeholders;
using namespace Catch::Matchers;

namespace tests {

class ZBLPotTest {
public:
  ZBLPotTest()
      : params{std::make_shared<Parameters>()},
        matter{nullptr},
        pot_zbl{nullptr},
        threshold{1e-6}
  {
    params->potential = PotType::ZBL;
    params->zbl_options.cut_inner = 2.0;
    params->zbl_options.cut_global = 2.5;

    pot_zbl = helper_functions::makePotential(params->potential, params);
    matter = std::make_shared<Matter>(pot_zbl, params);

    const std::string confile("pos.con");
    const bool file_read_ok = matter->con2matter(confile);
    REQUIRE(file_read_ok);
  }

  ~ZBLPotTest() = default;

protected:
  std::shared_ptr<Parameters> params;
  std::shared_ptr<Matter> matter;
  std::shared_ptr<Potential> pot_zbl;
  double threshold;
};

TEST_CASE_METHOD(ZBLPotTest, "ZBL Potential against LAMMPS", "[PotTest][ZBL]") {
  // Expected reference values from LAMMPS
  // See the ZBL/wip_exploration for how these are generated
  const double expected_energy = 0.38537731;
  AtomMatrix expected_forces(2, 3);
  expected_forces.row(0) << -2.37926, -2.57753, -2.7758; // Si
  expected_forces.row(1) << 2.37926, 2.57753, 2.7758;    // Au

  double calculated_energy = 0.0;
  AtomMatrix calculated_forces =
      Eigen::MatrixXd::Zero(matter->numberOfAtoms(), 3);

  pot_zbl->force(matter->numberOfAtoms(), matter->getPositions().data(),
                 matter->getAtomicNrs().data(), calculated_forces.data(),
                 &calculated_energy,
                 nullptr, // no variance
                 matter->getCell().data());

  REQUIRE_THAT(calculated_energy, WithinAbs(expected_energy, threshold));

  auto matEq =
      std::bind(helper_functions::eigenEquality<AtomMatrix>, _1, _2, threshold);
  REQUIRE(matEq(calculated_forces, expected_forces));
}

} // namespace tests

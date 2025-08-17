#include "../MatrixHelpers.hpp"
#include "Matter.h"
#include "catch2/catch_amalgamated.hpp"
#include <memory>
#include <string>

using namespace std::placeholders;
using namespace Catch::Matchers;

namespace tests {

class SocketNWChemPotTest {
public:
  SocketNWChemPotTest()
      : params{std::make_shared<Parameters>()},
        matter{nullptr},
        pot_socket{nullptr},
        threshold{1e-5}
  {
    params->potential = PotType::SocketNWChem;
    params->socket_nwchem_options.unix_socket_mode = true;
    params->socket_nwchem_options.unix_socket_path = "eon_nwchem_test_socket";

    pot_socket = helper_functions::makePotential(params->potential, params);
    matter = std::make_shared<Matter>(pot_socket, params);

    const std::string confile("pos.con");
    const bool file_read_ok = matter->con2matter(confile);
    REQUIRE(file_read_ok);
    REQUIRE(matter->numberOfAtoms() == 9);
  }

  ~SocketNWChemPotTest() = default;

protected:
  std::shared_ptr<Parameters> params;
  std::shared_ptr<Matter> matter;
  std::shared_ptr<Potential> pot_socket;
  double threshold;
};

TEST_CASE_METHOD(SocketNWChemPotTest, "SocketNWChemPot against reference data",
                 "[PotTest][SocketNWChem]") {
  const double expected_energy = -7078.99177203005;
  AtomMatrix expected_forces(9, 3);
  expected_forces << -3.35282162e+00, 2.62823327e+00, -1.69692821e-03,
      1.02335056e+00, 7.20474582e-01, -5.14220671e-05, 1.49756486e+00,
      -3.49562070e+00, 1.23412961e-03, -9.80361709e-01, -1.12871437e-01,
      1.13128548e-03, 1.12028115e+00, -4.64032733e-01, -1.23412961e-03,
      6.30074588e-01, -9.88846350e-02, -5.48159235e-01, 6.29766056e-01,
      -9.81647261e-02, 5.47336482e-01, -2.84106921e-01, 4.60073234e-01,
      -5.10672548e-01, -2.83695544e-01, 4.60793143e-01, 5.12060944e-01;

  double calculated_energy = 0.0;
  AtomMatrix calculated_forces =
      Eigen::MatrixXd::Zero(matter->numberOfAtoms(), 3);

  pot_socket->force(matter->numberOfAtoms(), matter->getPositions().data(),
                    matter->getAtomicNrs().data(), calculated_forces.data(),
                    &calculated_energy, nullptr, matter->getCell().data());

  REQUIRE_THAT(calculated_energy, WithinAbs(expected_energy, threshold));

  auto matEq =
      std::bind(helper_functions::eigenEquality<AtomMatrix>, _1, _2, threshold);
  REQUIRE(matEq(calculated_forces, expected_forces));
}

} // namespace tests

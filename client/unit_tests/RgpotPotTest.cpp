#include "../MatrixHelpers.hpp"
#include "Matter.h"
#include "catch2/catch_amalgamated.hpp"
#include <cmath>
#include <cstdlib>
#include <memory>

using namespace Catch::Matchers;

namespace tests {

namespace {

bool env_nonempty(const char *name) {
  const char *v = std::getenv(name);
  return v != nullptr && v[0] != '\0';
}

} // namespace

TEST_CASE("RgpotPot in-process nwchemc force (no potserv)",
          "[PotTest][RGPOT][nwchemc]") {
#ifndef WITH_RGPOT
  SKIP("built without WITH_RGPOT");
#else
  REQUIRE(env_nonempty("NWCHEMC_LIBRARY") ||
          env_nonempty("RGPOT_NWCHEMC_ENGINE") ||
          env_nonempty("RGPOT_NWCHEM_ENGINE"));

  Parameters params{};
  params.potential_options.potential = PotType::RGPOT;
  params.rgpot_options.backend = "nwchemc";
  params.rgpot_options.basis = "sto-3g";
  params.rgpot_options.theory = "scf";
  params.rgpot_options.scf_type = "rhf";
  params.rgpot_options.charge = 0;
  params.rgpot_options.multiplicity = 1;

  auto pot =
      eonc::helpers::makePotential(params.potential_options.potential, params);
  REQUIRE(pot != nullptr);
  REQUIRE(pot->getType() == PotType::RGPOT);

  auto matter = std::make_shared<Matter>(pot, params);
  REQUIRE(matter->con2matter("pos.con"));
  REQUIRE(matter->numberOfAtoms() == 9);

  double energy = 0.0;
  AtomMatrix forces = MatrixXd::Zero(matter->numberOfAtoms(), 3);
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy, nullptr,
             matter->getCell().data());

  REQUIRE(std::isfinite(energy));
  REQUIRE(std::abs(energy) > 1e-6);

  // Two successive forces must both succeed (warm multi-call path)
  double energy2 = 0.0;
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy2, nullptr,
             matter->getCell().data());
  REQUIRE(std::isfinite(energy2));
  REQUIRE(std::abs(energy - energy2) < 1e-4);
#endif
}

} // namespace tests

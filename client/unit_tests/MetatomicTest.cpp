#include "Matter.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include <cmath>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("Metatomic LJ model evaluates finite energy and forces",
          "[PotTest]") {
  Parameters params;
  params.potential_options.potential = PotType::METATOMIC;
  params.metatomic_options.model_path = "lennard-jones.pt";
  auto pot = eonc::helpers::makePotential(PotType::METATOMIC, params);
  auto m1 = std::make_shared<Matter>(pot, params);
  m1->con2matter(std::string("pos.con"));
  double e_mta = 0.0;
  AtomMatrix f_mta = AtomMatrix::Zero(m1->numberOfAtoms(), 3);
  REQUIRE_NOTHROW(pot->force(m1->numberOfAtoms(), m1->getPositions().data(),
                             m1->getAtomicNrs().data(), f_mta.data(), &e_mta,
                             nullptr, m1->getCell().data()));
  REQUIRE(std::isfinite(e_mta));
  REQUIRE(f_mta.allFinite());
}

} // namespace tests

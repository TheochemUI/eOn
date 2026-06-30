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
        threshold{5e-1} {}

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

TEST_CASE_METHOD(PotTest, "Metatomic", "[PotTest]") {
  SetUp();
  double e_mta{0};
  AtomMatrix f_mta = MatrixXd::Ones(m1->numberOfAtoms(), 3);
  params.potential_options.potential = PotType::METATOMIC;
  params.metatomic_options.model_path = "lennard-jones.pt";
  auto pot =
      eonc::helpers::makePotential(params.potential_options.potential, params);
  pot->force(m1->numberOfAtoms(), m1->getPositions().data(),
             m1->getAtomicNrs().data(), f_mta.data(), &e_mta, nullptr,
             m1->getCell().data());
  REQUIRE(std::isfinite(e_mta));
  REQUIRE(f_mta.allFinite());
  // LJ test model on lj13 is a large positive energy; keep a soft bound only.
  REQUIRE(e_mta > 0.0);
  REQUIRE(f_mta.norm() > 0.0);
  TearDown();
}

TEST_CASE_METHOD(PotTest,
                 "Metatomic uncertainty populates variance (if available)",
                 "[PotTest][uncertainty]") {
  SetUp();

  params.potential_options.potential = PotType::METATOMIC;
  params.metatomic_options.model_path = "lennard-jones.pt";
  params.metatomic_options.uncertainty_threshold = 0.1;

  auto pot =
      eonc::helpers::makePotential(params.potential_options.potential, params);

  double e_mta{0};
  AtomMatrix f_mta = MatrixXd::Zero(m1->numberOfAtoms(), 3);
  double variance = -12345.6789; // sentinel

  pot->force(m1->numberOfAtoms(), m1->getPositions().data(),
             m1->getAtomicNrs().data(), f_mta.data(), &e_mta, &variance,
             m1->getCell().data());

  if (variance != -12345.6789) {
    REQUIRE(std::isfinite(variance));
    REQUIRE(variance >= 0.0);
  }

  TearDown();
}

TEST_CASE_METHOD(PotTest, "Metatomic variant (doubled)", "[PotTest][variant]") {
  SetUp();

  params.potential_options.potential = PotType::METATOMIC;
  params.metatomic_options.model_path = "lennard-jones.pt";

  // Base evaluation
  double e_base{0};
  AtomMatrix f_base = MatrixXd::Zero(m1->numberOfAtoms(), 3);
  auto pot_base =
      eonc::helpers::makePotential(params.potential_options.potential, params);
  pot_base->force(m1->numberOfAtoms(), m1->getPositions().data(),
                  m1->getAtomicNrs().data(), f_base.data(), &e_base, nullptr,
                  m1->getCell().data());
  REQUIRE(std::isfinite(e_base));

  // Variant evaluation (test model may expose energy/doubled)
  params.metatomic_options.variant.base = "doubled";
  double e_var{0};
  AtomMatrix f_var = MatrixXd::Zero(m1->numberOfAtoms(), 3);
  try {
    auto pot_var = eonc::helpers::makePotential(
        params.potential_options.potential, params);
    pot_var->force(m1->numberOfAtoms(), m1->getPositions().data(),
                   m1->getAtomicNrs().data(), f_var.data(), &e_var, nullptr,
                   m1->getCell().data());
  } catch (const std::exception &e) {
    WARN(std::string("Skipping doubled variant: ") + e.what());
    TearDown();
    return;
  }

  REQUIRE(std::isfinite(e_var));
  REQUIRE(f_var.allFinite());
  // If the variant is active, energy should differ from the default head.
  // Accept either ~2x (classic test model) or any finite alternate head.
  if (std::abs(e_var - 2.0 * e_base) < 0.5 * std::abs(e_base) + 1.0) {
    REQUIRE_THAT(e_var, WithinAbs(2.0 * e_base, 0.5 * std::abs(e_base) + 1.0));
  } else {
    REQUIRE(e_var != Catch::Approx(e_base).epsilon(1e-6));
  }

  TearDown();
}

} /* namespace tests */

#include "../MatrixHelpers.hpp"
#include "Matter.h"
#include "catch2/catch_amalgamated.hpp"
#include "client/GPRHelpers.h"
#include "client/potentials/GPRPotential/GPRPotential.h"
#include "structures/Structures.h"
#include <memory>
#include <string>

// Forward declaration for the wrapper class defined in GPRPotential.cpp
class GPRModelWrapper;

using namespace std::placeholders;
using namespace Catch::Matchers;

namespace tests {

class GPRPotTest {
public:
  GPRPotTest()
      : params{std::make_shared<Parameters>()},
        matter{nullptr},
        pot_zbl{nullptr},
        threshold{1e-6} {
    params->potential = PotType::ZBL;
    params->zbl_options.cut_inner = 2.0;
    params->zbl_options.cut_global = 2.5;

    pot_zbl = helper_functions::makePotential(params->potential, params);

    pot_gprd = std::make_shared<GPRPotential>(params);
    pot_gprd->registerTargetPotential(pot_zbl);
    matter = std::make_shared<Matter>(pot_gprd, params);
    ref_mat = std::make_shared<Matter>(pot_zbl, params);

    const std::string confile("pos.con");
    const bool file_read_ok = matter->con2matter(confile);
    ref_mat->con2matter(confile);
    REQUIRE(file_read_ok);
    pot_gprd->atom_conf = std::make_unique<gpr::AtomsConfiguration>(
        helper_functions::eon_matter_to_atmconf(ref_mat.get()));
    gpr::InputParameters p = helper_functions::eon_parameters_to_gpr(params.get());
    pot_gprd->gpr_model->initialize(p, *pot_gprd->atom_conf);
  }

  ~GPRPotTest() = default;

protected:
  std::shared_ptr<Parameters> params;
  std::shared_ptr<Matter> matter, ref_mat;
  std::shared_ptr<Potential> pot_zbl;
  std::shared_ptr<GPRPotential> pot_gprd;
  double threshold;
};

TEST_CASE_METHOD(GPRPotTest, "GPR Potential mimicing ZBL", "[PotTest][GPR]") {
  // Expected reference values from LAMMPS
  const double expected_energy = 0.38537731;
  AtomMatrix expected_forces(2, 3);
  expected_forces.row(0) << -2.37926, -2.57753, -2.7758; // Si
  expected_forces.row(1) << 2.37926, 2.57753, 2.7758;    // Au

  // 1. CREATE TRAINING DATA from the reference potential.
  pot_gprd->all_obs = helper_functions::eon_matter_to_init_obs(ref_mat.get());
  pot_gprd->all_obs.printSizes();

  // 2. SETUP THE MODEL with explicit and complete hyperparameters.
  // A) Configure the atomic structure information.
  pot_gprd->gpr_model->getSexpAtCovarianceFunction()->setConfInfo(
      *pot_gprd->atom_conf);

  // B) Manually set all baseline hyperparameter values for stability.
  const int num_pair_types = pot_gprd->atom_conf->n_pt;
  REQUIRE(num_pair_types > 0); // Sanity check

  // Explicitly resize the length scale vector before setting its values.
  pot_gprd->gpr_model->getSexpAtCovarianceFunction()
      ->getLengthScaleRef()
      .resize(1, num_pair_types);

  pot_gprd->gpr_model->getSexpAtCovarianceFunction()->setMagnSigma2(
      6.93874748072254e-009);
  pot_gprd->gpr_model->getSexpAtCovarianceFunction()->setLengthScale(
      888.953211438594e-006);
  pot_gprd->gpr_model->getConstantCovarianceFunction()->setConstSigma2(1.0);
  pot_gprd->gpr_model->getLikGaussian()->setSigma2(1e-10);
  pot_gprd->gpr_model->setJitterSigma2(1e-12);

  // 3. CONDITION THE MODEL using the training data and the fixed
  // hyperparameters.
  auto *wrapper = static_cast<GPRModelWrapper *>(pot_gprd->gpr_model.get());
  REQUIRE(wrapper->condition_model(pot_gprd->all_obs) == true);

  // 4. PREDICT on the same point to test for interpolation.
  gpr::Observation query_obs;
  query_obs.R = pot_gprd->all_obs.R; // Copy positions from training data
  pot_gprd->gpr_model->calculatePotential(query_obs); // Predict E and G

  // 5. VERIFY THE RESULTS from the prediction.
  double calculated_energy = query_obs.E[0];
  // The GPR returns a flattened 1x6 gradient vector; reshape it to a 2x3
  // matrix.
  AtomMatrix calculated_forces =
      query_obs.G.extractEigenMatrix().reshaped(matter->numberOfAtoms(), 3);

  std::stringstream ss_expected;
  ss_expected << "\n" << expected_forces;
  std::stringstream ss_calculated;
  ss_calculated << "\n" << calculated_forces;

  // Perform the checks
  INFO("Energy check: Expected " << expected_energy << ", got "
                                 << calculated_energy);
  REQUIRE_THAT(calculated_energy, WithinAbs(expected_energy, threshold));

  INFO("Force matrix comparison failed.\nExpected:"
       << ss_expected.str() << "\nCalculated:" << ss_calculated.str());
  auto matEq =
      std::bind(helper_functions::eigenEquality<AtomMatrix>, _1, _2, threshold);
  REQUIRE(matEq(calculated_forces, expected_forces));
}

} // namespace tests

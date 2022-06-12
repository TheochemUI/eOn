/*
 * GPRTrainTest.cpp
 *
 *  Created on: 27 Mar 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "GPRTrainTest.h"

#include <algorithm>

namespace tests {

GPRTrainTest::GPRTrainTest() {
  this->reactantFilename = helper_functions::getRelevantFile("reactant.con");
  this->productFilename = helper_functions::getRelevantFile("product.con");

  this->parameters = std::make_unique<Parameters>();
  this->parameters->potential = "morse_pt";
  this->parameters->nebImages = 7;
  this->parameters->LogPotential = false;
  log_init(parameters.get(), (char *)"test.log");

  this->initmatter = std::make_unique<Matter>(parameters.get());
  this->finalmatter = std::make_unique<Matter>(parameters.get());

  this->initmatter->con2matter(reactantFilename);
  this->finalmatter->con2matter(productFilename);

  this->gprfunc = std::make_unique<gpr::GaussianProcessRegression>();

  // Constants
  this->init_eref = this->initmatter.get()->getPotentialEnergy();
  this->init_frcsref = this->initmatter.get()->getForcesFree();
  this->config_data = helper_functions::eon_matter_to_frozen_conf_info(
      this->initmatter.get(), 5);
  auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
  auto R_init = std::get<gpr::Coord>(config_data);
  // Setup the observations
  this->imgArray = helper_functions::prepInitialPath(this->parameters.get());
  this->obspath = helper_functions::prepInitialObs(imgArray);
  // Setup GPR
  auto eondat = std::make_pair(*this->parameters, *this->initmatter);
  *this->gprfunc = helper_functions::initializeGPR(*this->gprfunc, atoms_config,
                                                   obspath, eondat);
  this->gprfunc->setHyperparameters(obspath, atoms_config);
  this->gprfunc->optimize(obspath);
  // Comparer
  this->comparer = [&](const gpr::EigenMatrix &lhs,
                       const gpr::EigenMatrix &rhs) -> bool {
    return lhs.isApprox(rhs, 1e-3); // FIXME: lowered because of the gradients
  };
};

GPRTrainTest::~GPRTrainTest(){
    // TODO Auto-generated destructor stub
};

TEST_F(GPRTrainTest, TestSinglePoint) {
  gpr::Observation o;
  o.clear();
  o.R.resize(1, init_frcsref.rows() * init_frcsref.cols());
  o.G.resize(1, init_frcsref.rows() * init_frcsref.cols());
  gpr::EigenMatrix freePos = this->initmatter.get()->getPositionsFree();
  for (size_t idx{0}; idx < init_frcsref.rows() * init_frcsref.cols(); ++idx) {
    o.R(0, idx) = freePos.reshaped<Eigen::RowMajor>()[idx];
  }
  o.E.resize(1);
  this->gprfunc->calculatePotential(o);
  ASSERT_NEAR(o.E.extractEigenMatrix()(0), init_eref, this->threshold * 1e3)
      << "Energy does not match";
  // The first element is actually almost zero but isn't comparing correctly
  ASSERT_PRED2(comparer, o.G.extractEigenMatrix().reshaped<Eigen::RowMajor>(),
               init_frcsref.reshaped<Eigen::RowMajor>().array()*-1)
      << "Gradients don't match";
}

TEST_F(GPRTrainTest, TestMultiPoint) {
  // Multiple observations
  gpr::Observation oo = obspath;
  this->gprfunc->calculatePotential(oo);
  // Reshape them to the same size
  ASSERT_PRED2(comparer, oo.E.extractEigenMatrix().reshaped<Eigen::RowMajor>(),
               obspath.E.extractEigenMatrix().reshaped<Eigen::RowMajor>())
      << "Energies don't match";
  ASSERT_PRED2(comparer, oo.G.extractEigenMatrix().reshaped<Eigen::RowMajor>(),
               obspath.G.extractEigenMatrix().reshaped<Eigen::RowMajor>())
      << "Gradients don't match";
};

// TEST_F(GPRTrainTest, FunctionCalls) {
//   // Constants
//   const auto init_eref = this->initmatter.get()->getPotentialEnergy();
//   const auto init_frcsref = this->initmatter.get()->getForcesFree();
//   aux::ProblemSetUp problem_setup;
//   auto config_data = helper_functions::eon_matter_to_frozen_conf_info(
//       this->initmatter.get(), 5);
//   auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
//   auto R_init = std::get<gpr::Coord>(config_data);
//   // Setup GPR
//   auto eondat = std::make_pair(*this->parameters, *this->initmatter);
//   *this->gprfunc = helper_functions::initializeGPR(
//       *this->gprfunc, atoms_config, obspath, eondat);
//   this->gprfunc->setHyperparameters(obspath, atoms_config);
//   this->gprfunc->optimize(obspath);

//   // Function calls
//   GPRPotential pot{this->parameters.get()};
//   pot.registerGPRObject(this->gprfunc.get());
//   auto egf_gpro =
//       helper_functions::gpr_energy_and_forces(this->initmatter.get(), &pot);
//   auto energy_gpro = std::get<double>(egf_gpro);
//   auto forces_gpro = std::get<AtomMatrix>(egf_gpro);
//   ASSERT_NEAR(energy_gpro, init_eref, this->threshold * 1e3)
//       << "Energy does not match";
//   double testForcesNew = (forces_gpro - init_frcsref).norm();
//   ASSERT_NEAR(testForces, 0, 1e-3) << "Forces do not match";
// }

} /* namespace tests */

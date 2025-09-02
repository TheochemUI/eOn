#ifndef GPRPOT_INTERFACE
#define GPRPOT_INTERFACE

#include "client/Matter.h"
#include "client/Potential.h"
#include "structures/Structures.h"
#include "subprojects/gpr_optim/gpr/ml/GaussianProcessRegression.h"
#include "subprojects/gpr_optim/gpr/auxiliary/ProblemSetUp.h"
#include "subprojects/gpr_optim/gpr/auxiliary/AdditionalFunctionality.h"
#include <memory>

// A wrapper class to expose the protected conditioning methods from the base
// class. This is necessary because the public API does not provide a standalone
// "condition" method.
class GPRModelWrapper : public gpr::GaussianProcessRegression {
public:
  // This function performs the sequence of protected calls required to
  // condition the model. This logic is copied directly from the library's own
  // reference tests.
  bool condition_model(const gpr::Observation &training_data) {
    aux::AuxiliaryFunctionality aux_func;

    // 1. Convert the training data into the internal matrix format.
    aux_func.assembleMatrixOfRepetitiveCoordinates(
        training_data.R, this->R_matrix, this->R_indices);
    aux_func.assembleVectorFromEnergyAndGradient(training_data,
                                                 this->energy_and_gradient);

    // 2. Build and decompose the covariance matrix.
    if (!this->decomposeCovarianceMatrix(this->R_matrix, this->R_indices)) {
      return false; // Return false if decomposition fails.
    }

    // 3. Solve for the internal model weights.
    this->calculateMeanPrediction(this->energy_and_gradient);
    this->calculatePosteriorMeanPrediction();

    return true;
  }
};

class GPRPotential : public Potential {

private:
  // True or target potential
  bool is_initialized;

public:
  explicit GPRPotential(std::shared_ptr<Parameters> p);
  ~GPRPotential() = default;

  void registerTargetPotential(std::shared_ptr<Potential> _tpot);
  void train_optimize();

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
  // TODO(rg): These must move
  std::unique_ptr<gpr::AtomsConfiguration> atom_conf;
  std::unique_ptr<GPRModelWrapper> gpr_model;
  gpr::Observation all_obs;
  std::shared_ptr<Potential> tpot;
  // void setAtomsConfig(std::shared_ptr<Matter> _mat);
};
#endif

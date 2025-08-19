#include "GPRPotential.h"
#include "subprojects/gpr_optim/structures/Structures.h"

GPRPotential::GPRPotential(std::shared_ptr<Parameters> p)
    : Potential(p),
      is_initialized(false) {
  gpr_model = std::make_unique<GPRModelWrapper>();
  gpr::GPRSetup gpr_parameters;
  gpr_model->setParameters(gpr_parameters);
  // Actually
  // gpr_model->initialize(parameters, atoms_config);
}

void GPRPotential::registerTargetPotential(std::shared_ptr<Potential> _tpot) {
  this->tpot = _tpot;
  if (atom_conf) {
    is_initialized = true;
  }
}

// void GPRPotential::setAtomsConfig(std::shared_ptr<Matter> _mat) {
//   atom_conf = std::make_unique<gpr::AtomsConfiguration>(
//       helper_functions::eon_matter_to_atmconf(_mat.get()));
//   all_obs = helper_functions::eon_matter_to_init_obs(_mat.get());
//   if (tpot) {
//     is_initialized = true;
//   }
// }

void GPRPotential::train_optimize() {
  gpr_model->getSexpAtCovarianceFunction()->setConfInfo(*atom_conf);
  gpr_model->setHyperparameters(all_obs, *atom_conf);
  gpr_model->optimize(all_obs);
}

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void GPRPotential::force(long N, const double *R, const int *atomicNrs,
                         double *F, double *U, double *variance,
                         const double *box) {
  variance = nullptr;
  gpr::Observation observation;

  // Copy R points. Note, R should correspond to the moving atoms only.
  observation.R.resize(1, N * 3);
  for (int i = 0; i < N; i++) {
    observation.R.set(i, {R[3 * i], R[3 * i + 1], R[3 * i + 2]});
  }

  // Note, the following functions should be called before calling for
  // gpr_model->calculatePotential() gpr_model->decomposeCovarianceMatrix(R,
  // ind) - takes covariance matrix and vector of repetitive indices
  // gpr_model->calculateMeanPrediction() - takes a vector of combined energy
  // and force gpr_model->calculatePosteriorMeanPrediction() - no arguments
  gpr_model->getSexpAtCovarianceFunction()->setConfInfo(*atom_conf);
  gpr_model->calculatePotential(observation);

  for (int i = 0; i < N; i++) {
    F[3 * i] = observation.G[3 * i];
    F[3 * i + 1] = observation.G[3 * i + 1];
    F[3 * i + 2] = observation.G[3 * i + 2];
  }

  // FIXME: Test conversion, E should only have one element here
  *U = observation.E[0];
}

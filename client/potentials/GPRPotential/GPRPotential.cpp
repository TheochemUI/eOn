//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "GPRPotential.h"

void GPRPotential::prepare(std::shared_ptr<Matter> a_matter) {
  auto [initAtmConf, RinitObs] =
      helpers::gproptim::input::eon_matter_to_frozen_conf_info(
          a_matter, m_params->gprActiveRadius);
  m_atmconf = initAtmConf;
  for (int i = 0; i < 9; i++) {
    m_inp_gpp.cell_dimensions.value[i] = a_matter->getCell()(i);
  }
  m_gprm.initialize(m_inp_gpp, m_atmconf);
}

void GPRPotential::train_optimize(Eigen::MatrixXd a_features,
                                  Eigen::MatrixXd a_targets) {
  const size_t n_rows(a_features.rows()), n_feature_cols(a_features.cols());
  auto energies{a_targets.block(0, 0, n_rows, 1)};
  auto gradients{a_targets.block(0, 1, n_rows, n_feature_cols)};
  gpr::Observation obs;
  obs.R.resize(n_rows, n_feature_cols);
  obs.G.resize(n_rows, n_feature_cols);
  obs.E.resize(n_rows);
  obs.E.assignFromEigenMatrix(energies);
  obs.R.assignFromEigenMatrix(a_features);
  obs.G.assignFromEigenMatrix(gradients);
  m_gprm.setHyperparameters(obs, m_atmconf);
  m_gprm.optimize(obs);
  return;
}

void GPRPotential::force(long N, const double *R, const int *atomicNrs,
                         double *F, double *U, double *variance,
                         const double *box) {
  gpr::Observation eg_obs, var_obs;

  // Copy R points. Note, R should correspond to the moving atoms only.
  eg_obs.R.resize(1, N * 3);
  for (int i = 0; i < N; i++) {
    eg_obs.R.set(i, {R[3 * i], R[3 * i + 1], R[3 * i + 2]});
  }

  // Note, the following functions should be called before calling for
  // gpr_model->calculatePotential() gpr_model->decomposeCovarianceMatrix(R,
  // ind) - takes covariance matrix and vector of repetitive indices
  // gpr_model->calculateMeanPrediction() - takes a vector of combined energy
  // and force gpr_model->calculatePosteriorMeanPrediction() - no arguments
  m_gprm.calculatePotential(eg_obs);
  m_gprm.calculateVariance(var_obs);

  for (int i = 0; i < N; i++) {
    F[3 * i] = eg_obs.G[3 * i];
    F[3 * i + 1] = eg_obs.G[3 * i + 1];
    F[3 * i + 2] = eg_obs.G[3 * i + 2];
  }

  // FIXME: Test conversion, E should only have one element here
  *U = eg_obs.E[0];
  *variance = var_obs.E[0];
}

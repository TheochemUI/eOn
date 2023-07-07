//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef GPRPOT_INTERFACE
#define GPRPOT_INTERFACE

#include "../../GPRHelpers.h"
#include "../../SurrogatePotential.h"
#include "../../subprojects/gpr_optim/gpr/Enums.h"
#include "../../subprojects/gpr_optim/gpr/auxiliary/AdditionalFunctionality.h"
#include "../../subprojects/gpr_optim/gpr/ml/GaussianProcessRegression.h"
#include "../../subprojects/gpr_optim/structures/Structures.h"

/** GPR Potential from GP Optim. */
class GPRPotential final : public SurrogatePotential {

private:
  gpr::GaussianProcessRegression m_gprm;
  gpr::AtomsConfiguration m_atmconf;
  std::shared_ptr<spdlog::logger> m_log;
  gpr::InputParameters m_inp_gpp;

public:
  GPRPotential(std::shared_ptr<Parameters> a_params)
      : SurrogatePotential(PotType::GPR_Optim, a_params) {
    m_inp_gpp = helpers::gproptim::input::eon_parameters_to_gpr(a_params);
    if (spdlog::get("GPROptim")) {
      m_log = spdlog::get("GPROptim");
    } else {
      m_log = spdlog::basic_logger_st("GPROptim", "_GPROptim.log", true);
    }
    m_log->set_pattern("[%l] [GPROPTIM] %v");
    gpr::GPRSetup gpr_parameters;
    gpr_parameters.jitter_sigma2 = m_params->gprDimerJitterSigma2;
    gpr_parameters.sigma2 = m_params->gprDimerSigma2;
    gpr_parameters.optimization_alg = OptimizationAlgorithms::SCG_opt;
    m_gprm.setParameters(gpr_parameters);
  };

  void prepare(std::shared_ptr<Matter> a_gprm);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  void train_optimize(Eigen::MatrixXd features,
                      Eigen::MatrixXd targets) override;
};
#endif

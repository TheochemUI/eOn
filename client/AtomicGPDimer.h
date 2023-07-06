#ifndef ATOMICGPDIMER_H
#define ATOMICGPDIMER_H

#include "Eigen.h"
#include "GPRHelpers.h"
#include "HelperFunctions.h"
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"
#include <vector>

// Defined in the gprdimer target
#include "subprojects/gpr_optim/gpr/AtomicDimer.h"
#include "subprojects/gpr_optim/gpr/Enums.h"
#include "subprojects/gpr_optim/gpr/auxiliary/ProblemSetUp.h"
#include "subprojects/gpr_optim/gpr/covariance_functions/ConstantCF.h"
#include "subprojects/gpr_optim/gpr/covariance_functions/SexpatCF.h"
#include "subprojects/gpr_optim/gpr/ml/GaussianProcessRegression.h"
#include "subprojects/gpr_optim/structures/Structures.h"

// dimer method to find the lowest curvature mode
class AtomicGPDimer final : public LowestEigenmode {

public:
  // Optimization for the dimer
  static const char OPT_SCG[];
  static const char OPT_LBFGS[];

  AtomicGPDimer(std::shared_ptr<Matter> a_matter,
                std::shared_ptr<Parameters> a_params,
                std::shared_ptr<Potential> a_pot)
      : LowestEigenmode(a_pot, a_params),
        m_dimer_center{std::make_shared<Matter>(*a_matter)} {
    m_gp_params = helpers::gproptim::input::eon_parameters_to_gpr(params);
    for (int i = 0; i < 9; i++) {
      m_gp_params.cell_dimensions.value[i] = a_matter->getCell()(i);
    }
    m_log = spdlog::get("combi");
  }
  ~AtomicGPDimer() = default;

  void compute(std::shared_ptr<Matter> a_matter,
               AtomMatrix a_initialDirectionMatrix) override;
  double getEigenvalue() override;
  AtomMatrix getEigenvector() override;

private:
  std::shared_ptr<Matter> m_dimer_center; // initial center of the dimer
  AtomMatrix m_direction;                 // direction along the dimer
  AtomMatrix
      m_rotationalPlane; // direction normal to the plane of dimer rotation

  gpr::InputParameters m_gp_params;
  atmd::AtomicDimer m_atomic_dimer;
  aux::ProblemSetUp m_problem_setup;
  gpr::AtomsConfiguration m_atoms_config;
  gpr::Observation m_init_observations, m_init_middle_point;
  gpr::Coord m_orient_init, m_R_init;
  shared_ptr<spdlog::logger> m_log;
};

#endif

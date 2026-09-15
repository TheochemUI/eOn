/*
** Apply defaults originating from schema/eon_params.capnp (via codegen).
*/
#include "eon/ParametersSSOT.h"
#include "eon/BaseStructures.h"
#include "eon/generated/ParametersSSOTDefaults.h"

#include <string>
#include <unordered_set>

namespace eonc::config {
namespace {

using GD = eonc::params_ssot::GeneratedDefaults;

// Catalog field index for ssot_has_field (mirrors generated catalog sections).
const std::unordered_set<std::string> &field_index() {
  static const std::unordered_set<std::string> idx = {
#include "eon/generated/ParametersSSOTFieldIndex.inc"
  };
  return idx;
}

JobType job_from_ssot(std::string_view j) {
  if (j == "process_search")
    return JobType::Process_Search;
  if (j == "minimization")
    return JobType::Minimization;
  if (j == "saddle_search")
    return JobType::Saddle_Search;
  if (j == "basin_hopping")
    return JobType::Basin_Hopping;
  if (j == "parallel_replica" || j == "unbiased_parallel_replica")
    return JobType::Parallel_Replica;
  if (j == "nudged_elastic_band")
    return JobType::Nudged_Elastic_Band;
  if (j == "dynamics" || j == "molecular_dynamics")
    return JobType::Dynamics;
  if (j == "hessian")
    return JobType::Hessian;
  if (j == "point")
    return JobType::Point;
  if (j == "prefactor")
    return JobType::Prefactor;
  if (j == "monte_carlo")
    return JobType::Monte_Carlo;
  if (j == "structure_comparison")
    return JobType::Structure_Comparison;
  if (j == "gp_surrogate")
    return JobType::GP_Surrogate;
  if (j == "safe_hyperdynamics")
    return JobType::Safe_Hyperdynamics;
  if (j == "tad")
    return JobType::TAD;
  if (j == "replica_exchange")
    return JobType::Replica_Exchange;
  if (j == "finite_difference" || j == "finite_differences")
    return JobType::Finite_Difference;
  if (j == "global_optimization")
    return JobType::Global_Optimization;
  return JobType::Process_Search;
}

PotType pot_from_ssot(std::string_view pot) {
  if (pot == "lj")
    return PotType::LJ;
  if (pot == "eam_al")
    return PotType::EAM_AL;
  if (pot == "emt")
    return PotType::EMT;
  if (pot == "lammps")
    return PotType::LAMMPS;
  if (pot == "morse_pt")
    return PotType::MORSE_PT;
  if (pot == "metatomic")
    return PotType::METATOMIC;
  if (pot == "xtb")
    return PotType::XTB;
  if (pot == "rgpot")
    return PotType::RGPOT;
  if (pot == "ext_pot")
    return PotType::EXT_POT;
  return PotType::LJ;
}

OptType opt_from_ssot(std::string_view m) {
  if (m == "cg")
    return OptType::CG;
  if (m == "lbfgs")
    return OptType::LBFGS;
  if (m == "qm" || m == "quickmin")
    return OptType::QM;
  if (m == "sd")
    return OptType::SD;
  if (m == "fire")
    return OptType::FIRE;
  return OptType::CG;
}

} // namespace

void apply_ssot_defaults(Parameters &p) {
  ParametersLoadAccess::main_options(p).job = job_from_ssot(GD::MAIN_JOB);
  ParametersLoadAccess::main_options(p).randomSeed = GD::MAIN_RANDOM_SEED;
  ParametersLoadAccess::main_options(p).temperature = GD::MAIN_TEMPERATURE;
  ParametersLoadAccess::main_options(p).quiet = GD::MAIN_QUIET;
  ParametersLoadAccess::main_options(p).writeLog = GD::MAIN_WRITE_LOG;
  ParametersLoadAccess::main_options(p).checkpoint = GD::MAIN_CHECKPOINT;
  ParametersLoadAccess::main_options(p).iniFilename = std::string(GD::MAIN_INI_FILENAME);
  ParametersLoadAccess::main_options(p).conFilename = std::string(GD::MAIN_CON_FILENAME);
  ParametersLoadAccess::main_options(p).finiteDifference = GD::MAIN_FINITE_DIFFERENCE;
  ParametersLoadAccess::main_options(p).maxForceCalls = GD::MAIN_MAX_FORCE_CALLS;
  ParametersLoadAccess::main_options(p).removeNetForce = GD::MAIN_REMOVE_NET_FORCE;
  ParametersLoadAccess::main_options(p).writeConForces = GD::MAIN_WRITE_CON_FORCES;
  ParametersLoadAccess::main_options(p).parallel = GD::MAIN_PARALLEL;

  ParametersLoadAccess::potential_options(p).potential = pot_from_ssot(GD::POTENTIAL_POTENTIAL);
  ParametersLoadAccess::potential_options(p).MPIPollPeriod = GD::POTENTIAL_MPI_POLL_PERIOD;
  ParametersLoadAccess::potential_options(p).LAMMPSLogging = GD::POTENTIAL_LAMMPS_LOGGING;
  ParametersLoadAccess::potential_options(p).LAMMPSThreads = GD::POTENTIAL_LAMMPS_THREADS;
  ParametersLoadAccess::potential_options(p).EMTRasmussen = GD::POTENTIAL_EMT_RASMUSSEN;
  ParametersLoadAccess::potential_options(p).LogPotential = GD::POTENTIAL_LOG_POTENTIAL;
  ParametersLoadAccess::potential_options(p).extPotPath = std::string(GD::POTENTIAL_EXT_POT_PATH);
  ParametersLoadAccess::potential_options(p).potentialsPath =
      std::string(GD::POTENTIAL_POTENTIALS_PATH);

  ParametersLoadAccess::structure_comparison_options(p).distance_difference =
      GD::STRUCTURE_COMPARISON_DISTANCE_DIFFERENCE;
  ParametersLoadAccess::structure_comparison_options(p).neighbor_cutoff =
      GD::STRUCTURE_COMPARISON_NEIGHBOR_CUTOFF;
  ParametersLoadAccess::structure_comparison_options(p).check_rotation =
      GD::STRUCTURE_COMPARISON_CHECK_ROTATION;
  ParametersLoadAccess::structure_comparison_options(p).indistinguishable_atoms =
      GD::STRUCTURE_COMPARISON_INDISTINGUISHABLE_ATOMS;
  ParametersLoadAccess::structure_comparison_options(p).energy_difference =
      GD::STRUCTURE_COMPARISON_ENERGY_DIFFERENCE;
  ParametersLoadAccess::structure_comparison_options(p).remove_translation =
      GD::STRUCTURE_COMPARISON_REMOVE_TRANSLATION;

  ParametersLoadAccess::process_search_options(p).minimize_first = GD::PROCESS_SEARCH_MINIMIZE_FIRST;
  ParametersLoadAccess::process_search_options(p).minimization_offset =
      GD::PROCESS_SEARCH_MINIMIZATION_OFFSET;

  ParametersLoadAccess::optimizer_options(p).method = opt_from_ssot(GD::OPTIMIZER_OPT_METHOD);
  ParametersLoadAccess::optimizer_options(p).convergence_metric =
      std::string(GD::OPTIMIZER_CONVERGENCE_METRIC);
  ParametersLoadAccess::optimizer_options(p).max_iterations = GD::OPTIMIZER_MAX_ITERATIONS;
  ParametersLoadAccess::optimizer_options(p).max_move = GD::OPTIMIZER_MAX_MOVE;
  ParametersLoadAccess::optimizer_options(p).converged_force = GD::OPTIMIZER_CONVERGED_FORCE;
  ParametersLoadAccess::optimizer_options(p).time_step_input = GD::OPTIMIZER_TIME_STEP;
  ParametersLoadAccess::optimizer_options(p).max_time_step_input = GD::OPTIMIZER_MAX_TIME_STEP;

  ParametersLoadAccess::optimizer_options(p).lbfgs.memory = GD::OPTIMIZER_LBFGS_MEMORY;
  ParametersLoadAccess::optimizer_options(p).lbfgs.inverse_curvature =
      GD::OPTIMIZER_LBFGS_INVERSE_CURVATURE;
  ParametersLoadAccess::optimizer_options(p).lbfgs.max_inverse_curvature =
      GD::OPTIMIZER_LBFGS_MAX_INVERSE_CURVATURE;
  ParametersLoadAccess::optimizer_options(p).lbfgs.auto_scale = GD::OPTIMIZER_LBFGS_AUTO_SCALE;
  ParametersLoadAccess::optimizer_options(p).lbfgs.angle_reset = GD::OPTIMIZER_LBFGS_ANGLE_RESET;
  ParametersLoadAccess::optimizer_options(p).lbfgs.distance_reset = GD::OPTIMIZER_LBFGS_DISTANCE_RESET;

  ParametersLoadAccess::optimizer_options(p).cg.no_overshooting = GD::OPTIMIZER_CG_NO_OVERSHOOTING;
  ParametersLoadAccess::optimizer_options(p).cg.knock_out_max_move =
      GD::OPTIMIZER_CG_KNOCK_OUT_MAX_MOVE;
  ParametersLoadAccess::optimizer_options(p).cg.line_search = GD::OPTIMIZER_CG_LINE_SEARCH;
  ParametersLoadAccess::optimizer_options(p).cg.line_converged = GD::OPTIMIZER_CG_LINE_CONVERGED;
  ParametersLoadAccess::optimizer_options(p).cg.line_search_max_iter =
      GD::OPTIMIZER_CG_LINE_SEARCH_MAX_ITER;
  ParametersLoadAccess::optimizer_options(p).cg.max_iter_before_reset =
      GD::OPTIMIZER_CG_MAX_ITER_BEFORE_RESET;

  ParametersLoadAccess::optimizer_options(p).quickmin.steepest_descent =
      GD::OPTIMIZER_QUICKMIN_STEEPEST_DESCENT;
  ParametersLoadAccess::optimizer_options(p).sd.alpha = GD::OPTIMIZER_SD_ALPHA;
  ParametersLoadAccess::optimizer_options(p).sd.two_point = GD::OPTIMIZER_SD_TWO_POINT;
}

bool ssot_has_field(const char *section, const char *key) {
  if (!section || !key)
    return false;
  std::string id = std::string(section) + "." + key;
  return field_index().count(id) > 0;
}

} // namespace eonc::config

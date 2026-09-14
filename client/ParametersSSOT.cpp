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
  p.main_options.job = job_from_ssot(GD::MAIN_JOB);
  p.main_options.randomSeed = GD::MAIN_RANDOM_SEED;
  p.main_options.temperature = GD::MAIN_TEMPERATURE;
  p.main_options.quiet = GD::MAIN_QUIET;
  p.main_options.writeLog = GD::MAIN_WRITE_LOG;
  p.main_options.checkpoint = GD::MAIN_CHECKPOINT;
  p.main_options.iniFilename = std::string(GD::MAIN_INI_FILENAME);
  p.main_options.conFilename = std::string(GD::MAIN_CON_FILENAME);
  p.main_options.finiteDifference = GD::MAIN_FINITE_DIFFERENCE;
  p.main_options.maxForceCalls = GD::MAIN_MAX_FORCE_CALLS;
  p.main_options.removeNetForce = GD::MAIN_REMOVE_NET_FORCE;
  p.main_options.writeConForces = GD::MAIN_WRITE_CON_FORCES;
  p.main_options.parallel = GD::MAIN_PARALLEL;

  p.potential_options.potential = pot_from_ssot(GD::POTENTIAL_POTENTIAL);
  p.potential_options.MPIPollPeriod = GD::POTENTIAL_MPI_POLL_PERIOD;
  p.potential_options.LAMMPSLogging = GD::POTENTIAL_LAMMPS_LOGGING;
  p.potential_options.LAMMPSThreads = GD::POTENTIAL_LAMMPS_THREADS;
  p.potential_options.EMTRasmussen = GD::POTENTIAL_EMT_RASMUSSEN;
  p.potential_options.LogPotential = GD::POTENTIAL_LOG_POTENTIAL;
  p.potential_options.extPotPath = std::string(GD::POTENTIAL_EXT_POT_PATH);
  p.potential_options.potentialsPath =
      std::string(GD::POTENTIAL_POTENTIALS_PATH);

  p.structure_comparison_options.distance_difference =
      GD::STRUCTURE_COMPARISON_DISTANCE_DIFFERENCE;
  p.structure_comparison_options.neighbor_cutoff =
      GD::STRUCTURE_COMPARISON_NEIGHBOR_CUTOFF;
  p.structure_comparison_options.check_rotation =
      GD::STRUCTURE_COMPARISON_CHECK_ROTATION;
  p.structure_comparison_options.indistinguishable_atoms =
      GD::STRUCTURE_COMPARISON_INDISTINGUISHABLE_ATOMS;
  p.structure_comparison_options.energy_difference =
      GD::STRUCTURE_COMPARISON_ENERGY_DIFFERENCE;
  p.structure_comparison_options.remove_translation =
      GD::STRUCTURE_COMPARISON_REMOVE_TRANSLATION;

  p.process_search_options.minimize_first = GD::PROCESS_SEARCH_MINIMIZE_FIRST;
  p.process_search_options.minimization_offset =
      GD::PROCESS_SEARCH_MINIMIZATION_OFFSET;

  p.optimizer_options.method = opt_from_ssot(GD::OPTIMIZER_OPT_METHOD);
  p.optimizer_options.convergence_metric =
      std::string(GD::OPTIMIZER_CONVERGENCE_METRIC);
  p.optimizer_options.max_iterations = GD::OPTIMIZER_MAX_ITERATIONS;
  p.optimizer_options.max_move = GD::OPTIMIZER_MAX_MOVE;
  p.optimizer_options.converged_force = GD::OPTIMIZER_CONVERGED_FORCE;
  p.optimizer_options.time_step_input = GD::OPTIMIZER_TIME_STEP;
  p.optimizer_options.max_time_step_input = GD::OPTIMIZER_MAX_TIME_STEP;

  p.optimizer_options.lbfgs.memory = GD::OPTIMIZER_LBFGS_MEMORY;
  p.optimizer_options.lbfgs.inverse_curvature =
      GD::OPTIMIZER_LBFGS_INVERSE_CURVATURE;
  p.optimizer_options.lbfgs.max_inverse_curvature =
      GD::OPTIMIZER_LBFGS_MAX_INVERSE_CURVATURE;
  p.optimizer_options.lbfgs.auto_scale = GD::OPTIMIZER_LBFGS_AUTO_SCALE;
  p.optimizer_options.lbfgs.angle_reset = GD::OPTIMIZER_LBFGS_ANGLE_RESET;
  p.optimizer_options.lbfgs.distance_reset = GD::OPTIMIZER_LBFGS_DISTANCE_RESET;

  p.optimizer_options.cg.no_overshooting = GD::OPTIMIZER_CG_NO_OVERSHOOTING;
  p.optimizer_options.cg.knock_out_max_move =
      GD::OPTIMIZER_CG_KNOCK_OUT_MAX_MOVE;
  p.optimizer_options.cg.line_search = GD::OPTIMIZER_CG_LINE_SEARCH;
  p.optimizer_options.cg.line_converged = GD::OPTIMIZER_CG_LINE_CONVERGED;
  p.optimizer_options.cg.line_search_max_iter =
      GD::OPTIMIZER_CG_LINE_SEARCH_MAX_ITER;
  p.optimizer_options.cg.max_iter_before_reset =
      GD::OPTIMIZER_CG_MAX_ITER_BEFORE_RESET;

  p.optimizer_options.quickmin.steepest_descent =
      GD::OPTIMIZER_QUICKMIN_STEEPEST_DESCENT;
  p.optimizer_options.sd.alpha = GD::OPTIMIZER_SD_ALPHA;
  p.optimizer_options.sd.two_point = GD::OPTIMIZER_SD_TWO_POINT;
}

bool ssot_has_field(const char *section, const char *key) {
  if (!section || !key)
    return false;
  std::string id = std::string(section) + "." + key;
  return field_index().count(id) > 0;
}

} // namespace eonc::config

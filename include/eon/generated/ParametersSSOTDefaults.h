// AUTO-GENERATED from schema/eon_params.capnp — do not edit.
// Regenerate: python tools/params_ssot/codegen.py
#pragma once
#include <cstdint>
#include <string_view>
namespace eonc::params_ssot {
struct GeneratedDefaults {
  static constexpr auto MAIN_JOB = std::string_view{"process_search"};
  static constexpr auto MAIN_RANDOM_SEED = -1;
  static constexpr auto MAIN_TEMPERATURE = 300.0;
  static constexpr auto MAIN_QUIET = false;
  static constexpr auto MAIN_WRITE_LOG = true;
  static constexpr auto MAIN_CHECKPOINT = false;
  static constexpr auto MAIN_INI_FILENAME = std::string_view{"config.ini"};
  static constexpr auto MAIN_CON_FILENAME = std::string_view{"pos.con"};
  static constexpr auto MAIN_FINITE_DIFFERENCE = 0.01;
  static constexpr auto MAIN_MAX_FORCE_CALLS = 0;
  static constexpr auto MAIN_REMOVE_NET_FORCE = true;
  static constexpr auto MAIN_WRITE_CON_FORCES = false;
  static constexpr auto MAIN_PARALLEL = true;
  static constexpr auto POTENTIAL_POTENTIAL = std::string_view{"lj"};
  static constexpr auto POTENTIAL_MPI_POLL_PERIOD = 0.25;
  static constexpr auto POTENTIAL_LAMMPS_LOGGING = false;
  static constexpr auto POTENTIAL_LAMMPS_THREADS = 0;
  static constexpr auto POTENTIAL_EMT_RASMUSSEN = false;
  static constexpr auto POTENTIAL_LOG_POTENTIAL = false;
  static constexpr auto POTENTIAL_EXT_POT_PATH = std::string_view{"./ext_pot"};
  static constexpr auto POTENTIAL_POTENTIALS_PATH = std::string_view{""};
  static constexpr auto STRUCTURE_COMPARISON_DISTANCE_DIFFERENCE = 0.1;
  static constexpr auto STRUCTURE_COMPARISON_NEIGHBOR_CUTOFF = 3.3;
  static constexpr auto STRUCTURE_COMPARISON_CHECK_ROTATION = false;
  static constexpr auto STRUCTURE_COMPARISON_INDISTINGUISHABLE_ATOMS = true;
  static constexpr auto STRUCTURE_COMPARISON_ENERGY_DIFFERENCE = 0.01;
  static constexpr auto STRUCTURE_COMPARISON_REMOVE_TRANSLATION = true;
  static constexpr auto STRUCTURE_COMPARISON_USE_COVALENT = false;
  static constexpr auto STRUCTURE_COMPARISON_COVALENT_SCALE = 1.3;
  static constexpr auto STRUCTURE_COMPARISON_BRUTE_NEIGHBORS = false;
  static constexpr auto PROCESS_SEARCH_MINIMIZE_FIRST = true;
  static constexpr auto PROCESS_SEARCH_MINIMIZATION_OFFSET = 0.2;
  static constexpr auto OPTIMIZER_LBFGS_MEMORY = 20;
  static constexpr auto OPTIMIZER_LBFGS_INVERSE_CURVATURE = 0.01;
  static constexpr auto OPTIMIZER_LBFGS_MAX_INVERSE_CURVATURE = 0.0;
  static constexpr auto OPTIMIZER_LBFGS_AUTO_SCALE = true;
  static constexpr auto OPTIMIZER_LBFGS_ANGLE_RESET = true;
  static constexpr auto OPTIMIZER_LBFGS_DISTANCE_RESET = true;
  static constexpr auto OPTIMIZER_CG_NO_OVERSHOOTING = false;
  static constexpr auto OPTIMIZER_CG_KNOCK_OUT_MAX_MOVE = false;
  static constexpr auto OPTIMIZER_CG_LINE_SEARCH = false;
  static constexpr auto OPTIMIZER_CG_LINE_CONVERGED = 0.1;
  static constexpr auto OPTIMIZER_CG_LINE_SEARCH_MAX_ITER = 10;
  static constexpr auto OPTIMIZER_CG_MAX_ITER_BEFORE_RESET = 0;
  static constexpr auto OPTIMIZER_QUICKMIN_STEEPEST_DESCENT = false;
  static constexpr auto OPTIMIZER_SD_ALPHA = 0.1;
  static constexpr auto OPTIMIZER_SD_TWO_POINT = false;
  static constexpr auto OPTIMIZER_OPT_METHOD = std::string_view{"cg"};
  static constexpr auto OPTIMIZER_CONVERGENCE_METRIC = std::string_view{"norm"};
  static constexpr auto OPTIMIZER_MAX_ITERATIONS = 1000;
  static constexpr auto OPTIMIZER_MAX_MOVE = 0.2;
  static constexpr auto OPTIMIZER_CONVERGED_FORCE = 0.01;
  static constexpr auto OPTIMIZER_TIME_STEP = 1.0;
  static constexpr auto OPTIMIZER_MAX_TIME_STEP = 2.5;
};
} // namespace eonc::params_ssot

#include "eon/Parameters.h"

#include <magic_enum/magic_enum.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <stdexcept>
#include <string>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_parameters(nb::module_ &m) {
  nb::class_<eonc::Parameters>(m, "Parameters",
                               "Client runtime parameters (config.ini / JSON)")
      .def(nb::init<>())
      .def(
          "load",
          [](eonc::Parameters &self, const std::string &path) {
            if (self.load(path))
              throw std::runtime_error("Parameters.load failed for " + path);
          },
          nb::arg("path"))
      .def(
          "load_json",
          [](eonc::Parameters &self, const std::string &json_str) {
            if (self.load_json(json_str))
              throw std::runtime_error("Parameters.load_json failed");
          },
          nb::arg("json_str"))
      .def("to_json", &eonc::Parameters::to_json)
      // --- Main ---
      .def_prop_rw(
          "job", [](const eonc::Parameters &s) { return s.main_options.job; },
          [](eonc::Parameters &s, eonc::JobType j) { s.main_options.job = j; })
      .def_prop_rw(
          "potential",
          [](const eonc::Parameters &s) {
            return s.potential_options.potential;
          },
          [](eonc::Parameters &s, eonc::PotType p) {
            s.potential_options.potential = p;
          })
      .def_prop_rw(
          "temperature",
          [](const eonc::Parameters &s) { return s.main_options.temperature; },
          [](eonc::Parameters &s, double t) { s.main_options.temperature = t; })
      .def_prop_rw(
          "random_seed",
          [](const eonc::Parameters &s) { return s.main_options.randomSeed; },
          [](eonc::Parameters &s, long v) { s.main_options.randomSeed = v; })
      .def_prop_rw(
          "quiet",
          [](const eonc::Parameters &s) { return s.main_options.quiet; },
          [](eonc::Parameters &s, bool v) { s.main_options.quiet = v; })
      .def_prop_rw(
          "write_log",
          [](const eonc::Parameters &s) { return s.main_options.writeLog; },
          [](eonc::Parameters &s, bool v) { s.main_options.writeLog = v; })
      .def_prop_rw(
          "checkpoint",
          [](const eonc::Parameters &s) { return s.main_options.checkpoint; },
          [](eonc::Parameters &s, bool v) { s.main_options.checkpoint = v; })
      .def_prop_rw(
          "remove_net_force",
          [](const eonc::Parameters &s) {
            return s.main_options.removeNetForce;
          },
          [](eonc::Parameters &s, bool v) {
            s.main_options.removeNetForce = v;
          })
      .def_prop_rw(
          "write_con_forces",
          [](const eonc::Parameters &s) {
            return s.main_options.writeConForces;
          },
          [](eonc::Parameters &s, bool v) {
            s.main_options.writeConForces = v;
          })
      .def_prop_rw(
          "con_filename",
          [](const eonc::Parameters &s) { return s.main_options.conFilename; },
          [](eonc::Parameters &s, const std::string &v) {
            s.main_options.conFilename = v;
          })
      .def_prop_rw(
          "ini_filename",
          [](const eonc::Parameters &s) { return s.main_options.iniFilename; },
          [](eonc::Parameters &s, const std::string &v) {
            s.main_options.iniFilename = v;
          })
      // --- Potential ---
      .def_prop_rw(
          "potentials_path",
          [](const eonc::Parameters &s) {
            return s.potential_options.potentialsPath;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.potential_options.potentialsPath = v;
          })
      .def_prop_rw(
          "ext_pot_path",
          [](const eonc::Parameters &s) {
            return s.potential_options.extPotPath;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.potential_options.extPotPath = v;
          })
      // --- Optimizer ---
      .def_prop_rw(
          "opt_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.optimizer_options.max_iterations);
          },
          [](eonc::Parameters &s, long v) {
            s.optimizer_options.max_iterations = static_cast<size_t>(v);
          })
      .def_prop_rw(
          "opt_converged_force",
          [](const eonc::Parameters &s) {
            return s.optimizer_options.converged_force;
          },
          [](eonc::Parameters &s, double v) {
            s.optimizer_options.converged_force = v;
          })
      .def_prop_rw(
          "opt_max_move",
          [](const eonc::Parameters &s) {
            return s.optimizer_options.max_move;
          },
          [](eonc::Parameters &s, double v) {
            s.optimizer_options.max_move = v;
          })
      // --- Debug ---
      .def_prop_rw(
          "write_movies",
          [](const eonc::Parameters &s) {
            return s.debug_options.write_movies;
          },
          [](eonc::Parameters &s, bool v) { s.debug_options.write_movies = v; })
      // --- Metatomic ---
      .def_prop_rw(
          "metatomic_model_path",
          [](const eonc::Parameters &s) {
            return s.metatomic_options.model_path;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.metatomic_options.model_path = v;
          })
      .def_prop_rw(
          "metatomic_device",
          [](const eonc::Parameters &s) { return s.metatomic_options.device; },
          [](eonc::Parameters &s, const std::string &v) {
            s.metatomic_options.device = v;
          })
      .def_prop_rw(
          "metatomic_length_unit",
          [](const eonc::Parameters &s) {
            return s.metatomic_options.length_unit;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.metatomic_options.length_unit = v;
          })
      .def_prop_rw(
          "metatomic_deterministic",
          [](const eonc::Parameters &s) {
            return s.metatomic_options.deterministic;
          },
          [](eonc::Parameters &s, bool v) {
            s.metatomic_options.deterministic = v;
          })
      .def_prop_rw(
          "metatomic_deterministic_strict",
          [](const eonc::Parameters &s) {
            return s.metatomic_options.deterministic_strict;
          },
          [](eonc::Parameters &s, bool v) {
            s.metatomic_options.deterministic_strict = v;
          })
      // --- NEB ---
      .def_prop_rw(
          "neb_images",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.neb_options.image_count);
          },
          [](eonc::Parameters &s, long v) { s.neb_options.image_count = v; })
      .def_prop_rw(
          "neb_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.neb_options.max_iterations);
          },
          [](eonc::Parameters &s, long v) { s.neb_options.max_iterations = v; })
      .def_prop_rw(
          "neb_force_tolerance",
          [](const eonc::Parameters &s) {
            return s.neb_options.force_tolerance;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.force_tolerance = v;
          })
      .def_prop_rw(
          "neb_init_method",
          [](const eonc::Parameters &s) {
            return s.neb_options.initialization.method;
          },
          [](eonc::Parameters &s, eonc::NEBInit v) {
            s.neb_options.initialization.method = v;
          })
      .def_prop_rw(
          "neb_initial_path",
          [](const eonc::Parameters &s) {
            return s.neb_options.initialization.input_path;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.neb_options.initialization.input_path = v;
          })
      .def_prop_rw(
          "neb_minimize_endpoints",
          [](const eonc::Parameters &s) {
            return s.neb_options.endpoints.minimize;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.endpoints.minimize = v;
          })
      .def_prop_rw(
          "neb_climbing_image",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.enabled;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.climbing_image.enabled = v;
          })
      .def_prop_rw(
          "neb_climbing_converged_only",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.converged_only;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.climbing_image.converged_only = v;
          })
      .def_prop_rw(
          "neb_ci_after",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.trigger_force;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.climbing_image.trigger_force = v;
          })
      .def_prop_rw(
          "neb_ci_after_rel",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.trigger_factor;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.climbing_image.trigger_factor = v;
          })
      .def_prop_rw(
          "neb_energy_weighted",
          [](const eonc::Parameters &s) {
            return s.neb_options.spring.weighting.enabled;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.spring.weighting.enabled = v;
          })
      .def_prop_rw(
          "neb_ew_ksp_min",
          [](const eonc::Parameters &s) {
            return s.neb_options.spring.weighting.k_min;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.spring.weighting.k_min = v;
          })
      .def_prop_rw(
          "neb_ew_ksp_max",
          [](const eonc::Parameters &s) {
            return s.neb_options.spring.weighting.k_max;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.spring.weighting.k_max = v;
          })
      .def_prop_rw(
          "neb_ci_mmf",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.ocineb.use_mmf;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.climbing_image.ocineb.use_mmf = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_after",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.ocineb.trigger_force;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.climbing_image.ocineb.trigger_force = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_after_rel",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.ocineb.trigger_factor;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.climbing_image.ocineb.trigger_factor = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_nsteps",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                s.neb_options.climbing_image.ocineb.max_steps);
          },
          [](eonc::Parameters &s, long v) {
            s.neb_options.climbing_image.ocineb.max_steps = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_angle",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.ocineb.angle_tol;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.climbing_image.ocineb.angle_tol = v;
          })
      .def_prop_rw(
          "neb_mmf_peaks",
          [](const eonc::Parameters &s) {
            return s.neb_options.mmf_peaks.enabled;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.mmf_peaks.enabled = v;
          })
      // --- RgpotPot (incl. backend=metatomic → libmetatomic_engine) ---
      .def_prop_rw(
          "rgpot_backend",
          [](const eonc::Parameters &s) { return s.rgpot_options.backend; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.backend = v;
          })
      .def_prop_rw(
          "rgpot_engine_path",
          [](const eonc::Parameters &s) { return s.rgpot_options.engine_path; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.engine_path = v;
          })
      .def_prop_rw(
          "rgpot_engine_library",
          [](const eonc::Parameters &s) {
            return s.rgpot_options.engine_library;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.engine_library = v;
          })
      .def_prop_rw(
          "rgpot_model_path",
          [](const eonc::Parameters &s) { return s.rgpot_options.model_path; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.model_path = v;
          })
      .def_prop_rw(
          "rgpot_device",
          [](const eonc::Parameters &s) { return s.rgpot_options.device; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.device = v;
          })
      .def_prop_rw(
          "rgpot_length_unit",
          [](const eonc::Parameters &s) { return s.rgpot_options.length_unit; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.length_unit = v;
          })
      .def_prop_rw(
          "rgpot_check_consistency",
          [](const eonc::Parameters &s) {
            return s.rgpot_options.check_consistency;
          },
          [](eonc::Parameters &s, bool v) {
            s.rgpot_options.check_consistency = v;
          })
      .def_prop_rw(
          "rgpot_uncertainty_threshold",
          [](const eonc::Parameters &s) {
            return s.rgpot_options.uncertainty_threshold;
          },
          [](eonc::Parameters &s, double v) {
            s.rgpot_options.uncertainty_threshold = v;
          })
      .def_prop_rw(
          "rgpot_torch_determinism_strict",
          [](const eonc::Parameters &s) {
            return s.rgpot_options.torch_determinism_strict;
          },
          [](eonc::Parameters &s, bool v) {
            s.rgpot_options.torch_determinism_strict = v;
          })
      .def_prop_rw(
          "rgpot_basis",
          [](const eonc::Parameters &s) { return s.rgpot_options.basis; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.basis = v;
          })
      .def_prop_rw(
          "rgpot_theory",
          [](const eonc::Parameters &s) { return s.rgpot_options.theory; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.theory = v;
          })

      // --- Saddle search ---
      .def_prop_rw(
          "saddle_method",
          [](const eonc::Parameters &s) {
            return s.saddle_search_options.method;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.saddle_search_options.method = v;
          })
      .def_prop_rw(
          "saddle_minmode_method",
          [](const eonc::Parameters &s) {
            return s.saddle_search_options.minmode_method;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.saddle_search_options.minmode_method = v;
          },
          "Min-mode solver: dimer | lanczos | davidson | gprdimer")
      .def_prop_rw(
          "saddle_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.saddle_search_options.max_iterations);
          },
          [](eonc::Parameters &s, long v) {
            s.saddle_search_options.max_iterations = v;
          })
      .def_prop_rw(
          "saddle_max_energy",
          [](const eonc::Parameters &s) {
            return s.saddle_search_options.max_energy;
          },
          [](eonc::Parameters &s, double v) {
            s.saddle_search_options.max_energy = v;
          })
      .def_prop_rw(
          "saddle_displace_magnitude",
          [](const eonc::Parameters &s) {
            return s.saddle_search_options.displace_magnitude;
          },
          [](eonc::Parameters &s, double v) {
            s.saddle_search_options.displace_magnitude = v;
          })
      .def_prop_rw(
          "saddle_displace_type",
          [](const eonc::Parameters &s) {
            return s.saddle_search_options.displace_type;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.saddle_search_options.displace_type = v;
          })
      .def_prop_rw(
          "saddle_converged_force",
          [](const eonc::Parameters &s) {
            return s.saddle_search_options.converged_force;
          },
          [](eonc::Parameters &s, double v) {
            s.saddle_search_options.converged_force = v;
          })
      .def_prop_rw(
          "saddle_remove_rotation",
          [](const eonc::Parameters &s) {
            return s.saddle_search_options.remove_rotation;
          },
          [](eonc::Parameters &s, bool v) {
            s.saddle_search_options.remove_rotation = v;
          })
      // --- Dimer ---
      .def_prop_rw(
          "dimer_improved",
          [](const eonc::Parameters &s) { return s.dimer_options.improved; },
          [](eonc::Parameters &s, bool v) { s.dimer_options.improved = v; },
          "True → ImprovedDimer, False → classic Dimer")
      .def_prop_rw(
          "dimer_rotation_angle",
          [](const eonc::Parameters &s) {
            return s.dimer_options.rotation_angle;
          },
          [](eonc::Parameters &s, double v) {
            s.dimer_options.rotation_angle = v;
          })
      .def_prop_rw(
          "dimer_converged_angle",
          [](const eonc::Parameters &s) {
            return s.dimer_options.converged_angle;
          },
          [](eonc::Parameters &s, double v) {
            s.dimer_options.converged_angle = v;
          })
      .def_prop_rw(
          "dimer_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.dimer_options.max_iterations);
          },
          [](eonc::Parameters &s, long v) {
            s.dimer_options.max_iterations = v;
          })
      .def_prop_rw(
          "dimer_opt_method",
          [](const eonc::Parameters &s) { return s.dimer_options.opt_method; },
          [](eonc::Parameters &s, const std::string &v) {
            s.dimer_options.opt_method = v;
          },
          "Rotation optimizer: cg | lbfgs | sd")
      .def_prop_rw(
          "dimer_rotations_max",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.dimer_options.rotations_max);
          },
          [](eonc::Parameters &s, long v) {
            s.dimer_options.rotations_max = v;
          })
      .def_prop_rw(
          "dimer_rotations_min",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.dimer_options.rotations_min);
          },
          [](eonc::Parameters &s, long v) {
            s.dimer_options.rotations_min = v;
          })
      .def_prop_rw(
          "dimer_torque_max",
          [](const eonc::Parameters &s) { return s.dimer_options.torque_max; },
          [](eonc::Parameters &s, double v) { s.dimer_options.torque_max = v; })
      .def_prop_rw(
          "dimer_torque_min",
          [](const eonc::Parameters &s) { return s.dimer_options.torque_min; },
          [](eonc::Parameters &s, double v) { s.dimer_options.torque_min = v; })
      .def_prop_rw(
          "dimer_remove_rotation",
          [](const eonc::Parameters &s) {
            return s.dimer_options.remove_rotation;
          },
          [](eonc::Parameters &s, bool v) {
            s.dimer_options.remove_rotation = v;
          })
      // --- Lanczos / Davidson (PHVA mobile set for Krylov; default All) ---
      .def_prop_rw(
          "lanczos_phva_atoms",
          [](const eonc::Parameters &s) {
            return s.lanczos_options.phva_atoms;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.lanczos_options.phva_atoms = v;
          },
          "PHVA mobile/active atoms for Lanczos Krylov space, or All = free "
          "(not free/fixed)")
      .def_prop_rw(
          "davidson_phva_atoms",
          [](const eonc::Parameters &s) {
            return s.davidson_options.phva_atoms;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.davidson_options.phva_atoms = v;
          },
          "PHVA mobile/active atoms for Davidson Ritz space, or All = free "
          "(not free/fixed)")
      // --- Hessian ---
      .def_prop_rw(
          "hessian_phva_atoms",
          [](const eonc::Parameters &s) {
            return s.hessian_options.phva_atoms;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.hessian_options.phva_atoms = v;
          },
          "PHVA mobile/active atoms for dense FD Hessian, or All = free "
          "(not free/fixed)")
      .def_prop_rw(
          "hessian_zero_freq_value",
          [](const eonc::Parameters &s) {
            return s.hessian_options.zero_freq_value;
          },
          [](eonc::Parameters &s, double v) {
            s.hessian_options.zero_freq_value = v;
          })
      // --- Prefactor ---
      .def_prop_rw(
          "prefactor_rate",
          [](const eonc::Parameters &s) { return s.prefactor_options.rate; },
          [](eonc::Parameters &s, const std::string &v) {
            s.prefactor_options.rate = v;
          },
          "htst | qqhtst")
      .def_prop_rw(
          "prefactor_min_value",
          [](const eonc::Parameters &s) {
            return s.prefactor_options.min_value;
          },
          [](eonc::Parameters &s, double v) {
            s.prefactor_options.min_value = v;
          })
      .def_prop_rw(
          "prefactor_max_value",
          [](const eonc::Parameters &s) {
            return s.prefactor_options.max_value;
          },
          [](eonc::Parameters &s, double v) {
            s.prefactor_options.max_value = v;
          })
      .def_prop_rw(
          "prefactor_within_radius",
          [](const eonc::Parameters &s) {
            return s.prefactor_options.within_radius;
          },
          [](eonc::Parameters &s, double v) {
            s.prefactor_options.within_radius = v;
          })
      .def_prop_rw(
          "prefactor_min_displacement",
          [](const eonc::Parameters &s) {
            return s.prefactor_options.min_displacement;
          },
          [](eonc::Parameters &s, double v) {
            s.prefactor_options.min_displacement = v;
          })
      .def_prop_rw(
          "prefactor_filter_scheme",
          [](const eonc::Parameters &s) {
            return s.prefactor_options.filter_scheme;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.prefactor_options.filter_scheme = v;
          },
          "cutoff | fraction")
      .def_prop_rw(
          "prefactor_filter_fraction",
          [](const eonc::Parameters &s) {
            return s.prefactor_options.filter_fraction;
          },
          [](eonc::Parameters &s, double v) {
            s.prefactor_options.filter_fraction = v;
          })
      .def_prop_rw(
          "prefactor_all_free_atoms",
          [](const eonc::Parameters &s) {
            return s.prefactor_options.all_free_atoms;
          },
          [](eonc::Parameters &s, bool v) {
            s.prefactor_options.all_free_atoms = v;
          })

      // --- Dynamics / MC / BH ---
      .def_prop_rw(
          "dynamics_steps",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.dynamics_options.steps);
          },
          [](eonc::Parameters &s, long v) { s.dynamics_options.steps = v; })
      .def_prop_rw(
          "dynamics_time_step",
          [](const eonc::Parameters &s) {
            return s.dynamics_options.time_step_input;
          },
          [](eonc::Parameters &s, double v) {
            s.dynamics_options.time_step_input = v;
            s.dynamics_options.time_step = v / s.constants.timeUnit;
          })
      .def_prop_rw(
          "monte_carlo_steps",
          [](const eonc::Parameters &s) { return s.monte_carlo_options.steps; },
          [](eonc::Parameters &s, int v) { s.monte_carlo_options.steps = v; })
      .def_prop_rw(
          "monte_carlo_step_size",
          [](const eonc::Parameters &s) {
            return s.monte_carlo_options.step_size;
          },
          [](eonc::Parameters &s, double v) {
            s.monte_carlo_options.step_size = v;
          })
      .def_prop_rw(
          "basin_hopping_steps",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.basin_hopping_options.steps);
          },
          [](eonc::Parameters &s, long v) {
            s.basin_hopping_options.steps = v;
          })
      .def_prop_rw(
          "basin_hopping_displacement",
          [](const eonc::Parameters &s) {
            return s.basin_hopping_options.displacement;
          },
          [](eonc::Parameters &s, double v) {
            s.basin_hopping_options.displacement = v;
          })
      .def_prop_rw(
          "process_search_minimize_first",
          [](const eonc::Parameters &s) {
            return s.process_search_options.minimize_first;
          },
          [](eonc::Parameters &s, bool v) {
            s.process_search_options.minimize_first = v;
          })

      // --- Structure comparison ---

      .def_prop_rw(
          "comp_eps_r",
          [](const eonc::Parameters &s) {
            return s.structure_comparison_options.distance_difference;
          },
          [](eonc::Parameters &s, double v) {
            s.structure_comparison_options.distance_difference = v;
          })
      .def("__repr__", [](const eonc::Parameters &self) {
        return "<Parameters job=" +
               std::string(magic_enum::enum_name(self.main_options.job)) +
               " pot=" +
               std::string(
                   magic_enum::enum_name(self.potential_options.potential)) +
               ">";
      });
}

} // namespace eonc::pybind

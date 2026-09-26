#include "eon/ConFileIO.h"
#include "eon/HelperFunctions.h"
#include "eon/Parameters.h"

#include <magic_enum/magic_enum.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <cctype>
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
          "load_ini_text",
          [](eonc::Parameters &self, const std::string &ini_text) {
            if (self.load_ini_text(ini_text))
              throw std::runtime_error("Parameters.load_ini_text failed");
          },
          nb::arg("ini_text"))
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
          "job",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::main_options(s).job;
          },
          [](eonc::Parameters &s, eonc::JobType j) {
            eonc::ParametersLoadAccess::main_options(s).job = j;
          })
      .def_prop_rw(
          "potential",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::potential_options(s).potential;
          },
          [](eonc::Parameters &s, eonc::PotType p) {
            eonc::ParametersLoadAccess::potential_options(s).potential = p;
          })
      .def_prop_rw(
          "dftd_functional",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dftd_options(s).functional;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::dftd_options(s).functional = v;
          },
          "s-dftd3 / dftd4 method key (default pbe)")
      .def_prop_rw(
          "dftd_atm",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dftd_options(s).atm;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::dftd_options(s).atm = v;
          },
          "Axilrod-Teller-Muto three-body term")
      .def_prop_rw(
          "d3_damping",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dftd_options(s).d3_damping;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::dftd_options(s).d3_damping = v;
          },
          "bj | zero")
      .def_prop_rw(
          "d4_charge",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dftd_options(s).d4_charge;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::dftd_options(s).d4_charge = v;
          })
      .def_prop_rw(
          "expr_expression",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::expr_options(s).expression;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::expr_options(s).expression = v;
          },
          "Lepton expression over named terms, e.g. 0.5*lj + d3")
      .def_prop_rw(
          "expr_terms",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::expr_options(s).terms;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::expr_options(s).terms = v;
          },
          "Comma-separated term names matching the expression")
      .def_prop_rw(
          "mopac_charge",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::mopac_options(s).charge;
          },
          [](eonc::Parameters &s, int v) {
            eonc::ParametersLoadAccess::mopac_options(s).charge = v;
          })
      .def_prop_rw(
          "mopac_spin",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::mopac_options(s).spin;
          },
          [](eonc::Parameters &s, int v) {
            eonc::ParametersLoadAccess::mopac_options(s).spin = v;
          })
      .def_prop_rw(
          "mopac_model",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::mopac_options(s).model;
          },
          [](eonc::Parameters &s, int v) {
            eonc::ParametersLoadAccess::mopac_options(s).model = v;
          },
          "OpenMOPAC model id (4 is AM1)")
      .def_prop_rw(
          "mopac_engine_path",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::mopac_options(s).engine_path;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::mopac_options(s).engine_path = v;
          },
          "libmopacc path; empty uses default search")
      .def_prop_rw(
          "temperature",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::main_options(s).temperature;
          },
          [](eonc::Parameters &s, double t) {
            eonc::ParametersLoadAccess::main_options(s).temperature = t;
          })
      .def_prop_rw(
          "random_seed",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::main_options(s).randomSeed;
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::main_options(s).randomSeed = v;
          })
      .def_prop_rw(
          "quiet",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::main_options(s).quiet;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::main_options(s).quiet = v;
          })
      .def_prop_rw(
          "write_log",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::main_options(s).writeLog;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::main_options(s).writeLog = v;
          })
      .def_prop_rw(
          "checkpoint",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::main_options(s).checkpoint;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::main_options(s).checkpoint = v;
          })
      .def_prop_rw(
          "remove_net_force",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::main_options(s).removeNetForce;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::main_options(s).removeNetForce = v;
          })
      .def_prop_rw(
          "write_con_forces",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::main_options(s).writeConForces;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::main_options(s).writeConForces = v;
          })
      .def_prop_rw(
          "con_filename",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::main_options(s).conFilename;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::main_options(s).conFilename = v;
          })
      .def_prop_rw(
          "ini_filename",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::main_options(s).iniFilename;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::main_options(s).iniFilename = v;
          })
      // --- Potential ---
      .def_prop_rw(
          "emt_rasmussen",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::potential_options(s)
                .EMTRasmussen;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::potential_options(s).EMTRasmussen = v;
          })
      .def_prop_rw(
          "potentials_path",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::potential_options(s)
                .potentialsPath;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::potential_options(s).potentialsPath = v;
          })
      .def_prop_rw(
          "ext_pot_path",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::potential_options(s).extPotPath;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::potential_options(s).extPotPath = v;
          })
      .def_prop_rw(
          "pot_thread_safe",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::potential_options(s).thread_safe;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::potential_options(s).thread_safe = v;
          },
          "When false, NEB/dimer never share one Potential across threads")
      // --- Optimizer ---
      .def_prop_rw(
          "opt_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::optimizer_options(s)
                    .max_iterations);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::optimizer_options(s).max_iterations =
                static_cast<size_t>(v);
          })
      .def_prop_rw(
          "opt_converged_force",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::optimizer_options(s)
                .converged_force;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::optimizer_options(s).converged_force =
                v;
          })
      .def_prop_rw(
          "opt_max_move",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::optimizer_options(s).max_move;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::optimizer_options(s).max_move = v;
          })
      .def_prop_rw(
          "opt_convergence_metric",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::optimizer_options(s)
                .convergence_metric;
          },
          [](eonc::Parameters &s, const std::string &v) {
            if (auto label = eonc::helpers::convergenceMetricLabel(v)) {
              eonc::ParametersLoadAccess::optimizer_options(s)
                  .convergence_metric = v;
              eonc::ParametersLoadAccess::optimizer_options(s)
                  .convergence_metric_label = std::string(*label);
            } else {
              throw std::invalid_argument("unknown opt_convergence_metric: " +
                                          v);
            }
          },
          "norm | max_atom | max_component")
      .def_prop_rw(
          "opt_method",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::optimizer_options(s).method;
          },
          [](eonc::Parameters &s, eonc::OptType v) {
            eonc::ParametersLoadAccess::optimizer_options(s).method = v;
          },
          "Minimization / default job optimizer (CG, LBFGS, FIRE, QM, SD)")
      .def_prop_rw(
          "refine_opt_method",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::optimizer_options(s)
                .refine.method;
          },
          [](eonc::Parameters &s, eonc::OptType v) {
            eonc::ParametersLoadAccess::optimizer_options(s).refine.method = v;
          },
          "Switch optimizer when force drops below refine_threshold; "
          "OptType.None_ disables")
      .def_prop_rw(
          "refine_threshold",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::optimizer_options(s)
                .refine.threshold;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::optimizer_options(s).refine.threshold =
                v;
          })
      // --- Debug ---
      .def_prop_rw(
          "write_movies",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::debug_options(s).write_movies;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::debug_options(s).write_movies = v;
          })
      // --- Metatomic ---
      .def_prop_rw(
          "metatomic_model_path",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::metatomic_options(s).model_path;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::metatomic_options(s).model_path = v;
          })
      .def_prop_rw(
          "metatomic_device",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::metatomic_options(s).device;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::metatomic_options(s).device = v;
          })
      .def_prop_rw(
          "metatomic_length_unit",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::metatomic_options(s).length_unit;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::metatomic_options(s).length_unit = v;
          })
      .def_prop_rw(
          "metatomic_deterministic",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::metatomic_options(s)
                .deterministic;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::metatomic_options(s).deterministic = v;
          })
      .def_prop_rw(
          "metatomic_deterministic_strict",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::metatomic_options(s)
                .deterministic_strict;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::metatomic_options(s)
                .deterministic_strict = v;
          })
      // --- NEB ---
      .def_prop_rw(
          "neb_images",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::neb_options(s).image_count);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::neb_options(s).image_count = v;
          })
      .def_prop_rw(
          "neb_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::neb_options(s).max_iterations);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::neb_options(s).max_iterations = v;
          })
      .def_prop_rw(
          "neb_opt_method",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s).opt_method;
          },
          [](eonc::Parameters &s, eonc::OptType v) {
            eonc::ParametersLoadAccess::neb_options(s).opt_method = v;
          },
          "NEB band optimizer (independent of opt_method)")
      .def_prop_rw(
          "neb_force_tolerance",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s).force_tolerance;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s).force_tolerance = v;
          })
      .def_prop_rw(
          "neb_init_method",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .initialization.method;
          },
          [](eonc::Parameters &s, eonc::NEBInit v) {
            eonc::ParametersLoadAccess::neb_options(s).initialization.method =
                v;
          })
      .def_prop_rw(
          "neb_initial_path",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .initialization.input_path;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .initialization.input_path = v;
          })
      .def_prop_rw(
          "neb_minimize_endpoints",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .endpoints.minimize;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::neb_options(s).endpoints.minimize = v;
          })
      .def_prop_rw(
          "neb_climbing_image",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.enabled;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::neb_options(s).climbing_image.enabled =
                v;
          })
      .def_prop_rw(
          "neb_climbing_converged_only",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.converged_only;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.converged_only = v;
          })
      .def_prop_rw(
          "neb_climbing_band_slack",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.band_slack;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.band_slack = v;
          })
      .def_prop_rw(
          "neb_ci_after",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.trigger_force;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.trigger_force = v;
          })
      .def_prop_rw(
          "neb_ci_after_rel",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.trigger_factor;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.trigger_factor = v;
          })
      .def_prop_rw(
          "neb_energy_weighted",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .spring.weighting.enabled;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .spring.weighting.enabled = v;
          })
      .def_prop_rw(
          "neb_ew_ksp_min",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .spring.weighting.k_min;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s).spring.weighting.k_min =
                v;
          })
      .def_prop_rw(
          "neb_ew_ksp_max",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .spring.weighting.k_max;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s).spring.weighting.k_max =
                v;
          })
      .def_prop_rw(
          "neb_ci_mmf",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.ocineb.use_mmf;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.ocineb.use_mmf = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_after",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.ocineb.trigger_force;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.ocineb.trigger_force = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_after_rel",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.ocineb.trigger_factor;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.ocineb.trigger_factor = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_nsteps",
          [](const eonc::Parameters &s) {
            return static_cast<long>(eonc::ParametersLoadAccess::neb_options(s)
                                         .climbing_image.ocineb.max_steps);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.ocineb.max_steps = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_angle",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.ocineb.angle_tol;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .climbing_image.ocineb.angle_tol = v;
          })
      .def_prop_rw(
          "neb_mmf_peaks",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s).mmf_peaks.enabled;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::neb_options(s).mmf_peaks.enabled = v;
          })
      .def_prop_rw(
          "neb_match_endpoints",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s).match_endpoints;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::neb_options(s).match_endpoints = v;
          },
          "IRA permute+rotate reactant onto product before NEB interpolation")
      .def_prop_rw(
          "neb_zoom",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s).zoom.enabled;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::neb_options(s).zoom.enabled = v;
          },
          "Pack images onto a window around the climbing image")
      .def_prop_rw(
          "neb_zoom_alpha",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s).zoom.alpha;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s).zoom.alpha = v;
          })
      .def_prop_rw(
          "neb_zoom_offset",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::neb_options(s).zoom.offset);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::neb_options(s).zoom.offset =
                static_cast<int>(v);
          })
      .def_prop_rw(
          "neb_zoom_mode",
          [](const eonc::Parameters &s) {
            auto name = std::string(magic_enum::enum_name(
                eonc::ParametersLoadAccess::neb_options(s).zoom.mode));
            for (char &c : name) {
              c = static_cast<char>(
                  std::tolower(static_cast<unsigned char>(c)));
            }
            return name;
          },
          [](eonc::Parameters &s, const std::string &v) {
            auto mode = magic_enum::enum_cast<
                eonc::neb_options_t::zoom_options_t::Mode>(
                v, magic_enum::case_insensitive);
            if (!mode) {
              throw std::invalid_argument(
                  "neb_zoom_mode must be auto or manual");
            }
            eonc::ParametersLoadAccess::neb_options(s).zoom.mode = *mode;
          })
      .def_prop_rw(
          "neb_zoom_after",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::neb_options(s)
                .zoom.activation_threshold;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::neb_options(s)
                .zoom.activation_threshold = v;
          },
          "Force threshold that arms Zoom-NEB; 0 means 10x converged_force")
      .def_prop_rw(
          "neb_zoom_interpolation",
          [](const eonc::Parameters &s) {
            auto name = std::string(magic_enum::enum_name(
                eonc::ParametersLoadAccess::neb_options(s).zoom.interpolation));
            for (char &c : name) {
              c = static_cast<char>(
                  std::tolower(static_cast<unsigned char>(c)));
            }
            return name;
          },
          [](eonc::Parameters &s, const std::string &v) {
            auto how = magic_enum::enum_cast<
                eonc::neb_options_t::zoom_options_t::Interpolation>(
                v, magic_enum::case_insensitive);
            if (!how) {
              throw std::invalid_argument(
                  "neb_zoom_interpolation must be cubic or linear");
            }
            eonc::ParametersLoadAccess::neb_options(s).zoom.interpolation =
                *how;
          })
      .def_prop_rw(
          "neb_zoom_ci_stability",
          [](const eonc::Parameters &s) {
            return static_cast<long>(eonc::ParametersLoadAccess::neb_options(s)
                                         .zoom.stability_count);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::neb_options(s).zoom.stability_count =
                static_cast<int>(v);
          })
      .def_prop_rw(
          "neb_zoom_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::neb_options(s).zoom.max_iterations);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::neb_options(s).zoom.max_iterations =
                static_cast<int>(v);
          },
          "Iterations after zoom; 0 keeps the band iteration cap")
      // --- RgpotPot (incl. backend=metatomic → libmetatomic_engine) ---
      .def_prop_rw(
          "rgpot_backend",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s).backend;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::rgpot_options(s).backend = v;
          })
      .def_prop_rw(
          "rgpot_engine_path",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s).engine_path;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::rgpot_options(s).engine_path = v;
          })
      .def_prop_rw(
          "rgpot_engine_library",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s).engine_library;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::rgpot_options(s).engine_library = v;
          })
      .def_prop_rw(
          "rgpot_model_path",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s).model_path;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::rgpot_options(s).model_path = v;
          })
      .def_prop_rw(
          "rgpot_device",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s).device;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::rgpot_options(s).device = v;
          })
      .def_prop_rw(
          "rgpot_length_unit",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s).length_unit;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::rgpot_options(s).length_unit = v;
          })
      .def_prop_rw(
          "rgpot_check_consistency",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s)
                .check_consistency;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::rgpot_options(s).check_consistency = v;
          })
      .def_prop_rw(
          "rgpot_uncertainty_threshold",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s)
                .uncertainty_threshold;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::rgpot_options(s).uncertainty_threshold =
                v;
          })
      .def_prop_rw(
          "rgpot_torch_determinism_strict",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s)
                .torch_determinism_strict;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::rgpot_options(s)
                .torch_determinism_strict = v;
          })
      .def_prop_rw(
          "rgpot_basis",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s).basis;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::rgpot_options(s).basis = v;
          })
      .def_prop_rw(
          "rgpot_theory",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::rgpot_options(s).theory;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::rgpot_options(s).theory = v;
          })

      // --- Saddle search ---
      .def_prop_rw(
          "saddle_method",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::saddle_search_options(s).method;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::saddle_search_options(s).method = v;
          })
      .def_prop_rw(
          "saddle_minmode_method",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::saddle_search_options(s)
                .minmode_method;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::saddle_search_options(s)
                .minmode_method = v;
          },
          "Min-mode solver: dimer | lanczos | davidson | gprdimer")
      .def_prop_rw(
          "saddle_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::saddle_search_options(s)
                    .max_iterations);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::saddle_search_options(s)
                .max_iterations = v;
          })
      .def_prop_rw(
          "saddle_max_energy",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::saddle_search_options(s)
                .max_energy;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::saddle_search_options(s).max_energy = v;
          })
      .def_prop_rw(
          "saddle_displace_magnitude",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::saddle_search_options(s)
                .displace_magnitude;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::saddle_search_options(s)
                .displace_magnitude = v;
          })
      .def_prop_rw(
          "saddle_displace_type",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::saddle_search_options(s)
                .displace_type;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::saddle_search_options(s).displace_type =
                v;
          })
      .def_prop_rw(
          "saddle_converged_force",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::saddle_search_options(s)
                .converged_force;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::saddle_search_options(s)
                .converged_force = v;
          })
      .def_prop_rw(
          "saddle_remove_rotation",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::saddle_search_options(s)
                .remove_rotation;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::saddle_search_options(s)
                .remove_rotation = v;
          })
      // --- Dimer ---
      .def_prop_rw(
          "dimer_improved",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dimer_options(s).improved;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::dimer_options(s).improved = v;
          },
          "True → ImprovedDimer, False → classic Dimer")
      .def_prop_rw(
          "dimer_rotation_angle",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dimer_options(s).rotation_angle;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::dimer_options(s).rotation_angle = v;
          })
      .def_prop_rw(
          "dimer_converged_angle",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dimer_options(s).converged_angle;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::dimer_options(s).converged_angle = v;
          })
      .def_prop_rw(
          "dimer_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::dimer_options(s).max_iterations);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::dimer_options(s).max_iterations = v;
          })
      .def_prop_rw(
          "dimer_opt_method",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dimer_options(s).opt_method;
          },
          [](eonc::Parameters &s, eonc::OptType v) {
            eonc::ParametersLoadAccess::dimer_options(s).opt_method = v;
          },
          "Rotation optimizer (CG, LBFGS, SD)")
      .def_prop_rw(
          "dimer_rotations_max",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::dimer_options(s).rotations_max);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::dimer_options(s).rotations_max = v;
          })
      .def_prop_rw(
          "dimer_rotations_min",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::dimer_options(s).rotations_min);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::dimer_options(s).rotations_min = v;
          })
      .def_prop_rw(
          "dimer_torque_max",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dimer_options(s).torque_max;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::dimer_options(s).torque_max = v;
          })
      .def_prop_rw(
          "dimer_torque_min",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dimer_options(s).torque_min;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::dimer_options(s).torque_min = v;
          })
      .def_prop_rw(
          "dimer_remove_rotation",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dimer_options(s).remove_rotation;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::dimer_options(s).remove_rotation = v;
          })
      // --- Lanczos / Davidson (PHVA mobile set for Krylov; default All) ---
      .def_prop_rw(
          "lanczos_phva_atoms",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::lanczos_options(s).phva_atoms;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::lanczos_options(s).phva_atoms = v;
          },
          "PHVA mobile/active atoms for Lanczos Krylov space, or All = free "
          "(not free/fixed)")
      .def_prop_rw(
          "davidson_phva_atoms",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::davidson_options(s).phva_atoms;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::davidson_options(s).phva_atoms = v;
          },
          "PHVA mobile/active atoms for Davidson Ritz space, or All = free "
          "(not free/fixed)")
      // --- Hessian ---
      .def_prop_rw(
          "hessian_phva_atoms",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::hessian_options(s).phva_atoms;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::hessian_options(s).phva_atoms = v;
          },
          "PHVA mobile/active atoms for dense FD Hessian, or All = free "
          "(not free/fixed)")
      .def_prop_rw(
          "hessian_zero_freq_value",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::hessian_options(s)
                .zero_freq_value;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::hessian_options(s).zero_freq_value = v;
          })
      // --- Prefactor ---
      .def_prop_rw(
          "prefactor_rate",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::prefactor_options(s).rate;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::prefactor_options(s).rate = v;
          },
          "htst | qqhtst")
      .def_prop_rw(
          "prefactor_min_value",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::prefactor_options(s).min_value;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::prefactor_options(s).min_value = v;
          })
      .def_prop_rw(
          "prefactor_max_value",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::prefactor_options(s).max_value;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::prefactor_options(s).max_value = v;
          })
      .def_prop_rw(
          "prefactor_within_radius",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::prefactor_options(s)
                .within_radius;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::prefactor_options(s).within_radius = v;
          })
      .def_prop_rw(
          "prefactor_min_displacement",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::prefactor_options(s)
                .min_displacement;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::prefactor_options(s).min_displacement =
                v;
          })
      .def_prop_rw(
          "prefactor_filter_scheme",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::prefactor_options(s)
                .filter_scheme;
          },
          [](eonc::Parameters &s, const std::string &v) {
            eonc::ParametersLoadAccess::prefactor_options(s).filter_scheme = v;
          },
          "cutoff | fraction")
      .def_prop_rw(
          "prefactor_filter_fraction",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::prefactor_options(s)
                .filter_fraction;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::prefactor_options(s).filter_fraction =
                v;
          })
      .def_prop_rw(
          "prefactor_all_free_atoms",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::prefactor_options(s)
                .all_free_atoms;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::prefactor_options(s).all_free_atoms = v;
          })

      // --- Dynamics / MC / BH ---
      .def_prop_rw(
          "dynamics_steps",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::dynamics_options(s).steps);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::dynamics_options(s).steps = v;
          })
      .def_prop_rw(
          "dynamics_time_step",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::dynamics_options(s)
                .time_step_input;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::dynamics_options(s).time_step_input = v;
            eonc::ParametersLoadAccess::dynamics_options(s).time_step =
                v / eonc::ParametersLoadAccess::constants(s).timeUnit;
          })
      .def_prop_rw(
          "monte_carlo_steps",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::monte_carlo_options(s).steps;
          },
          [](eonc::Parameters &s, int v) {
            eonc::ParametersLoadAccess::monte_carlo_options(s).steps = v;
          })
      .def_prop_rw(
          "monte_carlo_step_size",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::monte_carlo_options(s).step_size;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::monte_carlo_options(s).step_size = v;
          })
      .def_prop_rw(
          "basin_hopping_steps",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                eonc::ParametersLoadAccess::basin_hopping_options(s).steps);
          },
          [](eonc::Parameters &s, long v) {
            eonc::ParametersLoadAccess::basin_hopping_options(s).steps = v;
          })
      .def_prop_rw(
          "basin_hopping_displacement",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::basin_hopping_options(s)
                .displacement;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::basin_hopping_options(s).displacement =
                v;
          })
      .def_prop_rw(
          "process_search_minimize_first",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::process_search_options(s)
                .minimize_first;
          },
          [](eonc::Parameters &s, bool v) {
            eonc::ParametersLoadAccess::process_search_options(s)
                .minimize_first = v;
          })

      // --- Structure comparison ---

      .def_prop_rw(
          "comp_eps_r",
          [](const eonc::Parameters &s) {
            return eonc::ParametersLoadAccess::structure_comparison_options(s)
                .distance_difference;
          },
          [](eonc::Parameters &s, double v) {
            eonc::ParametersLoadAccess::structure_comparison_options(s)
                .distance_difference = v;
          })
      .def("__repr__", [](const eonc::Parameters &self) {
        return "<Parameters job=" +
               std::string(magic_enum::enum_name(self.main_options().job)) +
               " pot=" +
               std::string(
                   magic_enum::enum_name(self.potential_options().potential)) +
               ">";
      });
}

} // namespace eonc::pybind

/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "Parameters.h"
#include "BaseStructures.h"
#include "BondBoost.h"
#include "Dynamics.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "Hessian.h"
#include "INIFile.h"
#include "ImprovedDimer.h"
#include "Job.h"
#include "NudgedElasticBand.h"
#include "Potential.h"
#include "Prefactor.h"
#include "PrefactorJob.h"
#include "ReplicaExchangeJob.h"
#include "magic_enum/magic_enum.hpp"
#include <errno.h>
#include <float.h>
#include <sstream>
#include <stdexcept>
#include <time.h>

Parameters::Parameters() {

  constants.kB = 8.6173324e-5;     // eV/K
  constants.timeUnit = 10.1805055; // fs

  // [Main] //
  main_options.job = JobType::Process_Search;
  main_options.randomSeed = -1;
  main_options.temperature = 300.0;
  main_options.checkpoint = false;
  main_options.quiet = false;
  main_options.writeLog = true;
  main_options.iniFilename = "config.ini";
  main_options.conFilename = "pos.con";
  main_options.finiteDifference = 0.01;
  main_options.maxForceCalls = 0;
  main_options.removeNetForce = true;

  // [Prefactor] //
  prefactor_options.default_value = 0.0;
  prefactor_options.max_value = 1e+21;
  prefactor_options.min_value = 1e+9;
  prefactor_options.within_radius = 3.3;
  prefactor_options.min_displacement = 0.25;
  prefactor_options.rate = Prefactor::RATE_HTST;
  prefactor_options.configuration = PrefactorJob::PREFACTOR_REACTANT;
  prefactor_options.all_free_atoms = false;
  prefactor_options.filter_scheme = Prefactor::FILTER_FRACTION;
  prefactor_options.filter_fraction = 0.90;

  // [Potential] //
  potential_options.potential = PotType::LJ;
  potential_options.MPIPollPeriod = 0.25; // seconds
  potential_options.MPIPotentialRank = -1;
  potential_options.LogPotential = false;
  potential_options.LAMMPSLogging = false;
  potential_options.LAMMPSThreads = 0;
  potential_options.EMTRasmussen = false;
  potential_options.extPotPath = "./ext_pot";

  // [AMS] //
  ams_options.engine = "";
  ams_options.forcefield = "";
  ams_options.model = "";
  ams_options.xc = "";
  ams_options.basis = "";
  ams_options.resources = "";

  // [AMS_ENV] //
  ams_options.env.amshome = "";
  ams_options.env.scm_tmpdir = "";
  ams_options.env.scm_pythondir = "";
  ams_options.env.amsbin = "";
  ams_options.env.scmlicense = "";
  ams_options.env.amsresources = "";

  // [XTBPot] //
  xtb_options.paramset = "GFNFF";
  xtb_options.acc = 1.0;
  xtb_options.elec_temperature = 0.0;
  xtb_options.maxiter = 250;
  xtb_options.charge = 0.0;
  xtb_options.uhf = 0;

  // [ZBLPot] //
  // NOTE(rg): No good defaults TBH
  zbl_options.cut_inner = 2.0;
  zbl_options.cut_global = 2.5;

  // [SocketNWChemPot] //
  socket_nwchem_options.host = "127.0.0.1";
  socket_nwchem_options.port = 9999;
  socket_nwchem_options.mem_in_gb = 2;
  socket_nwchem_options.nwchem_settings = "nwchem_settings.nwi";
  // expands to /tmp/ipi_eon_nwchem as per spec
  socket_nwchem_options.unix_socket_path = "eon_nwchem";
  socket_nwchem_options.unix_socket_mode = false;
  socket_nwchem_options.make_template_input = true;

  // [Structure Comparison] //
  structure_comparison_options.distance_difference = 0.1;
  structure_comparison_options.neighbor_cutoff = 3.3;
  structure_comparison_options.check_rotation = false;
  structure_comparison_options.indistinguishable_atoms = true;
  structure_comparison_options.energy_difference = 0.01;
  structure_comparison_options.remove_translation = true;

  // [Debug] //
  debug_options.write_movies = false;
  debug_options.write_movies_interval = 1;
  debug_options.estimate_neb_eigenvalues = false;
  debug_options.neb_mmf = LowestEigenmode::MINMODE_DIMER;

  // [Saddle Search] //
  saddle_search_options.displace_type = EpiCenters::DISP_LOAD;
  saddle_search_options.method = "min_mode";
  saddle_search_options.minmode_method = LowestEigenmode::MINMODE_DIMER;
  saddle_search_options.max_energy = 20.0;
  saddle_search_options.max_iterations = 1000;
  saddle_search_options.displace_radius = 4.0;
  saddle_search_options.displace_magnitude = 0.1;
  saddle_search_options.max_single_displace = 10.;
  saddle_search_options.nonnegative_displacement_abort = false;
  saddle_search_options.nonlocal_count_abort = 0;
  saddle_search_options.nonlocal_distance_abort = 0.0;
  saddle_search_options.remove_rotation = false;
  saddle_search_options.perp_force_ratio = 0.0;
  saddle_search_options.confine_positive.enabled = false;
  saddle_search_options.confine_positive.bowl_breakout = false;
  saddle_search_options.confine_positive.bowl_active = 20;
  saddle_search_options.confine_positive.min_force = 0.5;
  saddle_search_options.confine_positive.scale_ratio = 0.9;
  saddle_search_options.confine_positive.boost = 10.;
  saddle_search_options.confine_positive.min_active = 30;
  saddle_search_options.dynamics.temperature = 0.0;
  saddle_search_options.dynamics.state_check_interval_input = 100.0;
  saddle_search_options.dynamics.state_check_interval =
      saddle_search_options.dynamics.state_check_interval_input /
      constants.timeUnit;
  saddle_search_options.dynamics.record_interval_input = 10.0;
  saddle_search_options.dynamics.linear_interpolation = true;
  saddle_search_options.dynamics.max_init_curvature = 0.0;
  saddle_search_options.zero_mode_abort_curvature = 0.0;

  // [Optimizers] //
  optimizer_options.method = OptType::CG;
  optimizer_options.convergence_metric = "norm"s;
  optimizer_options.refine.method = OptType::None;
  optimizer_options.refine.threshold = 0.5;
  optimizer_options.max_iterations = 1000;
  optimizer_options.converged_force = 0.01;
  optimizer_options.max_move = 0.2;
  optimizer_options.time_step_input = 1.0;
  optimizer_options.max_time_step_input = 2.5;
  optimizer_options.max_time_step =
      optimizer_options.max_time_step_input / constants.timeUnit;

  optimizer_options.lbfgs.memory = 20;
  optimizer_options.lbfgs.inverse_curvature =
      0.01; // assumes stiffest curvature at minimum is 100 eV/A^2
  optimizer_options.lbfgs.auto_scale = true;
  optimizer_options.lbfgs.angle_reset = true;
  optimizer_options.lbfgs.distance_reset = true;

  optimizer_options.quickmin.steepest_descent = false;
  optimizer_options.cg.no_overshooting = false;
  optimizer_options.cg.knock_out_max_move = false;
  optimizer_options.cg.line_converged = 0.1;
  optimizer_options.cg.line_search = false;
  optimizer_options.cg.max_iter_before_reset = 0;
  optimizer_options.cg.line_search_max_iter = 10;
  optimizer_options.sd.alpha = 0.1;
  optimizer_options.sd.two_point = false;

  // [Process Search] //
  process_search_options.minimize_first = true;
  process_search_options.minimization_offset = optimizer_options.max_move;

  // [Dimer] //
  dimer_options.rotation_angle = 0.005;
  dimer_options.improved = true;
  dimer_options.converged_angle = 5.0;
  dimer_options.max_iterations = 1000;
  dimer_options.opt_method = ImprovedDimer::OPT_CG;
  dimer_options.torque_min = 0.1;
  dimer_options.torque_max = 1.0;
  dimer_options.rotations_min = 1;
  dimer_options.rotations_max = 10;
  dimer_options.remove_rotation = false;

  // [ASE_ORCA] //
  ase_orca_options.path = ""s;
  ase_orca_options.nproc = "1"s;

  // [ASE_NWCHEM] //
  ase_nwchem_options.path = ""s;
  ase_nwchem_options.nproc = "1"s;
  ase_nwchem_options.multiplicity = ""s;
  ase_nwchem_options.scf_thresh = 1e-5;
  ase_nwchem_options.scf_maxiter = 200;

  // [Metatomic] //
  metatomic_options.model_path = ""s;
  metatomic_options.device = "cpu"s;
  metatomic_options.length_unit = "angstrom"s;
  metatomic_options.extensions_directory = ""s;
  metatomic_options.check_consistency = false;
  metatomic_options.uncertainty_threshold = -1.0;
  metatomic_options.variant.base = ""s;
  metatomic_options.variant.energy = ""s;
  metatomic_options.variant.energy_uncertainty = ""s;

  // [Lanczos] //
  lanczos_options.tolerance = 0.01;
  lanczos_options.max_iterations = 20;
  lanczos_options.quit_early = true;

  // [GPR Dimer] //
  gpr_dimer_options.rotation_angle = 0.005;
  gpr_dimer_options.converged_angle = 0.08;
  gpr_dimer_options.relax_conv_angle = 0.001;
  gpr_dimer_options.init_rotations_max = 6;
  gpr_dimer_options.relax_rotations_max = 10;
  gpr_dimer_options.divisor_t_dimer_gp = 10;
  gpr_dimer_options.max_outer_iterations = 300;
  gpr_dimer_options.max_inner_iterations = 1000;
  gpr_dimer_options.midpoint_max_disp = 0.5;
  gpr_dimer_options.rot_opt_method = "lbfgs";
  gpr_dimer_options.trans_opt_method = "lbfgs";
  gpr_dimer_options.active_radius = 5.0;
  gpr_dimer_options.dimer_sep = 0.01;
  gpr_dimer_options.conv_step = 0.1;
  gpr_dimer_options.max_step = 0.1;
  gpr_dimer_options.ratio_at_limit = 0.66667;
  gpr_dimer_options.init_rot_gp = 0;
  gpr_dimer_options.init_trans_gp = 0;
  gpr_dimer_options.many_iterations = true;
  // GPR Params
  gpr_dimer_options.gpr_params.hyper_opt_method = "scg";
  gpr_dimer_options.gpr_params.sigma2 = 1e-8;
  gpr_dimer_options.gpr_params.jitter_sigma2 = 0;
  gpr_dimer_options.gpr_params.noise_sigma2 = 1e-8;
  gpr_dimer_options.gpr_params.prior_mu = 0;
  gpr_dimer_options.gpr_params.prior_sigma2 = 1;
  gpr_dimer_options.gpr_params.prior_nu = 20;
  // GPR Optimization Parameters
  gpr_dimer_options.opt_params.check_derivatives = false;
  gpr_dimer_options.opt_params.max_iterations = 400;
  gpr_dimer_options.opt_params.tol_func = 1e-4;
  gpr_dimer_options.opt_params.tol_sol = 1e-4;
  gpr_dimer_options.opt_params.lambda_limit = 1e17;
  gpr_dimer_options.opt_params.lambda_init = 10.0;
  // GPR Prune Parameters
  gpr_dimer_options.prune_params.use_prune = false;
  gpr_dimer_options.prune_params.begin = 8;
  gpr_dimer_options.prune_params.n_vals = 3;
  gpr_dimer_options.prune_params.threshold = 0.5;
  // GPR Debugging Parameters
  gpr_dimer_options.debug_params.report_level = 1;
  gpr_dimer_options.debug_params.debug_level = 2;
  gpr_dimer_options.debug_params.out_dir = "output";
  gpr_dimer_options.debug_params.pos_file = "position";
  gpr_dimer_options.debug_params.energy_file = "energy";
  gpr_dimer_options.debug_params.grad_file = "gradient";
  gpr_dimer_options.debug_params.out_ext = "dat";
  gpr_dimer_options.debug_params.offset_mid_point = 3.;
  gpr_dimer_options.debug_params.dy = 0.1;
  gpr_dimer_options.debug_params.dz = 0.1;

  // GP Surrogate Parameters
  gp_surrogate_options.enabled = false;
  gp_surrogate_options.sub_job = JobType::Unknown;
  gp_surrogate_options.uncertainty = 0.05;
  gp_surrogate_options.linear_path_always = false;
  gp_surrogate_options.potential = PotType::CatLearn;

  // [Hessian] //
  hessian_options.atom_list = string("All");
  hessian_options.zero_freq_value = 1e-6;

  // [Nudged Elastic Band] //
  // General optimization and path parameters
  neb_options.image_count = 5;
  neb_options.max_iterations = 1000;
  neb_options.opt_method = OptType::LBFGS;
  neb_options.force_tolerance = optimizer_options.converged_force;
  // Post-run peak handling
  neb_options.mmf_peaks.enabled = true;
  neb_options.mmf_peaks.tolerance = 0.05;

  // Inter-image spring dynamics
  neb_options.spring.constant = 5.0;
  neb_options.spring.use_elastic_band = false;
  neb_options.spring.doubly_nudged = false;
  neb_options.spring.use_switching = false;

  // Energy-weighted spring adjustments for high-curvature regions
  neb_options.spring.weighting.enabled = false;
  neb_options.spring.weighting.trigger = 10.0;
  neb_options.spring.weighting.k_max = 9.7;
  neb_options.spring.weighting.k_min = 0.97;

  // Onsager-Machlup settings
  neb_options.spring.om.enabled = false;
  neb_options.spring.om.optimize_k = true;
  neb_options.spring.om.k_scale = 1.0;
  neb_options.spring.om.k_min = 0.1;
  neb_options.spring.om.k_max = 100;

  // Climbing Image (CI-NEB) settings for barrier identification
  neb_options.climbing_image.enabled = true;
  neb_options.climbing_image.converged_only = true;
  neb_options.climbing_image.use_old_tangent = false;
  neb_options.climbing_image.trigger_force =
      std::numeric_limits<double>::infinity();
  neb_options.climbing_image.trigger_factor = 0.0;

  // Hybrid NEB-Dimer (RONEB) parameters using Min-Mode Following (MMF)
  neb_options.climbing_image.roneb.use_mmf = false;
  neb_options.climbing_image.roneb.trigger_force = 0.1;
  neb_options.climbing_image.roneb.trigger_factor = 0.0;
  // Use the angle criteria instead
  neb_options.climbing_image.roneb.max_steps = 1000;
  neb_options.climbing_image.roneb.ci_stability_count = 5;
  neb_options.climbing_image.roneb.angle_tol = 0.8;
  neb_options.climbing_image.roneb.penalty.base = 0.1;
  neb_options.climbing_image.roneb.penalty.strength = 0.5;

  // Initial path guess and pre-optimization
  neb_options.initialization.method = NEBInit::LINEAR;
  neb_options.initialization.input_path = ""s;
  neb_options.initialization.max_iterations = 5000;
  neb_options.initialization.nsteps = 250;
  neb_options.initialization.max_move = 0.1;
  neb_options.initialization.force_tolerance = 0.001;
  neb_options.initialization.sidpp_alpha = 0.33;
  neb_options.initialization.opt_method = OptType::LBFGS;
  neb_options.initialization.oversampling = false;
  neb_options.initialization.oversampling_factor = 3;

  // Boundary condition preferences
  neb_options.endpoints.minimize = true;
  neb_options.endpoints.use_path_file = false;

  // [Dynamics] //
  dynamics_options.time_step_input = 1.0;
  dynamics_options.time_input = 1000.0;
  dynamics_options.time_step =
      dynamics_options.time_step_input / constants.timeUnit;
  dynamics_options.time = dynamics_options.time_input / constants.timeUnit;
  dynamics_options.steps =
      long(floor(dynamics_options.time / dynamics_options.time_step + 0.5));

  // [Thermostat] //
  thermostat_options.kind = Dynamics::NONE;
  thermostat_options.andersen_alpha = 1.0;
  thermostat_options.andersen_tcol_input = 100.0;
  thermostat_options.nose_mass = 1.0;
  thermostat_options.langevin_friction_input = 0.01;
  thermostat_options.langevin_friction =
      thermostat_options.langevin_friction_input * constants.timeUnit;

  // [Parallel Replica] //
  parallel_replica_options.refine_transition = true;
  parallel_replica_options.auto_stop = false;
  parallel_replica_options.dephase_loop_stop = false;
  parallel_replica_options.dephase_time_input = 1000.0;
  parallel_replica_options.dephase_loop_max = 5;
  parallel_replica_options.state_check_interval_input = 1000.0;
  parallel_replica_options.record_interval_input = 50.0;
  parallel_replica_options.corr_time_input = 1000.0;
  parallel_replica_options.dephase_time =
      parallel_replica_options.dephase_time_input / constants.timeUnit;
  parallel_replica_options.state_check_interval =
      parallel_replica_options.state_check_interval_input / constants.timeUnit;
  parallel_replica_options.record_interval =
      parallel_replica_options.record_interval_input / constants.timeUnit;
  parallel_replica_options.corr_time =
      parallel_replica_options.corr_time_input / constants.timeUnit;

  // [Temperature Accelerated Dynamics] //
  tad_options.low_temperature = 300.0;
  tad_options.min_prefactor = 0.001;
  tad_options.confidence = 0.001;

  // [Replica Exchange] //
  replica_exchange_options.temperature_distribution = "exponential";
  replica_exchange_options.replicas = 10;
  replica_exchange_options.exchange_trials = replica_exchange_options.replicas;
  replica_exchange_options.sampling_time_input = 1000.0;
  replica_exchange_options.sampling_time =
      replica_exchange_options.sampling_time_input / constants.timeUnit;
  replica_exchange_options.temperature_low = 0.0;
  replica_exchange_options.temperature_high = 0.0;
  replica_exchange_options.exchange_period = 100.0;

  // [Hyperdynamics] //
  hyperdynamics_options.bias_potential = Hyperdynamics::NONE;
  hyperdynamics_options.boost_atom_list = string("All");
  hyperdynamics_options.dvmax = 0.0;
  hyperdynamics_options.qrr = 0.2;
  hyperdynamics_options.prr = 0.95;
  hyperdynamics_options.qcut = 3.0;
  hyperdynamics_options.rmd_time_input = 100.0;
  hyperdynamics_options.rmd_time =
      hyperdynamics_options.rmd_time_input / constants.timeUnit;

  // [Basin Hopping] //
  basin_hopping_options.displacement = 0.5;
  basin_hopping_options.push_apart_distance = 0.4;
  basin_hopping_options.initial_random_structure_probability = 0.0;
  basin_hopping_options.steps = 10000;
  basin_hopping_options.quenching_steps = 0;
  basin_hopping_options.single_atom_displace = false;
  basin_hopping_options.significant_structure = true;
  basin_hopping_options.displacement_algorithm = "standard";
  basin_hopping_options.displacement_distribution = "uniform";
  basin_hopping_options.swap_probability = 0.0;
  basin_hopping_options.jump_max = 10;
  basin_hopping_options.jump_steps = 0;
  basin_hopping_options.adjust_displacement = true;
  basin_hopping_options.adjust_period = 10;
  basin_hopping_options.adjust_fraction = 0.05;
  basin_hopping_options.target_ratio = 0.5;
  basin_hopping_options.write_unique = false;
  basin_hopping_options.stop_energy = -DBL_MAX;

  // [Global Optimization] //
  global_optimization_options.move_method = "md";
  global_optimization_options.decision_method = "npew";
  global_optimization_options.steps = 10000;
  global_optimization_options.beta = 1.05;
  global_optimization_options.alpha = 1.02;
  global_optimization_options.mdmin = 3;
  global_optimization_options.target_energy = -1.E50;

  // [Monte Carlo] //
  monte_carlo_options.step_size = 0.005;
  monte_carlo_options.steps = 1000;

  // [BGSD] //
  bgsd_options.alpha = 10.0;
  bgsd_options.beta = 0.2;
  bgsd_options.gradient_finite_difference = 0.000001;
  bgsd_options.h_force_convergence = 0.01;
  bgsd_options.grad2energy_convergence = 0.000001;
  bgsd_options.grad2force_convergence = 0.0001;

  // [Serve] //
  serve_options.host = "localhost";
  serve_options.port = 12345;
  serve_options.replicas = 1;
  serve_options.gateway_port = 0;
  serve_options.endpoints = "";

  // [CatLearn] //
  // No reasonable default for catlearn_options.path
  catlearn_options.model = "gp";
  catlearn_options.prior = "median";
  catlearn_options.use_deriv = true;
  catlearn_options.use_fingerprint = false;
  catlearn_options.parallel = false;
}

string Parameters::toLowerCase(string s) {
  for (string::size_type i = 0; i < s.length(); ++i) {
    s[i] = tolower(s[i]);
  }
  return s;
}

int Parameters::load(string filename) {
  FILE *fh;

  fh = fopen(filename.c_str(), "rb");
  if (fh == NULL) {
    fprintf(stderr, "error: %s\n", strerror(errno));
    return 1;
  }
  int error = load(fh);
  fclose(fh);
  return error;
}

int Parameters::load(FILE *file) {

  CIniFile ini;
  ini.CaseInsensitive();
  int error = 0;

  if (ini.ReadFile(file)) {

    // [Main] //

    main_options.job =
        magic_enum::enum_cast<JobType>(ini.GetValue("Main", "job"),
                                       magic_enum::case_insensitive)
            .value_or(JobType::Unknown);
    main_options.temperature =
        ini.GetValueF("Main", "temperature", main_options.temperature);
    main_options.randomSeed =
        ini.GetValueL("Main", "random_seed", main_options.randomSeed);
    main_options.checkpoint =
        ini.GetValueB("Main", "checkpoint", main_options.checkpoint);
    main_options.quiet = ini.GetValueB("Main", "quiet", main_options.quiet);
    main_options.writeLog =
        ini.GetValueB("Main", "write_log", main_options.writeLog);
    main_options.finiteDifference = ini.GetValueF(
        "Main", "finite_difference", main_options.finiteDifference);
    // Initialize random generator
    if (main_options.randomSeed < 0) {
      unsigned i = time(NULL);
      main_options.randomSeed = i;
      helper_functions::random(i);
    } else {
      helper_functions::random(main_options.randomSeed);
    }
    main_options.maxForceCalls =
        ini.GetValueL("Main", "max_force_calls", main_options.maxForceCalls);
    main_options.removeNetForce =
        ini.GetValueB("Main", "remove_net_force", main_options.removeNetForce);

    // [Potential] //

    potential_options.potential =
        magic_enum::enum_cast<PotType>(ini.GetValue("Potential", "potential"),
                                       magic_enum::case_insensitive)
            .value_or(PotType::UNKNOWN);
    potential_options.MPIPollPeriod = ini.GetValueF(
        "Potential", "mpi_poll_period", potential_options.MPIPollPeriod);
    potential_options.LAMMPSLogging = ini.GetValueB(
        "Potential", "lammps_logging", potential_options.LAMMPSLogging);
    potential_options.LAMMPSThreads = (int)ini.GetValueL(
        "Potential", "lammps_threads", potential_options.LAMMPSThreads);
    potential_options.EMTRasmussen = ini.GetValueB(
        "Potential", "emt_rasmussen", potential_options.EMTRasmussen);
    potential_options.extPotPath =
        ini.GetValue("Potential", "ext_pot_path", potential_options.extPotPath);

    if (potential_options.potential == PotType::MPI ||
        potential_options.potential == PotType::VASP ||
        potential_options.potential == PotType::BOPFOX ||
        potential_options.potential == PotType::BOP) {
      potential_options.LogPotential = true;
    } else {
      potential_options.LogPotential = false;
    }
    potential_options.LogPotential = ini.GetValueB(
        "Potential", "log_potential", potential_options.LogPotential);

    // [AMS]
    if (potential_options.potential == PotType::AMS) {
      ams_options.engine = ini.GetValue("AMS", "engine", ams_options.engine);
      ams_options.forcefield =
          ini.GetValue("AMS", "forcefield", ams_options.forcefield);
      ams_options.resources =
          ini.GetValue("AMS", "resources", ams_options.resources);
      ams_options.model = ini.GetValue("AMS", "model", ams_options.model);
      ams_options.xc = ini.GetValue("AMS", "xc", ams_options.xc);
      ams_options.basis = ini.GetValue("AMS", "basis", ams_options.basis);
    }
    // [AMS_IO]
    if (potential_options.potential == PotType::AMS_IO) {
      ams_options.engine = ini.GetValue("AMS_IO", "engine", ams_options.engine);
      ams_options.forcefield =
          ini.GetValue("AMS_IO", "forcefield", ams_options.forcefield);
      ams_options.model = ini.GetValue("AMS_IO", "model", ams_options.model);
      ams_options.xc = ini.GetValue("AMS_IO", "xc", ams_options.xc);
    }
    // [AMS_ENV]
    // This is only needed if the regular calls do not work
    // e.g. on a MacOS machine
    if (potential_options.potential == PotType::AMS_IO ||
        potential_options.potential == PotType::AMS) {
      ams_options.env.amshome =
          ini.GetValue("AMS_ENV", "amshome", ams_options.env.amshome);
      ams_options.env.scm_tmpdir =
          ini.GetValue("AMS_ENV", "scm_tmpdir", ams_options.env.scm_tmpdir);
      ams_options.env.scmlicense =
          ini.GetValue("AMS_ENV", "scmlicense", ams_options.env.scmlicense);
      ams_options.env.scm_pythondir = ini.GetValue(
          "AMS_ENV", "scm_pythondir", ams_options.env.scm_pythondir);
      ams_options.env.amsbin =
          ini.GetValue("AMS_ENV", "amsbin", ams_options.env.amsbin);
      ams_options.env.amsresources =
          ini.GetValue("AMS_ENV", "amsresources", ams_options.env.amsresources);
    }
    // [XTBPot]
    if (potential_options.potential == PotType::XTB) {
      xtb_options.paramset =
          ini.GetValue("XTBPot", "paramset", xtb_options.paramset);
      xtb_options.acc = ini.GetValueF("XTBPot", "accuracy", xtb_options.acc);
      xtb_options.elec_temperature = ini.GetValueF(
          "XTBPot", "electronic_temperature", xtb_options.elec_temperature);
      xtb_options.maxiter =
          ini.GetValueL("XTBPot", "max_iterations", xtb_options.maxiter);
      xtb_options.uhf = ini.GetValueL("XTBPot", "uhf", xtb_options.uhf);
      xtb_options.charge =
          ini.GetValueF("XTBPot", "charge", xtb_options.charge);
    }
    // [ZBLPot]
    if (potential_options.potential == PotType::ZBL) {
      zbl_options.cut_inner =
          ini.GetValueF("ZBLPot", "cut_inner", zbl_options.cut_inner);
      zbl_options.cut_global =
          ini.GetValueF("ZBLPot", "cut_global", zbl_options.cut_global);
      if (zbl_options.cut_inner > zbl_options.cut_global) {
        throw std::runtime_error(
            "Switching function must begin before the global cutoff!");
      }
    }
    // [SocketNWChemPot]
    if (potential_options.potential == PotType::SocketNWChem) {
      socket_nwchem_options.host =
          ini.GetValue("SocketNWChemPot", "host", socket_nwchem_options.host);
      socket_nwchem_options.port =
          ini.GetValueL("SocketNWChemPot", "port", socket_nwchem_options.port);
      socket_nwchem_options.mem_in_gb = ini.GetValueL(
          "SocketNWChemPot", "mem_in_gb", socket_nwchem_options.mem_in_gb);
      socket_nwchem_options.nwchem_settings =
          ini.GetValue("SocketNWChemPot", "nwchem_settings",
                       socket_nwchem_options.nwchem_settings);
      socket_nwchem_options.unix_socket_path =
          ini.GetValue("SocketNWChemPot", "unix_socket_path",
                       socket_nwchem_options.unix_socket_path);
      socket_nwchem_options.unix_socket_mode =
          ini.GetValueB("SocketNWChemPot", "unix_socket_mode",
                        socket_nwchem_options.unix_socket_mode);
      socket_nwchem_options.make_template_input =
          ini.GetValueB("SocketNWChemPot", "make_template_input",
                        socket_nwchem_options.make_template_input);
    }

    // [Debug] //

    debug_options.write_movies =
        ini.GetValueB("Debug", "write_movies", debug_options.write_movies);
    debug_options.write_movies_interval = ini.GetValueL(
        "Debug", "write_movies_interval", debug_options.write_movies_interval);
    debug_options.estimate_neb_eigenvalues =
        ini.GetValueB("Debug", "estimate_neb_eigenvalues",
                      debug_options.estimate_neb_eigenvalues);
    debug_options.neb_mmf = toLowerCase(
        ini.GetValue("Debug", "neb_mmf_estimator", debug_options.neb_mmf));

    // [Structure Comparison] //

    structure_comparison_options.distance_difference =
        ini.GetValueF("Structure Comparison", "distance_difference",
                      structure_comparison_options.distance_difference);
    structure_comparison_options.neighbor_cutoff =
        ini.GetValueF("Structure Comparison", "neighbor_cutoff",
                      structure_comparison_options.neighbor_cutoff);
    structure_comparison_options.check_rotation =
        ini.GetValueB("Structure Comparison", "check_rotation",
                      structure_comparison_options.check_rotation);
    structure_comparison_options.energy_difference =
        ini.GetValueF("Structure Comparison", "energy_difference",
                      structure_comparison_options.energy_difference);
    structure_comparison_options.indistinguishable_atoms =
        ini.GetValueB("Structure Comparison", "indistinguishable_atoms",
                      structure_comparison_options.indistinguishable_atoms);
    structure_comparison_options.remove_translation =
        ini.GetValueB("Structure Comparison", "remove_translation",
                      structure_comparison_options.remove_translation);

    // [Process Search] //

    process_search_options.minimize_first =
        ini.GetValueB("Process Search", "minimize_first",
                      process_search_options.minimize_first);
    process_search_options.minimization_offset =
        ini.GetValueF("Process Search", "minimization_offset",
                      process_search_options.minimization_offset);

    // [Optimizers] //
    auto inp_optMethod = magic_enum::enum_cast<OptType>(
                             ini.GetValue("Optimizer", "opt_method", "none"),
                             magic_enum::case_insensitive)
                             .value_or(OptType::Unknown);
    if (inp_optMethod != OptType::None) {
      optimizer_options.method = inp_optMethod;
    }

    optimizer_options.convergence_metric =
        toLowerCase(ini.GetValue("Optimizer", "convergence_metric",
                                 optimizer_options.convergence_metric));
    if (optimizer_options.convergence_metric == "max_atom") {
      optimizer_options.convergence_metric_label = "Max atom force";
    } else if (optimizer_options.convergence_metric == "max_component") {
      optimizer_options.convergence_metric_label = "Max force comp";
    } else if (optimizer_options.convergence_metric == "norm") {
      // optimizer_options.convergence_metric_label = "\u2016Force\u2016";
      optimizer_options.convergence_metric_label = "||Force||";
    } else {
      fprintf(stderr, "unknown convergence_metric %s\n",
              optimizer_options.convergence_metric.c_str());
      exit(1);
    }

    if (ini.FindKey("Refine") != -1) {
      optimizer_options.refine.method =
          magic_enum::enum_cast<OptType>(ini.GetValue("Refine", "opt_method"),
                                         magic_enum::case_insensitive)
              .value_or(OptType::None);
      optimizer_options.refine.threshold = ini.GetValueF(
          "Refine", "threshold", optimizer_options.refine.threshold);
    }

    optimizer_options.converged_force = ini.GetValueF(
        "Optimizer", "converged_force", optimizer_options.converged_force);
    optimizer_options.max_iterations = static_cast<size_t>(ini.GetValueL(
        "Optimizer", "max_iterations", optimizer_options.max_iterations));
    optimizer_options.max_move =
        ini.GetValueF("Optimizer", "max_move", optimizer_options.max_move);
    process_search_options.minimization_offset = optimizer_options.max_move;
    // Handle each optimizer separately
    if (ini.FindKey("QuickMin") != -1) {
      optimizer_options.time_step_input = ini.GetValueF(
          "QuickMin", "time_step", optimizer_options.time_step_input);
      optimizer_options.time_step =
          optimizer_options.time_step_input / constants.timeUnit;
      optimizer_options.quickmin.steepest_descent =
          ini.GetValueB("Optimizer", "qm_steepest_descent",
                        optimizer_options.quickmin.steepest_descent);
    }
    if (ini.FindKey("FIRE") != -1) {
      SPDLOG_WARN("Overwriting QuickMin timestep with Fire timestep!!");
      optimizer_options.time_step_input =
          ini.GetValueF("FIRE", "time_step", optimizer_options.time_step_input);
      optimizer_options.time_step =
          optimizer_options.time_step_input / constants.timeUnit;
      optimizer_options.max_time_step_input = ini.GetValueF(
          "FIRE", "time_step_max", optimizer_options.max_time_step_input);
      optimizer_options.max_time_step =
          optimizer_options.max_time_step_input / constants.timeUnit;
    }
    if (ini.FindKey("LBFGS") != -1) {
      optimizer_options.lbfgs.memory = ini.GetValueL(
          "LBFGS", "lbfgs_memory", optimizer_options.lbfgs.memory);
      optimizer_options.lbfgs.inverse_curvature =
          ini.GetValueF("LBFGS", "lbfgs_inverse_curvature",
                        optimizer_options.lbfgs.inverse_curvature);
      optimizer_options.lbfgs.max_inverse_curvature =
          ini.GetValueF("LBFGS", "lbfgs_max_inverse_curvature",
                        optimizer_options.lbfgs.max_inverse_curvature);
      optimizer_options.lbfgs.auto_scale = ini.GetValueB(
          "LBFGS", "lbfgs_auto_scale", optimizer_options.lbfgs.auto_scale);
      optimizer_options.lbfgs.angle_reset = ini.GetValueB(
          "LBFGS", "lbfgs_angle_reset", optimizer_options.lbfgs.angle_reset);
      optimizer_options.lbfgs.distance_reset =
          ini.GetValueB("LBFGS", "lbfgs_distance_reset",
                        optimizer_options.lbfgs.distance_reset);
    }
    if (ini.FindKey("CG") != -1) {
      optimizer_options.cg.no_overshooting = ini.GetValueB(
          "CG", "cg_no_overshooting", optimizer_options.cg.no_overshooting);
      optimizer_options.cg.knock_out_max_move =
          ini.GetValueB("CG", "cg_knock_out_max_move",
                        optimizer_options.cg.knock_out_max_move);
      optimizer_options.cg.line_search = ini.GetValueB(
          "CG", "cg_line_search", optimizer_options.cg.line_search);
      optimizer_options.cg.line_converged = ini.GetValueF(
          "CG", "cg_line_converged", optimizer_options.cg.line_converged);
      optimizer_options.cg.max_iter_before_reset =
          ini.GetValueL("CG", "cg_max_iter_before_reset",
                        optimizer_options.cg.max_iter_before_reset);
      optimizer_options.cg.line_search_max_iter =
          ini.GetValueL("CG", "cg_max_iter_line_search",
                        optimizer_options.cg.line_search_max_iter);
    }
    if (ini.FindKey("SD") != -1) {
      optimizer_options.sd.alpha =
          ini.GetValueF("SD", "sd_alpha", optimizer_options.sd.alpha);
      optimizer_options.sd.two_point =
          ini.GetValueB("SD", "sd_twopoint", optimizer_options.sd.two_point);
    }

    // [Dimer] //

    dimer_options.rotation_angle =
        ini.GetValueF("Dimer", "finite_angle", dimer_options.rotation_angle);
    dimer_options.improved =
        ini.GetValueB("Dimer", "improved", dimer_options.improved);
    dimer_options.converged_angle = ini.GetValueF(
        "Dimer", "converged_angle", dimer_options.converged_angle);
    dimer_options.max_iterations =
        ini.GetValueL("Dimer", "max_iterations", dimer_options.max_iterations);
    dimer_options.opt_method = toLowerCase(
        ini.GetValue("Dimer", "opt_method", dimer_options.opt_method));
    dimer_options.rotations_min =
        ini.GetValueL("Dimer", "rotations_min", dimer_options.rotations_min);
    dimer_options.rotations_max =
        ini.GetValueL("Dimer", "rotations_max", dimer_options.rotations_max);
    dimer_options.torque_min =
        ini.GetValueF("Dimer", "torque_min", dimer_options.torque_min);
    dimer_options.torque_max =
        ini.GetValueF("Dimer", "torque_max", dimer_options.torque_max);
    dimer_options.remove_rotation = ini.GetValueB(
        "Dimer", "remove_rotation", dimer_options.remove_rotation);

    // GP Surrogate Parameters
    gp_surrogate_options.enabled =
        ini.GetValueB("Surrogate", "use_surrogate", false);
    // If use_surrogate is true, job -> gp_surrogate and sub_job->job
    if (gp_surrogate_options.enabled) {
      // TODO: What about other jobs
      gp_surrogate_options.sub_job = main_options.job;
      main_options.job = JobType::GP_Surrogate;
    }
    gp_surrogate_options.uncertainty = ini.GetValueF(
        "Surrogate", "gp_uncertainty", gp_surrogate_options.uncertainty);
    if (ini.FindKey("Surrogate") != -1) {
      gp_surrogate_options.potential =
          magic_enum::enum_cast<PotType>(ini.GetValue("Surrogate", "potential"),
                                         magic_enum::case_insensitive)
              .value_or(PotType::UNKNOWN);
      if (gp_surrogate_options.potential != PotType::CatLearn) {
        throw std::runtime_error("We only support catlearn for GP right now");
      }
    }
    // [CatLearn]
    if (ini.FindKey("CatLearn") != -1) {
      // Case sensitive!!
      catlearn_options.path = ini.GetValue("CatLearn", "catl_path");
      catlearn_options.model = ini.GetValue("CatLearn", "model", "catl_model");
      catlearn_options.prior = ini.GetValue("CatLearn", "prior", "catl_prior");
      catlearn_options.use_deriv =
          ini.GetValueB("CatLearn", "use_derivatives", "catl_deriv");
      catlearn_options.use_fingerprint =
          ini.GetValueB("CatLearn", "use_fingerprint", "catl_fingerprint");
      catlearn_options.parallel = ini.GetValueB(
          "CatLearn", "parallel_hyperparameter_opt", "catl_parallel");
    }
    // [ASE_ORCA]
    if (ini.FindKey("ASE_ORCA") != -1) {
      // Case sensitive!!
      // TODO: This should be handled in clienteon so you can still call
      // eonclient for single point calculations easily
      ase_orca_options.path = ini.GetValue("ASE_ORCA", "orca_path", "");
      ase_orca_options.nproc = ini.GetValue("ASE_ORCA", "nproc");
      ase_orca_options.simpleinput =
          ini.GetValue("ASE_ORCA", "simpleinput", "");
    }
    // [ASE_NWCHEM]
    if (ini.FindKey("ASE_NWCHEM") != -1) {
      // Case sensitive!!
      // TODO: This should be handled in clienteon so you can still call
      // eonclient for single point calculations easily
      ase_nwchem_options.path = ini.GetValue("ASE_NWCHEM", "nwchem_path", "");
      ase_nwchem_options.nproc = ini.GetValue("ASE_NWCHEM", "nproc");
      ase_nwchem_options.multiplicity =
          ini.GetValue("ASE_NWCHEM", "multiplicity");
      ase_nwchem_options.scf_thresh =
          ini.GetValueF("ASE_NWCHEM", "scf_thresh", 1e-5);
      ase_nwchem_options.scf_maxiter =
          ini.GetValueL("ASE_NWCHEM", "scf_maxiter", 200);
    }
    // [Metatomic]
    if (ini.FindKey("Metatomic") != -1) {
      // Case sensitive!!
      // TODO: This should be handled in clienteon so you can still call
      // eonclient for single point calculations easily
      metatomic_options.model_path =
          ini.GetValue("Metatomic", "model_path", "");
      metatomic_options.device = ini.GetValue("Metatomic", "device", "cpu");
      metatomic_options.length_unit =
          ini.GetValue("Metatomic", "length_unit", "angstrom");
      metatomic_options.extensions_directory =
          ini.GetValue("Metatomic", "extensions_directory", "");
      metatomic_options.check_consistency =
          ini.GetValueB("Metatomic", "check_consistency", false);
      metatomic_options.uncertainty_threshold =
          ini.GetValueF("Metatomic", "uncertainty_threshold", -1.0);
      auto &_variant = metatomic_options.variant;
      _variant.base = ini.GetValue("Metatomic", "variant_base", "");
      _variant.energy = ini.GetValue("Metatomic", "variant_energy", "");
      _variant.energy_uncertainty =
          ini.GetValue("Metatomic", "variant_energy_uncertainty", "");
    }
    // [Serve]
    if (ini.FindKey("Serve") != -1) {
      serve_options.host = ini.GetValue("Serve", "host", serve_options.host);
      serve_options.port = static_cast<uint16_t>(
          ini.GetValueL("Serve", "port", serve_options.port));
      serve_options.replicas = static_cast<size_t>(
          ini.GetValueL("Serve", "replicas", serve_options.replicas));
      serve_options.gateway_port = static_cast<uint16_t>(
          ini.GetValueL("Serve", "gateway_port", serve_options.gateway_port));
      serve_options.endpoints =
          ini.GetValue("Serve", "endpoints", serve_options.endpoints);
    }

    // GP_NEB only
    gp_surrogate_options.linear_path_always =
        ini.GetValueB("Surrogate", "gp_linear_path_always",
                      gp_surrogate_options.linear_path_always);
    // [Lanczos] //

    lanczos_options.tolerance =
        ini.GetValueF("Lanczos", "tolerance", lanczos_options.tolerance);
    lanczos_options.max_iterations = ini.GetValueL(
        "Lanczos", "max_iterations", lanczos_options.max_iterations);
    lanczos_options.quit_early =
        ini.GetValueB("Lanczos", "quit_early", lanczos_options.quit_early);

    // [GPR Dimer] //
    gpr_dimer_options.rotation_angle = ini.GetValueF(
        "GPR Dimer", "finite_angle", gpr_dimer_options.rotation_angle);
    gpr_dimer_options.converged_angle = ini.GetValueF(
        "GPR Dimer", "converged_angle", gpr_dimer_options.converged_angle);
    gpr_dimer_options.relax_conv_angle =
        ini.GetValueF("GPR Dimer", "relaxation_converged_angle",
                      gpr_dimer_options.relax_conv_angle);
    gpr_dimer_options.init_rotations_max =
        (int)ini.GetValueL("GPR Dimer", "max_initial_rotation_iterations",
                           gpr_dimer_options.init_rotations_max);
    gpr_dimer_options.relax_rotations_max =
        (int)ini.GetValueL("GPR Dimer", "max_relaxation_rotation_iterations",
                           gpr_dimer_options.relax_rotations_max);
    gpr_dimer_options.divisor_t_dimer_gp = (int)ini.GetValueL(
        "GPR Dimer", "divisor_t_dimer", gpr_dimer_options.divisor_t_dimer_gp);
    gpr_dimer_options.max_outer_iterations =
        (int)ini.GetValueL("GPR Dimer", "max_outer_iterations",
                           gpr_dimer_options.max_outer_iterations);
    gpr_dimer_options.max_inner_iterations =
        (int)ini.GetValueL("GPR Dimer", "max_inner_iterations",
                           gpr_dimer_options.max_inner_iterations);
    gpr_dimer_options.midpoint_max_disp =
        ini.GetValueF("GPR Dimer", "max_midpoint_displacement",
                      gpr_dimer_options.midpoint_max_disp);
    gpr_dimer_options.rot_opt_method = ini.GetValue(
        "GPR Dimer", "rotation_opt_method", gpr_dimer_options.rot_opt_method);
    gpr_dimer_options.trans_opt_method =
        ini.GetValue("GPR Dimer", "translation_opt_method",
                     gpr_dimer_options.trans_opt_method);
    gpr_dimer_options.active_radius = ini.GetValueF(
        "GPR Dimer", "active_radius", gpr_dimer_options.active_radius);
    gpr_dimer_options.dimer_sep = ini.GetValueF("GPR Dimer", "dimer_separation",
                                                gpr_dimer_options.dimer_sep);
    gpr_dimer_options.conv_step = ini.GetValueF(
        "GPR Dimer", "convex_region_step_size", gpr_dimer_options.conv_step);
    gpr_dimer_options.max_step =
        ini.GetValueF("GPR Dimer", "max_step_size", gpr_dimer_options.max_step);
    gpr_dimer_options.ratio_at_limit = ini.GetValueF(
        "GPR Dimer", "ratio_at_limit", gpr_dimer_options.ratio_at_limit);
    gpr_dimer_options.init_rot_gp = ini.GetValueB(
        "GPR Dimer", "nogp_initial_rotations", gpr_dimer_options.init_rot_gp);
    gpr_dimer_options.init_trans_gp = ini.GetValueB(
        "GPR Dimer", "nogp_init_translations", gpr_dimer_options.init_trans_gp);
    gpr_dimer_options.many_iterations = ini.GetValueB(
        "GPR Dimer", "has_many_iterations", gpr_dimer_options.many_iterations);
    // GPR Params
    gpr_dimer_options.gpr_params.hyper_opt_method =
        ini.GetValue("GPR Dimer", "hyperparameter_opt_method",
                     gpr_dimer_options.gpr_params.hyper_opt_method);
    gpr_dimer_options.gpr_params.sigma2 = ini.GetValueF(
        "GPR Dimer", "gpr_variance", gpr_dimer_options.gpr_params.sigma2);
    gpr_dimer_options.gpr_params.jitter_sigma2 =
        ini.GetValueF("GPR Dimer", "gpr_jitter_variance",
                      gpr_dimer_options.gpr_params.jitter_sigma2);
    gpr_dimer_options.gpr_params.noise_sigma2 =
        ini.GetValueF("GPR Dimer", "gpr_noise_variance",
                      gpr_dimer_options.gpr_params.noise_sigma2);
    gpr_dimer_options.gpr_params.prior_mu = ini.GetValueF(
        "GPR Dimer", "prior_mean", gpr_dimer_options.gpr_params.prior_mu);
    gpr_dimer_options.gpr_params.prior_sigma2 =
        ini.GetValueF("GPR Dimer", "prior_variance",
                      gpr_dimer_options.gpr_params.prior_sigma2);
    gpr_dimer_options.gpr_params.prior_nu =
        ini.GetValueF("GPR Dimer", "prior_degrees_of_freedom",
                      gpr_dimer_options.gpr_params.prior_nu);
    // GPR Optimization Parameters
    gpr_dimer_options.opt_params.check_derivatives =
        ini.GetValueB("GPR Dimer", "check_derivatives",
                      gpr_dimer_options.opt_params.check_derivatives);
    gpr_dimer_options.opt_params.max_iterations =
        (int)ini.GetValueL("GPR Dimer", "opt_max_iterations",
                           gpr_dimer_options.opt_params.max_iterations);
    gpr_dimer_options.opt_params.tol_func = ini.GetValueF(
        "GPR Dimer", "opt_tol_func", gpr_dimer_options.opt_params.tol_func);
    gpr_dimer_options.opt_params.tol_sol = ini.GetValueF(
        "GPR Dimer", "opt_tol_sol", gpr_dimer_options.opt_params.tol_sol);
    gpr_dimer_options.opt_params.lambda_limit =
        ini.GetValueF("GPR Dimer", "opt_lambda_limit",
                      gpr_dimer_options.opt_params.lambda_limit);
    gpr_dimer_options.opt_params.lambda_init =
        ini.GetValueF("GPR Dimer", "opt_lambda_init",
                      gpr_dimer_options.opt_params.lambda_init);
    // GPR Debugging Parameters
    gpr_dimer_options.debug_params.report_level =
        (int)ini.GetValueL("GPR Dimer", "report_level",
                           gpr_dimer_options.debug_params.report_level);
    gpr_dimer_options.debug_params.debug_level = (int)ini.GetValueL(
        "GPR Dimer", "debug_level", gpr_dimer_options.debug_params.debug_level);
    gpr_dimer_options.debug_params.out_dir =
        ini.GetValue("GPR Dimer", "debug_output_directory",
                     gpr_dimer_options.debug_params.out_dir);
    gpr_dimer_options.debug_params.pos_file =
        ini.GetValue("GPR Dimer", "debug_position_basename",
                     gpr_dimer_options.debug_params.pos_file);
    gpr_dimer_options.debug_params.energy_file =
        ini.GetValue("GPR Dimer", "debug_energy_basename",
                     gpr_dimer_options.debug_params.energy_file);
    gpr_dimer_options.debug_params.grad_file =
        ini.GetValue("GPR Dimer", "debug_gradient_basename",
                     gpr_dimer_options.debug_params.grad_file);
    gpr_dimer_options.debug_params.offset_mid_point =
        ini.GetValueF("GPR Dimer", "debug_midpoint_offset",
                      gpr_dimer_options.debug_params.offset_mid_point);
    gpr_dimer_options.debug_params.dy = ini.GetValueF(
        "GPR Dimer", "debug_y_step", gpr_dimer_options.debug_params.dy);
    gpr_dimer_options.debug_params.dz = ini.GetValueF(
        "GPR Dimer", "debug_z_step", gpr_dimer_options.debug_params.dz);
    // GPR Prune
    gpr_dimer_options.prune_params.use_prune = ini.GetValueB(
        "GPR Dimer", "use_prune", gpr_dimer_options.prune_params.use_prune);
    gpr_dimer_options.prune_params.begin = (int)ini.GetValueL(
        "GPR Dimer", "start_prune_at", gpr_dimer_options.prune_params.begin);
    gpr_dimer_options.prune_params.n_vals = (int)ini.GetValueL(
        "GPR Dimer", "nprune_vals", gpr_dimer_options.prune_params.n_vals);
    gpr_dimer_options.prune_params.threshold =
        ini.GetValueF("GPR Dimer", "prune_threshold",
                      gpr_dimer_options.prune_params.threshold);

    // [Prefactor] //

    prefactor_options.default_value = ini.GetValueF(
        "Prefactor", "default_value", prefactor_options.default_value);
    prefactor_options.max_value =
        ini.GetValueF("Prefactor", "max_value", prefactor_options.max_value);
    prefactor_options.min_value =
        ini.GetValueF("Prefactor", "min_value", prefactor_options.min_value);
    prefactor_options.within_radius = ini.GetValueF(
        "Prefactor", "within_radius", prefactor_options.within_radius);
    prefactor_options.min_displacement = ini.GetValueF(
        "Prefactor", "min_displacement", prefactor_options.min_displacement);
    prefactor_options.rate = toLowerCase(
        ini.GetValue("Prefactor", "rate_estimation", prefactor_options.rate));
    prefactor_options.configuration = toLowerCase(ini.GetValue(
        "Prefactor", "configuration", prefactor_options.configuration));
    prefactor_options.all_free_atoms = ini.GetValueB(
        "Prefactor", "all_free_atoms", prefactor_options.all_free_atoms);
    prefactor_options.filter_scheme = toLowerCase(ini.GetValue(
        "Prefactor", "filter_scheme", prefactor_options.filter_scheme));
    prefactor_options.filter_fraction = ini.GetValueF(
        "Prefactor", "filter_fraction", prefactor_options.filter_fraction);

    // [Hessian] //

    hessian_options.atom_list = toLowerCase(
        ini.GetValue("Hessian", "atom_list", hessian_options.atom_list));
    hessian_options.zero_freq_value = ini.GetValueF(
        "Hessian", "zero_freq_value", hessian_options.zero_freq_value);

    // [Nudged Elastic Band] //
    const std::string neb_section = "Nudged Elastic Band";

    // Core path parameters
    neb_options.image_count =
        ini.GetValueL(neb_section, "images", neb_options.image_count);
    neb_options.max_iterations = ini.GetValueL(
        neb_section, "max_iterations", optimizer_options.max_iterations);
    neb_options.force_tolerance = ini.GetValueF(
        neb_section, "converged_force", optimizer_options.converged_force);
    auto neb_optMethod = magic_enum::enum_cast<OptType>(
                             ini.GetValue(neb_section, "opt_method", "none"),
                             magic_enum::case_insensitive)
                             .value_or(OptType::Unknown);
    if (neb_optMethod != OptType::None) {
      neb_options.opt_method = neb_optMethod;
    }
    neb_options.mmf_peaks.enabled = ini.GetValueB(
        neb_section, "setup_mmf_peaks", neb_options.mmf_peaks.enabled);
    neb_options.mmf_peaks.tolerance = ini.GetValueF(
        neb_section, "mmf_peak_tolerance", neb_options.mmf_peaks.tolerance);

    // Inter-image spring configuration
    neb_options.spring.constant =
        ini.GetValueF(neb_section, "spring", neb_options.spring.constant);
    neb_options.spring.use_elastic_band = ini.GetValueB(
        neb_section, "elastic_band", neb_options.spring.use_elastic_band);
    neb_options.spring.doubly_nudged = ini.GetValueB(
        neb_section, "doubly_nudged", neb_options.spring.doubly_nudged);
    neb_options.spring.use_switching =
        ini.GetValueB(neb_section, "doubly_nudged_switching",
                      neb_options.spring.use_switching);

    // Energy weighting for resolution control
    neb_options.spring.weighting.enabled = ini.GetValueB(
        neb_section, "energy_weighted", neb_options.spring.weighting.enabled);
    neb_options.spring.weighting.trigger = ini.GetValueF(
        neb_section, "ew_trigger", neb_options.spring.weighting.trigger);
    neb_options.spring.weighting.k_min = ini.GetValueF(
        neb_section, "ew_ksp_min", neb_options.spring.weighting.k_min);
    neb_options.spring.weighting.k_max = ini.GetValueF(
        neb_section, "ew_ksp_max", neb_options.spring.weighting.k_max);

    // Onsager-Machlup settings
    neb_options.spring.om.enabled = ini.GetValueB(
        neb_section, "onsager_machlup", neb_options.spring.om.enabled);
    neb_options.spring.om.optimize_k = ini.GetValueB(
        neb_section, "om_optimize_k", neb_options.spring.om.optimize_k);
    neb_options.spring.om.k_scale =
        ini.GetValueF(neb_section, "om_k_scale", neb_options.spring.om.k_scale);
    neb_options.spring.om.k_min =
        ini.GetValueF(neb_section, "om_k_min", neb_options.spring.om.k_min);
    neb_options.spring.om.k_max =
        ini.GetValueF(neb_section, "om_k_max", neb_options.spring.om.k_max);

    // Climbing Image (CI-NEB) for saddle point identification
    neb_options.climbing_image.enabled =
        ini.GetValueB(neb_section, "climbing_image_method",
                      neb_options.climbing_image.enabled);
    neb_options.climbing_image.converged_only =
        ini.GetValueB(neb_section, "climbing_image_converged_only",
                      neb_options.climbing_image.converged_only);
    neb_options.climbing_image.use_old_tangent = ini.GetValueB(
        neb_section, "old_tangent", neb_options.climbing_image.use_old_tangent);
    neb_options.climbing_image.trigger_force = ini.GetValueF(
        neb_section, "ci_after", neb_options.climbing_image.trigger_force);
    neb_options.climbing_image.trigger_factor = ini.GetValueF(
        neb_section, "ci_after_rel", neb_options.climbing_image.trigger_factor);

    // Hybrid Dimer / Min-Mode Following refinement
    auto &roneb = neb_options.climbing_image.roneb;
    roneb.use_mmf = ini.GetValueB(neb_section, "ci_mmf", roneb.use_mmf);
    roneb.trigger_force =
        ini.GetValueF(neb_section, "ci_mmf_after", roneb.trigger_force);
    roneb.trigger_factor =
        ini.GetValueF(neb_section, "ci_mmf_after_rel", roneb.trigger_factor);
    roneb.max_steps =
        ini.GetValueL(neb_section, "ci_mmf_nsteps", roneb.max_steps);
    roneb.ci_stability_count = ini.GetValueL(
        neb_section, "ci_mmf_ci_stability_count", roneb.ci_stability_count);
    roneb.angle_tol =
        ini.GetValueF(neb_section, "ci_mmf_angle", roneb.angle_tol);
    roneb.penalty.strength = ini.GetValueF(
        neb_section, "ci_mmf_penalty_strength", roneb.penalty.strength);
    roneb.penalty.base =
        ini.GetValueF(neb_section, "ci_mmf_penalty_base", roneb.penalty.base);

    // Initial path construction (e.g., IDPP)
    auto &init = neb_options.initialization;
    init.method =
        magic_enum::enum_cast<NEBInit>(ini.GetValue(neb_section, "initializer"),
                                       magic_enum::case_insensitive)
            .value_or(NEBInit::LINEAR);
    init.input_path =
        ini.GetValue(neb_section, "initial_path_in", init.input_path);
    init.max_iterations =
        ini.GetValueL(neb_section, "init_max_iterations", init.max_iterations);
    init.nsteps = ini.GetValueL(neb_section, "init_nsteps", init.nsteps);
    init.max_move = ini.GetValueF(neb_section, "init_max_move", init.max_move);
    init.force_tolerance = ini.GetValueF(neb_section, "init_force_threshold",
                                         init.force_tolerance);
    init.sidpp_alpha =
        ini.GetValueF(neb_section, "sidpp_growth_alpha", init.sidpp_alpha);
    auto neb_ipath_optMethod =
        magic_enum::enum_cast<OptType>(
            ini.GetValue(neb_section, "ipath_opt_method", "none"),
            magic_enum::case_insensitive)
            .value_or(OptType::Unknown);
    if (neb_ipath_optMethod != OptType::None) {
      neb_options.initialization.opt_method = neb_ipath_optMethod;
    }
    init.oversampling =
        ini.GetValueB(neb_section, "oversampling", init.oversampling);
    init.oversampling_factor = ini.GetValueL(neb_section, "oversampling_factor",
                                             init.oversampling_factor);

    // Endpoint handling
    neb_options.endpoints.minimize = ini.GetValueB(
        neb_section, "minimize_endpoints", neb_options.endpoints.minimize);
    neb_options.endpoints.use_path_file =
        ini.GetValueB(neb_section, "minimize_endpoints_for_ipath",
                      neb_options.endpoints.use_path_file);

    // [Dynamics] //

    dynamics_options.time_step_input = ini.GetValueF(
        "Dynamics", "time_step", dynamics_options.time_step_input);
    dynamics_options.time_step =
        dynamics_options.time_step_input / constants.timeUnit;
    dynamics_options.time_input =
        ini.GetValueF("Dynamics", "time", dynamics_options.time_input);
    dynamics_options.time = dynamics_options.time_input / constants.timeUnit;
    dynamics_options.steps =
        long(floor(dynamics_options.time / dynamics_options.time_step + 0.5));
    thermostat_options.kind =
        toLowerCase(ini.GetValue("Dynamics", "thermostat", "andersen"));
    thermostat_options.andersen_alpha = ini.GetValueF(
        "Dynamics", "andersen_alpha", thermostat_options.andersen_alpha);
    thermostat_options.andersen_tcol_input =
        ini.GetValueF("Dynamics", "andersen_collision_period",
                      thermostat_options.andersen_tcol_input);
    thermostat_options.andersen_tcol =
        thermostat_options.andersen_tcol_input / constants.timeUnit;
    thermostat_options.nose_mass =
        ini.GetValueF("Dynamics", "nose_mass", thermostat_options.nose_mass);
    thermostat_options.langevin_friction_input =
        ini.GetValueF("Dynamics", "langevin_friction",
                      thermostat_options.langevin_friction_input);
    thermostat_options.langevin_friction =
        thermostat_options.langevin_friction_input * constants.timeUnit;

    // [Parallel Replica]

    parallel_replica_options.auto_stop =
        ini.GetValueB("Parallel Replica", "stop_after_transition",
                      parallel_replica_options.auto_stop);
    parallel_replica_options.refine_transition =
        ini.GetValueB("Parallel Replica", "refine_transition",
                      parallel_replica_options.refine_transition);
    parallel_replica_options.dephase_loop_stop =
        ini.GetValueB("Parallel Replica", "dephase_loop_stop",
                      parallel_replica_options.dephase_loop_stop);
    parallel_replica_options.dephase_time_input =
        ini.GetValueF("Parallel Replica", "dephase_time",
                      parallel_replica_options.dephase_time_input);
    parallel_replica_options.dephase_time =
        parallel_replica_options.dephase_time_input / constants.timeUnit;
    parallel_replica_options.dephase_loop_max =
        ini.GetValueL("Parallel Replica", "dephase_loop_max",
                      parallel_replica_options.dephase_loop_max);
    parallel_replica_options.state_check_interval_input =
        ini.GetValueF("Parallel Replica", "state_check_interval",
                      parallel_replica_options.state_check_interval_input);
    parallel_replica_options.state_check_interval =
        parallel_replica_options.state_check_interval_input /
        constants.timeUnit;
    parallel_replica_options.record_interval_input = ini.GetValueF(
        "Parallel Replica", "state_save_interval",
        0.1 * parallel_replica_options.state_check_interval_input);
    parallel_replica_options.record_interval =
        parallel_replica_options.record_interval_input / constants.timeUnit;
    parallel_replica_options.corr_time_input =
        ini.GetValueF("Parallel Replica", "post_transition_time",
                      parallel_replica_options.corr_time_input);
    parallel_replica_options.corr_time =
        parallel_replica_options.corr_time_input / constants.timeUnit;

    //[Temperature Accelerated Dynamics] //

    tad_options.low_temperature =
        ini.GetValueF("TAD", "low_temperature", tad_options.low_temperature);
    tad_options.min_prefactor =
        ini.GetValueF("TAD", "min_prefactor", tad_options.min_prefactor);
    tad_options.confidence =
        ini.GetValueF("TAD", "confidence", tad_options.confidence);

    // [Replica Exchange] //

    replica_exchange_options.temperature_distribution = toLowerCase(
        ini.GetValue("Replica Exchange", "temperature_distribution",
                     replica_exchange_options.temperature_distribution));
    replica_exchange_options.replicas = ini.GetValueL(
        "Replica Exchange", "replicas", replica_exchange_options.replicas);
    replica_exchange_options.exchange_trials =
        ini.GetValueL("Replica Exchange", "exchange_trials",
                      replica_exchange_options.exchange_trials);
    replica_exchange_options.sampling_time_input =
        ini.GetValueF("Replica Exchange", "sampling_time",
                      replica_exchange_options.sampling_time_input);
    replica_exchange_options.sampling_time =
        replica_exchange_options.sampling_time_input / constants.timeUnit;
    replica_exchange_options.temperature_low = ini.GetValueF(
        "Replica Exchange", "temperature_low", main_options.temperature);
    replica_exchange_options.temperature_high =
        ini.GetValueF("Replica Exchange", "temperature_high",
                      replica_exchange_options.temperature_high);
    replica_exchange_options.exchange_period_input =
        ini.GetValueF("Replica Exchange", "exchange_period",
                      replica_exchange_options.exchange_period_input);
    replica_exchange_options.exchange_period =
        replica_exchange_options.exchange_period_input / constants.timeUnit;

    // [Hyperdynamics] //

    hyperdynamics_options.rmd_time_input = ini.GetValueF(
        "Hyperdynamics", "bb_rmd_time", hyperdynamics_options.rmd_time_input);
    hyperdynamics_options.rmd_time =
        hyperdynamics_options.rmd_time_input / constants.timeUnit;
    hyperdynamics_options.boost_atom_list =
        toLowerCase(ini.GetValue("Hyperdynamics", "bb_boost_atomlist",
                                 hyperdynamics_options.boost_atom_list));
    hyperdynamics_options.dvmax =
        ini.GetValueF("Hyperdynamics", "bb_dvmax", hyperdynamics_options.dvmax);
    hyperdynamics_options.qrr = ini.GetValueF(
        "Hyperdynamics", "bb_stretch_threshold", hyperdynamics_options.qrr);
    hyperdynamics_options.prr = ini.GetValueF(
        "Hyperdynamics", "bb_ds_curvature", hyperdynamics_options.prr);
    hyperdynamics_options.qcut =
        ini.GetValueF("Hyperdynamics", "bb_rcut", hyperdynamics_options.qcut);
    hyperdynamics_options.bias_potential =
        toLowerCase(ini.GetValue("Hyperdynamics", "bias_potential",
                                 hyperdynamics_options.bias_potential));

    // [Saddle Search] //

    saddle_search_options.method = toLowerCase(
        ini.GetValue("Saddle Search", "method", saddle_search_options.method));
    saddle_search_options.minmode_method =
        toLowerCase(ini.GetValue("Saddle Search", "min_mode_method",
                                 saddle_search_options.minmode_method));
    saddle_search_options.displace_magnitude =
        ini.GetValueF("Saddle Search", "displace_magnitude",
                      saddle_search_options.displace_magnitude);
    saddle_search_options.displace_radius =
        ini.GetValueF("Saddle Search", "displace_radius",
                      saddle_search_options.displace_radius);
    saddle_search_options.max_energy = ini.GetValueF(
        "Saddle Search", "max_energy", saddle_search_options.max_energy);
    saddle_search_options.max_iterations = ini.GetValueL(
        "Saddle Search", "max_iterations", optimizer_options.max_iterations);
    saddle_search_options.nonnegative_displacement_abort =
        ini.GetValueB("Saddle Search", "nonnegative_displacement_abort",
                      saddle_search_options.nonnegative_displacement_abort);
    saddle_search_options.max_single_displace =
        ini.GetValueF("Saddle Search", "max_single_displace",
                      saddle_search_options.max_single_displace);
    // must be loaded after optimizer_options.converged_force
    saddle_search_options.converged_force = ini.GetValueF(
        "Saddle Search", "converged_force", optimizer_options.converged_force);
    saddle_search_options.perp_force_ratio =
        ini.GetValueF("Saddle Search", "perp_force_ratio",
                      saddle_search_options.perp_force_ratio);
    saddle_search_options.displace_type = toLowerCase(ini.GetValue(
        "Saddle Search", "client_displace_type", EpiCenters::DISP_LOAD));
    saddle_search_options.nonlocal_count_abort =
        ini.GetValueL("Saddle Search", "nonlocal_count_abort",
                      saddle_search_options.nonlocal_count_abort);
    saddle_search_options.nonlocal_distance_abort =
        ini.GetValueF("Saddle Search", "nonlocal_distance_abort",
                      saddle_search_options.nonlocal_distance_abort);
    if (saddle_search_options.displace_type !=
            EpiCenters::DISP_NOT_FCC_OR_HCP &&
        saddle_search_options.displace_type !=
            EpiCenters::DISP_MIN_COORDINATED &&
        saddle_search_options.displace_type != EpiCenters::DISP_LAST_ATOM &&
        saddle_search_options.displace_type != EpiCenters::DISP_RANDOM &&
        saddle_search_options.displace_type != EpiCenters::DISP_LISTED_ATOMS) {
      saddle_search_options.displace_type = EpiCenters::DISP_LOAD;
    }
    // Parse comma-separated atom list for listed_atoms displacement type
    {
      string atomListStr =
          ini.GetValue("Saddle Search", "displace_atom_list", "");
      if (!atomListStr.empty()) {
        std::stringstream ss(atomListStr);
        string token;
        while (std::getline(ss, token, ',')) {
          // Trim whitespace
          size_t start = token.find_first_not_of(" \t");
          size_t end = token.find_last_not_of(" \t");
          if (start != string::npos) {
            saddle_search_options.displace_atom_list.push_back(
                std::stol(token.substr(start, end - start + 1)));
          }
        }
      }
    }
    saddle_search_options.confine_positive.enabled =
        ini.GetValueB("Saddle Search", "confine_positive",
                      saddle_search_options.confine_positive.enabled);
    if (saddle_search_options.confine_positive.enabled) {
      saddle_search_options.confine_positive.bowl_breakout =
          ini.GetValueB("Saddle Search", "bowl_breakout",
                        saddle_search_options.confine_positive.bowl_breakout);
      saddle_search_options.confine_positive.bowl_active =
          ini.GetValueL("Saddle Search", "bowl_active_atoms",
                        saddle_search_options.confine_positive.bowl_active);
      saddle_search_options.confine_positive.min_force =
          ini.GetValueF("Saddle Search", "confine_positive_min_move",
                        saddle_search_options.confine_positive.min_force);
      saddle_search_options.confine_positive.scale_ratio =
          ini.GetValueF("Saddle Search", "confine_positive_scale_ratio",
                        saddle_search_options.confine_positive.scale_ratio);
      saddle_search_options.confine_positive.boost =
          ini.GetValueF("Saddle Search", "confine_positive_boost",
                        saddle_search_options.confine_positive.boost);
      saddle_search_options.confine_positive.min_active =
          ini.GetValueL("Saddle Search", "confine_positive_min_active",
                        saddle_search_options.confine_positive.min_active);
    }
    saddle_search_options.dynamics.temperature = main_options.temperature;
    saddle_search_options.dynamics.temperature =
        ini.GetValueF("Saddle Search", "dynamics_temperature",
                      saddle_search_options.dynamics.temperature);
    saddle_search_options.dynamics.state_check_interval_input = ini.GetValueF(
        "Saddle Search", "dynamics_state_check_interval",
        saddle_search_options.dynamics.state_check_interval_input);
    saddle_search_options.dynamics.state_check_interval =
        saddle_search_options.dynamics.state_check_interval_input /
        constants.timeUnit;
    saddle_search_options.dynamics.record_interval_input =
        ini.GetValueF("Saddle Search", "dynamics_record_interval",
                      saddle_search_options.dynamics.record_interval_input);
    saddle_search_options.dynamics.record_interval =
        saddle_search_options.dynamics.record_interval_input /
        constants.timeUnit;
    saddle_search_options.dynamics.linear_interpolation =
        ini.GetValueB("Saddle Search", "dynamics_linear_interpolation",
                      saddle_search_options.dynamics.linear_interpolation);
    saddle_search_options.remove_rotation =
        ini.GetValueB("Saddle Search", "remove_rotation",
                      saddle_search_options.remove_rotation);
    saddle_search_options.dynamics.max_init_curvature =
        ini.GetValueF("Saddle Search", "dynamics_max_init_curvature",
                      saddle_search_options.dynamics.max_init_curvature);
    saddle_search_options.zero_mode_abort_curvature =
        ini.GetValueF("Saddle Search", "zero_mode_abort_curvature",
                      saddle_search_options.zero_mode_abort_curvature);

    // [Basin Hopping] //

    basin_hopping_options.displacement = ini.GetValueF(
        "Basin Hopping", "displacement", basin_hopping_options.displacement);
    basin_hopping_options.push_apart_distance =
        ini.GetValueF("Basin Hopping", "push_apart_distance",
                      basin_hopping_options.push_apart_distance);
    basin_hopping_options.initial_random_structure_probability = ini.GetValueF(
        "Basin Hopping", "initial_random_structure_probability",
        basin_hopping_options.initial_random_structure_probability);
    basin_hopping_options.steps =
        ini.GetValueL("Basin Hopping", "steps", basin_hopping_options.steps);
    basin_hopping_options.quenching_steps =
        ini.GetValueL("Basin Hopping", "quenching_steps",
                      basin_hopping_options.quenching_steps);
    basin_hopping_options.single_atom_displace =
        ini.GetValueB("Basin Hopping", "single_atom_displace",
                      basin_hopping_options.single_atom_displace);
    basin_hopping_options.significant_structure =
        ini.GetValueB("Basin Hopping", "significant_structure",
                      basin_hopping_options.significant_structure);
    basin_hopping_options.displacement_algorithm =
        toLowerCase(ini.GetValue("Basin Hopping", "displacement_algorithm",
                                 basin_hopping_options.displacement_algorithm));
    basin_hopping_options.displacement_distribution = toLowerCase(
        ini.GetValue("Basin Hopping", "displacement_distribution",
                     basin_hopping_options.displacement_distribution));
    basin_hopping_options.swap_probability =
        ini.GetValueF("Basin Hopping", "swap_probability",
                      basin_hopping_options.swap_probability);
    basin_hopping_options.jump_max = ini.GetValueL(
        "Basin Hopping", "jump_max", basin_hopping_options.jump_max);
    basin_hopping_options.jump_steps = ini.GetValueL(
        "Basin Hopping", "jump_steps", basin_hopping_options.jump_steps);
    basin_hopping_options.adjust_displacement =
        ini.GetValueB("Basin Hopping", "adjust_displacement",
                      basin_hopping_options.adjust_displacement);
    basin_hopping_options.adjust_period = ini.GetValueL(
        "Basin Hopping", "adjust_period", basin_hopping_options.adjust_period);
    basin_hopping_options.adjust_fraction =
        ini.GetValueF("Basin Hopping", "adjust_fraction",
                      basin_hopping_options.adjust_fraction);
    basin_hopping_options.target_ratio = ini.GetValueF(
        "Basin Hopping", "target_ratio", basin_hopping_options.target_ratio);
    basin_hopping_options.write_unique = ini.GetValueB(
        "Basin Hopping", "write_unique", basin_hopping_options.write_unique);
    basin_hopping_options.stop_energy = ini.GetValueF(
        "Basin Hopping", "stop_energy", basin_hopping_options.stop_energy);

    // [Global Optimization] //

    global_optimization_options.move_method =
        toLowerCase(ini.GetValue("Global Optimization", "move_method",
                                 global_optimization_options.move_method));
    global_optimization_options.decision_method =
        toLowerCase(ini.GetValue("Global Optimization", "decision_method",
                                 global_optimization_options.decision_method));
    global_optimization_options.steps = ini.GetValueL(
        "Global Optimization", "steps", global_optimization_options.steps);
    global_optimization_options.beta = ini.GetValueF(
        "Global Optimization", "beta", global_optimization_options.beta);
    global_optimization_options.alpha = ini.GetValueF(
        "Global Optimization", "alpha", global_optimization_options.alpha);
    global_optimization_options.mdmin = ini.GetValueL(
        "Global Optimization", "mdmin", global_optimization_options.mdmin);
    global_optimization_options.target_energy =
        ini.GetValueF("Global Optimization", "target_energy",
                      global_optimization_options.target_energy);

    // [BGSD] //

    bgsd_options.alpha = ini.GetValueF("BGSD", "alpha", bgsd_options.alpha);
    bgsd_options.beta = ini.GetValueF("BGSD", "beta", bgsd_options.beta);
    bgsd_options.gradient_finite_difference =
        ini.GetValueF("BGSD", "gradientfinitedifference",
                      bgsd_options.gradient_finite_difference);
    bgsd_options.grad2energy_convergence = ini.GetValueF(
        "BGSD", "grad2energyconvergence", bgsd_options.grad2energy_convergence);
    bgsd_options.grad2force_convergence = ini.GetValueF(
        "BGSD", "grad2forceconvergence", bgsd_options.grad2force_convergence);

    // [Monte Carlo] //

    monte_carlo_options.step_size = ini.GetValueF(
        "Monte Carlo", "step_size", monte_carlo_options.step_size);
    monte_carlo_options.steps =
        ini.GetValueI("Monte Carlo", "steps", monte_carlo_options.steps);

    // Sanity Checks
    if (parallel_replica_options.state_check_interval > dynamics_options.time &&
        magic_enum::enum_name<JobType>(main_options.job) ==
            "parallel_replica") {
      SPDLOG_ERROR("[Parallel Replica] state_check_interval must be <= time");
      error = 1;
    }

    // Check if an initial path exists without a specific non-linear initializer
    if (!neb_options.initialization.input_path.empty() &&
        neb_options.initialization.method == NEBInit::LINEAR) {
      SPDLOG_WARN("[Nudged Elastic Band] 'initial_path_in' is provided, but "
                  "'initializer' defaults to linear. "
                  "Ensure this is intentional, as the loaded path will not be "
                  "used without initializer set to file.");
    }

    if (saddle_search_options.dynamics.record_interval_input >
        saddle_search_options.dynamics.state_check_interval_input) {
      SPDLOG_ERROR("[Saddle Search] dynamics_record_interval must be <= "
                   "dynamics_state_check_interval");
      error = 1;
    }

    if (potential_options.potential == PotType::AMS ||
        potential_options.potential == PotType::AMS_IO) {
      if (ams_options.forcefield.empty() && ams_options.model.empty() &&
          ams_options.xc.empty()) {
        SPDLOG_ERROR("[AMS] Must provide atleast forcefield or model or xc");
        error = 1;
      }

      if (!ams_options.forcefield.empty() && !ams_options.model.empty() &&
          !ams_options.xc.empty()) {
        SPDLOG_ERROR("[AMS] Must provide either forcefield or model");
        error = 1;
      }
    }

  } else {
    SPDLOG_ERROR("Couldn't parse the ini file");
    error = 1;
  }

  return error;
}

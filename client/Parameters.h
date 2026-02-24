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
#pragma once

#include "BaseStructures.h"
#include <cstring>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

#ifdef EONMPI
#include "mpi.h"
#endif

/** Contains all runtime parameters and results. No functionality just
 * bookkeeping.*/
class Parameters {

public:
  Parameters();
  ~Parameters() = default;
  // TODO: Is this complete?
  Parameters(const Parameters &) = default;
  int load(string filename);
  int load(FILE *file);

  /** string constants: declared here, defined in Parameters.cpp. **/

  // potentials //
  // jobs //

  // Physical Constants
  struct constants_t {
    double kB;
    double timeUnit;
  } constants;

  /** input parameters **/

  // [Main] //
  struct main_options_t {
    JobType job;
    long randomSeed;
    double temperature;
    bool quiet;
    bool writeLog;
    bool checkpoint;
    string iniFilename;
    string conFilename;
    double finiteDifference;
    long maxForceCalls;
    bool removeNetForce;
  } main_options;

  // [Potential] //
  struct potential_options_t {
    PotType potential;
    double MPIPollPeriod;
    bool LAMMPSLogging;
    int LAMMPSThreads;
    bool EMTRasmussen;
    bool LogPotential;
    string extPotPath;
    int MPIPotentialRank;
#ifdef EONMPI
    MPI_Comm MPIClientComm;
#endif
  } potential_options;

  // [AMS] and [AMS_IO] //
  struct ams_options_t {
    string engine;     // MOPAC, ADF, BAND, REAXFF, FORCEFIELD
    string forcefield; // OPt.ff etc. (REAXFF)
    string model;      // Model hamiltonian (MOPAC)
    string resources;  // DFTB
    string xc;         // Exchange (BAND, ADF)
    string basis;      // Basis (BAND, ADF)
    struct env_t {
      string amshome;
      string scm_tmpdir;
      string scmlicense;
      string scm_pythondir;
      string amsbin;
      string amsresources;
    } env;
  } ams_options;

  // [XTBPot] //
  struct xtb_options_t {
    string paramset;
    double elec_temperature;
    size_t maxiter;
    double acc;
    double charge;
    int uhf;
  } xtb_options;

  // [ZBLPot] //
  struct zbl_options_t {
    double cut_inner;
    double cut_global;
  } zbl_options;

  // [SocketNWChemPot] //
  struct socket_nwchem_options_t {
    std::string host;
    int port;
    int mem_in_gb;
    std::string nwchem_settings;
    std::string unix_socket_path;
    bool unix_socket_mode;
    bool make_template_input;
  } socket_nwchem_options;

  // [Structure Comparison] //
  struct structure_comparison_options_t {
    double distance_difference;
    double neighbor_cutoff;
    bool check_rotation;
    bool indistinguishable_atoms;
    double energy_difference;
    bool remove_translation;
  } structure_comparison_options;

  // [Process Search] //
  struct process_search_options_t {
    bool minimize_first;
    double minimization_offset;
  } process_search_options;

  // [Saddle Search] //
  struct saddle_search_options_t {
    long max_jump_attempts;
    long max_iterations;
    string method;
    string minmode_method;
    string displace_type;
    std::vector<long> displace_atom_list;
    double max_energy;
    double displace_magnitude;
    double max_single_displace;
    double displace_radius;
    double converged_force;
    double perp_force_ratio;
    bool nonnegative_displacement_abort;
    long nonlocal_count_abort;
    double nonlocal_distance_abort;
    bool remove_rotation;
    double zero_mode_abort_curvature;
    struct dynamics_t {
      double temperature;
      double state_check_interval_input;
      double state_check_interval;
      double record_interval_input;
      double record_interval;
      bool linear_interpolation;
      double max_init_curvature;
    } dynamics;
    struct confine_positive_t {
      bool enabled;
      bool bowl_breakout;
      long bowl_active;
      double min_force;
      double scale_ratio;
      double boost;
      long min_active;
    } confine_positive;
  } saddle_search_options;

  // [Optimizer] //
  struct optimizer_options_t {
    OptType method;
    string convergence_metric; // norm, max_atom, max_component
    string convergence_metric_label;
    size_t max_iterations;
    double max_move;
    double converged_force;
    double time_step_input;
    double time_step;
    double max_time_step_input;
    double max_time_step;
    struct refine_t {
      OptType method;
      double threshold;
    } refine;
    struct lbfgs_t {
      long memory;
      double inverse_curvature;
      double max_inverse_curvature;
      bool auto_scale;
      bool angle_reset;
      bool distance_reset;
    } lbfgs;
    struct cg_t {
      bool no_overshooting;
      bool knock_out_max_move;
      bool line_search;
      double line_converged;
      long line_search_max_iter;
      long max_iter_before_reset;
    } cg;
    struct quickmin_t {
      bool steepest_descent;
    } quickmin;
    struct sd_t {
      double alpha;
      bool two_point;
    } sd;
  } optimizer_options;

  // [Dimer] //
  struct dimer_options_t {
    double rotation_angle;
    bool improved;
    double converged_angle;
    long max_iterations;
    string opt_method;
    long rotations_max;
    long rotations_min;
    double torque_max;
    double torque_min;
    bool remove_rotation;
  } dimer_options;

  // [GPR Dimer] //
  struct gpr_dimer_options_t {
    double rotation_angle;
    double converged_angle;
    double relax_conv_angle;
    long init_rotations_max;
    long relax_rotations_max;
    long divisor_t_dimer_gp;
    long max_outer_iterations;
    long max_inner_iterations;
    double midpoint_max_disp;
    string rot_opt_method;
    string trans_opt_method;
    double active_radius;
    double dimer_sep;
    double conv_step;
    double max_step;
    double ratio_at_limit;
    bool init_rot_gp;
    bool init_trans_gp;
    bool many_iterations;
    struct gpr_params_t {
      string hyper_opt_method;
      double sigma2;
      double jitter_sigma2;
      double noise_sigma2;
      double prior_mu;
      double prior_sigma2;
      long prior_nu;
    } gpr_params;
    struct opt_params_t {
      bool check_derivatives;
      int max_iterations;
      double tol_func;
      double tol_sol;
      long lambda_limit;
      long lambda_init;
    } opt_params;
    struct prune_params_t {
      bool use_prune;
      int begin;
      int n_vals;
      double threshold;
    } prune_params;
    struct debug_params_t {
      int report_level;
      int debug_level;
      string out_dir;
      string pos_file;
      string energy_file;
      string grad_file;
      string out_ext;
      double offset_mid_point;
      double dy;
      double dz;
    } debug_params;
  } gpr_dimer_options;

  // [GP Surrogate] //
  struct gp_surrogate_options_t {
    bool enabled;
    JobType sub_job;
    double uncertainty;
    bool linear_path_always;
    PotType potential;
  } gp_surrogate_options;

  // [CatLearn] //
  struct catlearn_options_t {
    std::string path;
    std::string model;
    std::string prior;
    bool use_deriv;
    bool use_fingerprint;
    bool parallel;
  } catlearn_options;

  // [ASE_ORCA] //
  struct ase_orca_options_t {
    std::string path;
    std::string nproc;
    std::string simpleinput;
  } ase_orca_options;

  // [ASE_NWCHEM] //
  struct ase_nwchem_options_t {
    std::string path;
    std::string nproc;
    std::string multiplicity;
    double scf_thresh;
    long scf_maxiter;
  } ase_nwchem_options;

  // [Metatomic] //
  struct metatomic_options_t {
    std::string model_path;  // Path to the TorchScript model file.
    std::string device;      // "cpu", "cuda", "mps", or empty to auto-detect.
    std::string length_unit; // The unit of length used in the simulation (e.g.,
                             // "angstrom").
    std::string extensions_directory; // Path for TorchScript extensions.
    bool check_consistency;           // To enable model's internal checks.
    double uncertainty_threshold;     // Threshold for uncertainty reporting.
                                      // also used to populate the variance
                                      // -1 to disable, 100meV/atom default
    struct variants_t {
      std::string base;               // global key for energy and forces
      std::string energy;             // override for energy variant
      std::string energy_uncertainty; // override for variant on the
                                      // energy uncertainty
    } variant;
  } metatomic_options;

  // [Lanczos] //
  struct lanczos_options_t {
    double tolerance;
    long max_iterations;
    bool quit_early;
  } lanczos_options;

  // [Prefactor] //
  struct prefactor_options_t {
    double default_value;
    double max_value;
    double min_value;
    double within_radius;
    double min_displacement;
    string rate;
    string configuration;
    bool all_free_atoms;
    string filter_scheme;
    double filter_fraction;
  } prefactor_options;

  // [Hessian] //
  struct hessian_options_t {
    string atom_list;
    double zero_freq_value;
  } hessian_options;

  // [Nudged Elastic Band] //
  struct neb_options_t {
    // Core parameters for the chain-of-states simulation
    long image_count;    // Number of replicas along the reaction coordinate
    long max_iterations; // Maximum steps for the path optimization
    OptType opt_method; // Optimization algorithm (e.g., QuickMin, FIRE, L-BFGS)
    double
        force_tolerance; // Convergence criterion for the root-mean-square force
    struct mmf_peak_options_t {
      bool enabled;     // Initialize mode estimates for each peak
      double tolerance; // Cutoff for generating peaks
    } mmf_peaks;

    struct spring_options_t {
      double constant;       // The spring constant (k) connecting images
      bool use_elastic_band; // Toggle for the elastic band projection
      bool doubly_nudged;    // Inclusion of the perpendicular spring force
                             // components
      bool use_switching;    // Switching function for doubly nudged NEB

      struct energy_weighting_t {
        bool
            enabled; // Adjusts spring constants based on image potential energy
        double trigger; // Threshold to activate energy weighted springs
        double k_min;   // Minimum spring constant for low-energy regions
        double k_max; // Maximum spring constant for high-energy barrier regions
      } weighting;
      struct onsager_machlup_t {
        bool enabled; // Main toggle for OM-NEB

        // Adaptive Spring Constant Logic
        bool optimize_k; // Re-calculate k_om every step?
        double k_scale;  // Scaling factor for the heuristic (dimensionless)

        // Optional: Safety bounds to prevent numerical explosions
        double k_min;
        double k_max;
      } om;
    } spring;

    struct climbing_image_options_t {
      bool enabled; // Enables the Climbing Image (CI-NEB) modification
      bool
          converged_only; // Wait for initial MEP convergence before starting CI
      bool use_old_tangent;  // Algorithm choice for the local tangent estimate
      double trigger_force;  // Force threshold to activate the CI-NEB algorithm
      double trigger_factor; // relative factor

      struct hybrid_dimer_t {
        bool use_mmf; // Integrates Min-Mode Following (MMF) for saddle point
                      // refinement
        double trigger_force; // Force threshold to activate hybrid dimer search
        long max_steps;       // Maximum steps for the MMF refinement phase
        long ci_stability_count; // Number of stable iterations before settling
                                 // on a CI for MMF
        double angle_tol;        // Angular tolerance for the dimer rotation
        double trigger_factor;   // relative factor
        struct roneb_penalty_t {
          // penalty_factor = base + (strength * alignment);
          // std::abs(finalMode.normalized().dot(currentTangent.normalized()))
          double strength; // amount by which the threshold is reduced
          double base;     // baseline
        } penalty;
      } roneb;
    } climbing_image;

    struct path_initialization_t {
      NEBInit method; // Method for generating the initial guess (e.g., IDPP)
      std::string
          input_path;     // Path to a file containing initial image coordinates
      int max_iterations; // Maximum iterations for the path pre-optimizer
      int nsteps;         // Iterations for the path pre-optimizer
      double max_move;    // Maximum displacement per step during initialization
      double force_tolerance; // Convergence criterion for the initial path
                              // generation
      double sidpp_alpha; // Growth parameter for Steepest Intensity Decent Path
                          // Probing
      OptType opt_method; // for the IDPP
      bool oversampling;
      int oversampling_factor;
    } initialization;

    struct endpoint_options_t {
      bool minimize; // Flag to optimize the reactant and product geometries
                     // first
      bool use_path_file; // Pull endpoint geometries from the initial path file
    } endpoints;

  } neb_options;

  // [Molecular Dynamics] //
  struct dynamics_options_t {
    double time_step_input;
    double time_step;
    double time_input;
    double time;
    long steps;
  } dynamics_options;

  // [Parallel Replica] //
  struct parallel_replica_options_t {
    bool refine_transition;
    bool auto_stop;
    bool dephase_loop_stop;
    double dephase_time_input;
    double dephase_time;
    long dephase_loop_max;
    double state_check_interval_input;
    double state_check_interval;
    double record_interval_input;
    double record_interval;
    double corr_time_input;
    double corr_time;
  } parallel_replica_options;

  // [Temperature Accelerated Dynamics] //
  struct tad_options_t {
    double low_temperature;
    double min_prefactor;
    double confidence;
  } tad_options;

  // [Thermostat] //
  struct thermostat_options_t {
    string kind;
    double andersen_alpha;
    double andersen_tcol_input;
    double andersen_tcol;
    double nose_mass;
    double langevin_friction_input;
    double langevin_friction;
  } thermostat_options;

  // [Replica Exchange] //
  struct replica_exchange_options_t {
    string temperature_distribution;
    long replicas;
    long exchange_trials;
    double sampling_time_input;
    double sampling_time;
    double temperature_high;
    double temperature_low;
    double exchange_period_input;
    double exchange_period;
  } replica_exchange_options;

  // [Bond Boost / Hyperdynamics] //
  struct hyperdynamics_options_t {
    string bias_potential;
    string boost_atom_list;
    double rmd_time_input;
    double rmd_time;
    double dvmax;
    double qrr;
    double prr;
    double qcut;
  } hyperdynamics_options;

  // [Basin Hopping] //
  struct basin_hopping_options_t {
    double displacement;
    double initial_random_structure_probability;
    double push_apart_distance;
    long steps;
    long quenching_steps;
    bool significant_structure;
    bool single_atom_displace;
    string displacement_algorithm;
    string displacement_distribution;
    double swap_probability;
    long jump_max;
    long jump_steps;
    bool adjust_displacement;
    long adjust_period;
    double adjust_fraction;
    double target_ratio;
    bool write_unique;
    double stop_energy;
  } basin_hopping_options;

  // [Global Optimization] //
  struct global_optimization_options_t {
    string move_method;
    string decision_method;
    long steps;
    double beta;
    double alpha;
    long mdmin;
    double target_energy;
  } global_optimization_options;

  // [Monte Carlo] //
  struct monte_carlo_options_t {
    double step_size;
    int steps;
  } monte_carlo_options;

  // [BGSD] //
  struct bgsd_options_t {
    double alpha;
    double beta;
    double gradient_finite_difference;
    double h_force_convergence;
    double grad2energy_convergence;
    double grad2force_convergence;
  } bgsd_options;

  // [Serve] //
  struct serve_options_t {
    string host;
    uint16_t port;
    size_t replicas;
    uint16_t gateway_port; // 0 = disabled
    string endpoints;      // "pot:port,pot:host:port,..." spec string
  } serve_options;

  // [Debug] //
  struct debug_options_t {
    bool write_movies;
    long write_movies_interval;
    bool estimate_neb_eigenvalues;
    string neb_mmf;
  } debug_options;

private:
  string toLowerCase(string s);
};

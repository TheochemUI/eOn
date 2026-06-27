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
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#ifdef EONMPI
#include "mpi.h"
#endif

/** Contains all runtime parameters and results. No functionality just
 * bookkeeping.*/
namespace eonc {

class Parameters {

public:
  Parameters();
  ~Parameters() = default;
  Parameters(const Parameters &) = default;
  int load(std::string filename);
  int load(FILE *file);
  int load_json(const std::string &json_str);
  std::string to_json() const;

  // Physical Constants
  struct constants_t {
    double kB{8.6173324e-5};     // eV/K
    double timeUnit{10.1805055}; // fs
  } constants;

  // [Main] //
  struct main_options_t {
    JobType job{JobType::Process_Search};
    long randomSeed{-1};
    double temperature{300.0};
    bool quiet{false};
    bool writeLog{true};
    bool checkpoint{false};
    std::string iniFilename{"config.ini"};
    std::string conFilename{"pos.con"};
    double finiteDifference{0.01};
    long maxForceCalls{0};
    bool removeNetForce{true};
    bool parallel{true}; // parallel force evaluation via std::jthread
  } main_options;

  // [Potential] //
  struct potential_options_t {
    PotType potential{PotType::LJ};
    double MPIPollPeriod{0.25};
    bool LAMMPSLogging{false};
    int LAMMPSThreads{0};
    bool EMTRasmussen{false};
    bool LogPotential{false};
    std::string extPotPath{"./ext_pot"};
    std::string potentialsPath{
        ""}; // colon-separated dirs for Fortran potential .so files
    int MPIPotentialRank{-1};
#ifdef EONMPI
    MPI_Comm MPIClientComm;
#endif
  } potential_options;

  // [AMS] and [AMS_IO] //
  struct ams_options_t {
    std::string engine;
    std::string forcefield;
    std::string model;
    std::string resources;
    std::string xc;
    std::string basis;
    struct env_t {
      std::string amshome;
      std::string scm_tmpdir;
      std::string scmlicense;
      std::string scm_pythondir;
      std::string amsbin;
      std::string amsresources;
    } env;
  } ams_options;

  // [XTBPot] //
  struct xtb_options_t {
    std::string paramset{"GFNFF"};
    double elec_temperature{0.0};
    size_t maxiter{250};
    double acc{1.0};
    double charge{0.0};
    int uhf{0};
  } xtb_options;

  // [ZBLPot] //
  struct zbl_options_t {
    double cut_inner{2.0};
    double cut_global{2.5};
  } zbl_options;

  // [SocketNWChemPot] //
  struct socket_nwchem_options_t {
    std::string host{"127.0.0.1"};
    int port{9999};
    int mem_in_gb{2};
    std::string nwchem_settings{"nwchem_settings.nwi"};
    std::string unix_socket_path{"eon_nwchem"};
    bool unix_socket_mode{false};
    bool make_template_input{true};
  } socket_nwchem_options;

  // [Structure Comparison] //
  struct structure_comparison_options_t {
    double distance_difference{0.1};
    double neighbor_cutoff{3.3};
    bool check_rotation{false};
    bool indistinguishable_atoms{true};
    double energy_difference{0.01};
    bool remove_translation{true};
  } structure_comparison_options;

  // [Process Search] //
  struct process_search_options_t {
    bool minimize_first{true};
    double minimization_offset{0.2}; // resolved to optimizer_options.max_move
  } process_search_options;

  // [Saddle Search] //
  struct saddle_search_options_t {
    long max_jump_attempts{0};
    long max_iterations{1000};
    std::string method{"min_mode"};
    std::string minmode_method{"dimer"};
    std::string displace_type{"load"};
    std::vector<long> displace_atom_list;
    double max_energy{20.0};
    double displace_magnitude{0.1};
    double max_single_displace{10.0};
    double displace_radius{4.0};
    double converged_force{0.0};
    double perp_force_ratio{0.0};
    bool nonnegative_displacement_abort{false};
    long nonlocal_count_abort{0};
    double nonlocal_distance_abort{0.0};
    bool remove_rotation{false};
    double zero_mode_abort_curvature{0.0};
    struct dynamics_t {
      double temperature{0.0};
      double state_check_interval_input{100.0};
      double state_check_interval{0.0}; // computed: input / timeUnit
      double record_interval_input{10.0};
      double record_interval{0.0}; // computed
      bool linear_interpolation{true};
      double max_init_curvature{0.0};
    } dynamics;
    struct confine_positive_t {
      bool enabled{false};
      bool bowl_breakout{false};
      long bowl_active{20};
      double min_force{0.5};
      double scale_ratio{0.9};
      double boost{10.0};
      long min_active{30};
    } confine_positive;
  } saddle_search_options;

  // [Optimizer] //
  struct optimizer_options_t {
    OptType method{OptType::CG};
    std::string convergence_metric{"norm"};
    std::string convergence_metric_label;
    size_t max_iterations{1000};
    double max_move{0.2};
    double converged_force{0.01};
    double time_step_input{1.0};
    double time_step{0.0}; // computed: input / timeUnit
    double max_time_step_input{2.5};
    double max_time_step{0.0}; // computed: input / timeUnit
    struct refine_t {
      OptType method{OptType::None};
      double threshold{0.5};
    } refine;
    struct lbfgs_t {
      long memory{20};
      double inverse_curvature{0.01};
      double max_inverse_curvature{0.0};
      bool auto_scale{true};
      bool angle_reset{true};
      bool distance_reset{true};
    } lbfgs;
    struct cg_t {
      bool no_overshooting{false};
      bool knock_out_max_move{false};
      bool line_search{false};
      double line_converged{0.1};
      long line_search_max_iter{10};
      long max_iter_before_reset{0};
    } cg;
    struct quickmin_t {
      bool steepest_descent{false};
    } quickmin;
    struct sd_t {
      double alpha{0.1};
      bool two_point{false};
    } sd;
  } optimizer_options;

  // [Dimer] //
  struct dimer_options_t {
    double rotation_angle{0.005};
    bool improved{true};
    double converged_angle{5.0};
    long max_iterations{1000};
    std::string opt_method{"cg"};
    long rotations_max{10};
    long rotations_min{1};
    double torque_max{1.0};
    double torque_min{0.1};
    bool remove_rotation{false};
  } dimer_options;

  // [GPR Dimer] //
  struct gpr_dimer_options_t {
    double rotation_angle{0.005};
    double converged_angle{0.08};
    double relax_conv_angle{0.001};
    long init_rotations_max{6};
    long relax_rotations_max{10};
    long divisor_t_dimer_gp{10};
    long max_outer_iterations{300};
    long max_inner_iterations{1000};
    double midpoint_max_disp{0.5};
    std::string rot_opt_method{"lbfgs"};
    std::string trans_opt_method{"lbfgs"};
    double active_radius{5.0};
    double dimer_sep{0.01};
    double conv_step{0.1};
    double max_step{0.1};
    double ratio_at_limit{0.66667};
    bool init_rot_gp{false};
    bool init_trans_gp{false};
    bool many_iterations{true};
    struct gpr_params_t {
      std::string hyper_opt_method{"scg"};
      double sigma2{1e-8};
      double jitter_sigma2{0.0};
      double noise_sigma2{1e-8};
      double prior_mu{0.0};
      double prior_sigma2{1.0};
      long prior_nu{20};
    } gpr_params;
    struct opt_params_t {
      bool check_derivatives{false};
      int max_iterations{400};
      double tol_func{1e-4};
      double tol_sol{1e-4};
      long lambda_limit{static_cast<long>(1e17)};
      long lambda_init{10};
    } opt_params;
    struct prune_params_t {
      bool use_prune{false};
      int begin{8};
      int n_vals{3};
      double threshold{0.5};
    } prune_params;
    struct debug_params_t {
      int report_level{1};
      int debug_level{2};
      std::string out_dir{"output"};
      std::string pos_file{"position"};
      std::string energy_file{"energy"};
      std::string grad_file{"gradient"};
      std::string out_ext{"dat"};
      double offset_mid_point{3.0};
      double dy{0.1};
      double dz{0.1};
    } debug_params;
  } gpr_dimer_options;

  // [GP Surrogate] //
  struct gp_surrogate_options_t {
    bool enabled{false};
    JobType sub_job{JobType::Unknown};
    double uncertainty{0.05};
    bool linear_path_always{false};
    PotType potential{PotType::CatLearn};
  } gp_surrogate_options;

  // [CatLearn] //
  struct catlearn_options_t {
    std::string path;
    std::string model{"gp"};
    std::string prior{"median"};
    bool use_deriv{true};
    bool use_fingerprint{false};
    bool parallel{false};
  } catlearn_options;

  // [ASE_ORCA] //
  struct ase_orca_options_t {
    std::string path;
    std::string nproc{"1"};
    std::string simpleinput;
  } ase_orca_options;

  // [ASE_NWCHEM] //
  struct ase_nwchem_options_t {
    std::string path;
    std::string nproc{"1"};
    std::string multiplicity;
    double scf_thresh{1e-5};
    long scf_maxiter{200};
  } ase_nwchem_options;

  // [Metatomic] //
  struct metatomic_options_t {
    std::string model_path;
    std::string device{"cpu"};
    std::string length_unit{"angstrom"};
    std::string extensions_directory;
    bool check_consistency{false};
    double uncertainty_threshold{-1.0};
    struct variants_t {
      std::string base;
      std::string energy;
      std::string energy_uncertainty;
    } variant;
  } metatomic_options;

  // [Lanczos] //
  struct lanczos_options_t {
    double tolerance{0.01};
    long max_iterations{20};
    bool quit_early{true};
  } lanczos_options;

  // [Davidson] min-mode (alternative to dimer rotation / Lanczos)
  struct davidson_options_t {
    double tolerance{0.01};
    long max_iterations{20};
    /// Heuristic | (H v)_i / v_i | preconditioner (not true diag(H); off by
    /// default).
    bool diagonal_preconditioner{false};
  } davidson_options;

  // [Prefactor] //
  struct prefactor_options_t {
    double default_value{0.0};
    double max_value{1e+21};
    double min_value{1e+9};
    double within_radius{3.3};
    double min_displacement{0.25};
    std::string rate{"htst"};
    std::string configuration{"reactant"};
    bool all_free_atoms{false};
    std::string filter_scheme{"fraction"};
    double filter_fraction{0.90};
  } prefactor_options;

  // [Hessian] //
  struct hessian_options_t {
    std::string atom_list{"All"};
    double zero_freq_value{1e-6};
  } hessian_options;

  // [Nudged Elastic Band] //
  struct neb_options_t {
    long image_count{5};
    long max_iterations{1000};
    OptType opt_method{OptType::LBFGS};
    double force_tolerance{
        0.01}; // resolved to optimizer_options.converged_force
    struct mmf_peak_options_t {
      bool enabled{true};
      double tolerance{0.05};
    } mmf_peaks;

    struct spring_options_t {
      double constant{5.0};
      bool use_elastic_band{false};
      bool doubly_nudged{false};
      bool use_switching{false};

      struct energy_weighting_t {
        bool enabled{false};
        double trigger{10.0};
        double k_min{0.97};
        double k_max{9.7};
      } weighting;
      struct onsager_machlup_t {
        bool enabled{false};
        bool optimize_k{true};
        double k_scale{1.0};
        double k_min{0.1};
        double k_max{100.0};
      } om;
    } spring;

    struct climbing_image_options_t {
      bool enabled{true};
      bool converged_only{true};
      bool use_old_tangent{false};
      double trigger_force{std::numeric_limits<double>::infinity()};
      double trigger_factor{0.0};

      struct hybrid_dimer_t {
        bool use_mmf{false};
        double trigger_force{0.1};
        long max_steps{1000};
        long ci_stability_count{5};
        double angle_tol{0.7071}; // 1/sqrt(2): Householder stability bound
        double trigger_factor{0.0};
      } ocineb;
    } climbing_image;

    struct path_initialization_t {
      NEBInit method{NEBInit::LINEAR};
      std::string input_path;
      int max_iterations{5000};
      int nsteps{250};
      double max_move{0.1};
      double force_tolerance{0.001};
      double sidpp_alpha{0.33};
      double sidpp_frontier_tol{0.01}; ///< Convergence tol before adding images
      bool sidpp_reparam{true};        ///< Reparameterize after growth complete
      bool sidpp_ideal_ksp{false};     ///< Scale spring constant during growth
      OptType opt_method{OptType::LBFGS};
      bool oversampling{false};
      int oversampling_factor{3};
    } initialization;

    struct endpoint_options_t {
      bool minimize{true};
      bool use_path_file{false};
    } endpoints;

  } neb_options;

  // [Molecular Dynamics] //
  struct dynamics_options_t {
    double time_step_input{1.0};
    double time_step{0.0}; // computed: input / timeUnit
    double time_input{1000.0};
    double time{0.0}; // computed: input / timeUnit
    long steps{0};    // computed: time / time_step
  } dynamics_options;

  // [Parallel Replica] //
  struct parallel_replica_options_t {
    bool refine_transition{true};
    bool auto_stop{false};
    bool dephase_loop_stop{false};
    double dephase_time_input{1000.0};
    double dephase_time{0.0}; // computed: input / timeUnit
    long dephase_loop_max{5};
    double state_check_interval_input{1000.0};
    double state_check_interval{0.0}; // computed
    double record_interval_input{50.0};
    double record_interval{0.0}; // computed
    double corr_time_input{1000.0};
    double corr_time{0.0}; // computed
  } parallel_replica_options;

  // [Temperature Accelerated Dynamics] //
  struct tad_options_t {
    double low_temperature{300.0};
    double min_prefactor{0.001};
    double confidence{0.001};
  } tad_options;

  // [Thermostat] //
  struct thermostat_options_t {
    std::string kind{"none"};
    double andersen_alpha{1.0};
    double andersen_tcol_input{100.0};
    double andersen_tcol{0.0}; // computed
    double nose_mass{1.0};
    double langevin_friction_input{0.01};
    double langevin_friction{0.0}; // computed: input * timeUnit
  } thermostat_options;

  // [Replica Exchange] //
  struct replica_exchange_options_t {
    std::string temperature_distribution{"exponential"};
    long replicas{10};
    long exchange_trials{10}; // resolved to replicas
    double sampling_time_input{1000.0};
    double sampling_time{0.0}; // computed
    double temperature_high{0.0};
    double temperature_low{0.0};
    double exchange_period_input{100.0};
    double exchange_period{0.0}; // computed
  } replica_exchange_options;

  // [Bond Boost / Hyperdynamics] //
  struct hyperdynamics_options_t {
    std::string bias_potential{"none"};
    std::string boost_atom_list{"All"};
    double rmd_time_input{100.0};
    double rmd_time{0.0}; // computed
    double dvmax{0.0};
    double qrr{0.2};
    double prr{0.95};
    double qcut{3.0};
  } hyperdynamics_options;

  // [Basin Hopping] //
  struct basin_hopping_options_t {
    double displacement{0.5};
    double initial_random_structure_probability{0.0};
    double push_apart_distance{0.4};
    long steps{10000};
    long quenching_steps{0};
    bool significant_structure{true};
    bool single_atom_displace{false};
    std::string displacement_algorithm{"standard"};
    std::string displacement_distribution{"uniform"};
    double swap_probability{0.0};
    long jump_max{10};
    long jump_steps{0};
    bool adjust_displacement{true};
    long adjust_period{10};
    double adjust_fraction{0.05};
    double target_ratio{0.5};
    bool write_unique{false};
    double stop_energy{-std::numeric_limits<double>::max()};
  } basin_hopping_options;

  // [Global Optimization] //
  struct global_optimization_options_t {
    std::string move_method{"md"};
    std::string decision_method{"npew"};
    long steps{10000};
    double beta{1.05};
    double alpha{1.02};
    long mdmin{3};
    double target_energy{-1.0e50};
  } global_optimization_options;

  // [Monte Carlo] //
  struct monte_carlo_options_t {
    double step_size{0.005};
    int steps{1000};
  } monte_carlo_options;

  // [BGSD] //
  struct bgsd_options_t {
    double alpha{10.0};
    double beta{0.2};
    double gradient_finite_difference{0.000001};
    double h_force_convergence{0.01};
    double grad2energy_convergence{0.000001};
    double grad2force_convergence{0.0001};
  } bgsd_options;

  // [Serve] //
  struct serve_options_t {
    std::string host{"localhost"};
    uint16_t port{12345};
    size_t replicas{1};
    uint16_t gateway_port{0};
    std::string endpoints;
  } serve_options;

  // [ARTn] //
  struct artn_options_t {
    double push_step_size{0.3};   // maps to pARTn "push_step_size"
    double force_threshold{0.05}; // maps to pARTn "forc_thr"
    int max_iterations{500};
    int ninit{-1}; // maps to pARTn "ninit" (-1 = use pARTn default): number
                   // of initial push steps before Lanczos eigenmode
                   // estimation. 0 skips the push (eOn supplies the mode);
                   // larger values let pARTn explore further from the
                   // minimum before switching to Lanczos.
    std::string nperp_limitation{"default"}; // maps to pARTn "nperp_limitation"
    int lanczos_min_size{
        -1};         // maps to pARTn "lanczos_min_size" (-1 = use default)
    int nsmooth{-1}; // maps to pARTn "nsmooth" (-1 = use default)
    int nnewchance{
        3}; // maps to pARTn "nnewchance": retries when lowest eigval > 0
            // (convex region). pARTn default is 0 (immediate fail), which
            // is too aggressive for small clusters where Lanczos can
            // momentarily settle on a positive eigenvalue before the
            // unstable mode emerges. 3 retries is a reasonable default.
    std::string filin{""}; // maps to pARTn "filin": path to an artn.in input
                           // file. Empty preserves pARTn's post-create default
                           // (NAN_STR sentinel "BBBB" set by artn_create in
                           // artn_api.f90), which tells m_setup_artn to skip
                           // reading any file. If set, pARTn opens the file
                           // with status="old" and errors out when absent.
  } artn_options;

  // [IRA] //
  struct ira_options_t {
    double distance_threshold{0.3};
    double symmetry_threshold{0.1};
    bool use_pbc{false};
  } ira_options;

  // [Debug] //
  struct debug_options_t {
    bool write_movies{false};
    long write_movies_interval{1};
    bool write_deprecated_outs{false};
    bool estimate_neb_eigenvalues{false};
    std::string neb_mmf{"dimer"};
  } debug_options;
};

} // namespace eonc

using eonc::Parameters;

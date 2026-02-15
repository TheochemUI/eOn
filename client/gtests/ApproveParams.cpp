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
#include "ApprovalTests.hpp"
#include "Parameters.h"
#include "catch2/catch_amalgamated.hpp"
#include <iomanip>
#include <iostream>
#include <string>

std::ostream &operator<<(std::ostream &os, const Parameters &params) {
  os << std::setprecision(18);
  os << "Current Parameters are: " << std::endl;

  os << "\n[Main]" << std::endl;
  os << "kB: " << params.constants.kB << std::endl;
  os << "timeUnit: " << params.constants.timeUnit << std::endl;
  os << "job: " << magic_enum::enum_name(params.main_options.job) << std::endl;
  os << "randomSeed: " << params.main_options.randomSeed << std::endl;
  os << "temperature: " << params.main_options.temperature << std::endl;
  os << "quiet: " << std::boolalpha << params.main_options.quiet << std::endl;
  os << "writeLog: " << std::boolalpha << params.main_options.writeLog
     << std::endl;
  os << "checkpoint: " << std::boolalpha << params.main_options.checkpoint
     << std::endl;
  os << "iniFilename: " << params.main_options.iniFilename << std::endl;
  os << "conFilename: " << params.main_options.conFilename << std::endl;
  os << "finiteDifference: " << params.main_options.finiteDifference
     << std::endl;
  os << "maxForceCalls: " << params.main_options.maxForceCalls << std::endl;
  os << "removeNetForce: " << std::boolalpha
     << params.main_options.removeNetForce << std::endl;
  os << "\n[Potential]" << std::endl;
  os << "potential: "
     << magic_enum::enum_name(params.potential_options.potential) << std::endl;
  os << "MPIPollPeriod: " << params.potential_options.MPIPollPeriod
     << std::endl;
  os << "MPIPotentialRank: " << params.potential_options.MPIPotentialRank
     << std::endl;
  os << "LAMMPSLogging: " << std::boolalpha
     << params.potential_options.LAMMPSLogging << std::endl;
  os << "LAMMPSThreads: " << params.potential_options.LAMMPSThreads
     << std::endl;
  os << "EMTRasmussen: " << std::boolalpha
     << params.potential_options.EMTRasmussen << std::endl;
  os << "LogPotential: " << std::boolalpha
     << params.potential_options.LogPotential << std::endl;
  os << "extPotPath: " << params.potential_options.extPotPath << std::endl;

  os << "\n[AMS]" << std::endl;
  os << "engine: " << params.ams_options.engine << std::endl;
  os << "forcefield: " << params.ams_options.forcefield << std::endl;
  os << "model: " << params.ams_options.model << std::endl;
  os << "resources: " << params.ams_options.resources << std::endl;
  os << "xc: " << params.ams_options.xc << std::endl;
  os << "basis: " << params.ams_options.basis << std::endl;

  os << "\n[AMS_ENV]" << std::endl;
  os << "amshome: " << params.ams_options.env.amshome << std::endl;
  os << "scm_tmpdir: " << params.ams_options.env.scm_tmpdir << std::endl;
  os << "scmlicense: " << params.ams_options.env.scmlicense << std::endl;
  os << "scm_pythondir: " << params.ams_options.env.scm_pythondir << std::endl;
  os << "amsbin: " << params.ams_options.env.amsbin << std::endl;
  os << "amsresources: " << params.ams_options.env.amsresources << std::endl;

  os << "\n[XTBPot]" << std::endl;
  os << "xtb_paramset: " << params.xtb_options.paramset << std::endl;
  os << "xtb_elec_temperature: " << params.xtb_options.elec_temperature
     << std::endl;
  os << "xtb_maxiter: " << params.xtb_options.maxiter << std::endl;
  os << "xtb_acc: " << params.xtb_options.acc << std::endl;

  os << "\n[Structure Comparison]" << std::endl;
  os << "distanceDifference: "
     << params.structure_comparison_options.distance_difference << std::endl;
  os << "neighborCutoff: "
     << params.structure_comparison_options.neighbor_cutoff << std::endl;
  os << "checkRotation: " << std::boolalpha
     << params.structure_comparison_options.check_rotation << std::endl;
  os << "indistinguishableAtoms: " << std::boolalpha
     << params.structure_comparison_options.indistinguishable_atoms
     << std::endl;
  os << "energyDifference: "
     << params.structure_comparison_options.energy_difference << std::endl;
  os << "removeTranslation: " << std::boolalpha
     << params.structure_comparison_options.remove_translation << std::endl;

  os << "\n[Process Search]" << std::endl;
  os << "processSearchMinimizeFirst: " << std::boolalpha
     << params.process_search_options.minimize_first << std::endl;
  os << "processSearchMinimizationOffset: "
     << params.process_search_options.minimization_offset << std::endl;

  os << "\n[Saddle Search]" << std::endl;
  os << "saddleMaxJumpAttempts: "
     << params.saddle_search_options.max_jump_attempts << std::endl;
  os << "saddleMaxIterations: " << params.saddle_search_options.max_iterations
     << std::endl;
  os << "saddleMethod: " << params.saddle_search_options.method << std::endl;
  os << "saddleMinmodeMethod: " << params.saddle_search_options.minmode_method
     << std::endl;
  os << "saddleDisplaceType: " << params.saddle_search_options.displace_type
     << std::endl;
  os << "saddleMaxEnergy: " << params.saddle_search_options.max_energy
     << std::endl;
  os << "saddleDisplaceMagnitude: "
     << params.saddle_search_options.displace_magnitude << std::endl;
  os << "saddleMaxSingleDisplace: "
     << params.saddle_search_options.max_single_displace << std::endl;
  os << "saddleDisplaceRadius: " << params.saddle_search_options.displace_radius
     << std::endl;
  os << "saddleConvergedForce: " << params.saddle_search_options.converged_force
     << std::endl;
  os << "saddlePerpForceRatio: "
     << params.saddle_search_options.perp_force_ratio << std::endl;
  os << "saddleNonnegativeDisplacementAbort: " << std::boolalpha
     << params.saddle_search_options.nonnegative_displacement_abort
     << std::endl;
  os << "saddleNonlocalCountAbort: "
     << params.saddle_search_options.nonlocal_count_abort << std::endl;
  os << "saddleNonlocalDistanceAbort: "
     << params.saddle_search_options.nonlocal_distance_abort << std::endl;
  os << "saddleRemoveRotation: " << std::boolalpha
     << params.saddle_search_options.remove_rotation << std::endl;
  os << "saddleDynamicsTemperature: "
     << params.saddle_search_options.dynamics.temperature << std::endl;
  os << "saddleDynamicsStateCheckIntervalInput: "
     << params.saddle_search_options.dynamics.state_check_interval_input
     << std::endl;
  os << "saddleDynamicsStateCheckInterval: "
     << params.saddle_search_options.dynamics.state_check_interval << std::endl;
  os << "saddleDynamicsRecordIntervalInput: "
     << params.saddle_search_options.dynamics.record_interval_input
     << std::endl;
  os << "saddleDynamicsRecordInterval: "
     << params.saddle_search_options.dynamics.record_interval << std::endl;
  os << "saddleDynamicsLinearInterpolation: " << std::boolalpha
     << params.saddle_search_options.dynamics.linear_interpolation << std::endl;
  os << "saddleDynamicsMaxInitCurvature: "
     << params.saddle_search_options.dynamics.max_init_curvature << std::endl;
  os << "saddleConfinePositive: " << std::boolalpha
     << params.saddle_search_options.confine_positive.enabled << std::endl;
  os << "saddleBowlBreakout: " << std::boolalpha
     << params.saddle_search_options.confine_positive.bowl_breakout
     << std::endl;
  os << "saddleBowlActive: "
     << params.saddle_search_options.confine_positive.bowl_active << std::endl;
  os << "saddleConfinePositiveMinForce: "
     << params.saddle_search_options.confine_positive.min_force << std::endl;
  os << "saddleConfinePositiveScaleRatio: "
     << params.saddle_search_options.confine_positive.scale_ratio << std::endl;
  os << "saddleConfinePositiveBoost: "
     << params.saddle_search_options.confine_positive.boost << std::endl;
  os << "saddleConfinePositiveMinActive: "
     << params.saddle_search_options.confine_positive.min_active << std::endl;
  os << "saddleZeroModeAbortCurvature: "
     << params.saddle_search_options.zero_mode_abort_curvature << std::endl;

  os << "\n[Optimizers]" << std::endl;
  os << "optMethod: " << magic_enum::enum_name(params.optimizer_options.method)
     << std::endl;
  os << "optConvergenceMetric: " << params.optimizer_options.convergence_metric
     << std::endl;
  os << "optConvergenceMetricLabel: "
     << params.optimizer_options.convergence_metric_label << std::endl;
  os << "optMaxIterations: " << params.optimizer_options.max_iterations
     << std::endl;
  os << "optMaxMove: " << params.optimizer_options.max_move << std::endl;
  os << "optConvergedForce: " << params.optimizer_options.converged_force
     << std::endl;
  os << "optTimeStepInput: " << params.optimizer_options.time_step_input
     << std::endl;
  os << "optTimeStep: " << params.optimizer_options.time_step << std::endl;
  os << "optMaxTimeStepInput: " << params.optimizer_options.max_time_step_input
     << std::endl;
  os << "optMaxTimeStep: " << params.optimizer_options.max_time_step
     << std::endl;
  os << "optLBFGSMemory: " << params.optimizer_options.lbfgs.memory
     << std::endl;
  os << "optLBFGSInverseCurvature: "
     << params.optimizer_options.lbfgs.inverse_curvature << std::endl;
  os << "optLBFGSMaxInverseCurvature: "
     << params.optimizer_options.lbfgs.max_inverse_curvature << std::endl;
  os << "optLBFGSAutoScale: " << std::boolalpha
     << params.optimizer_options.lbfgs.auto_scale << std::endl;
  os << "optLBFGSAngleReset: " << std::boolalpha
     << params.optimizer_options.lbfgs.angle_reset << std::endl;
  os << "optLBFGSDistanceReset: " << std::boolalpha
     << params.optimizer_options.lbfgs.distance_reset << std::endl;
  os << "optQMSteepestDecent: " << std::boolalpha
     << params.optimizer_options.quickmin.steepest_descent << std::endl;
  os << "optCGNoOvershooting: " << std::boolalpha
     << params.optimizer_options.cg.no_overshooting << std::endl;
  os << "optCGKnockOutMaxMove: " << std::boolalpha
     << params.optimizer_options.cg.knock_out_max_move << std::endl;
  os << "optCGLineSearch: " << std::boolalpha
     << params.optimizer_options.cg.line_search << std::endl;
  os << "optCGLineConverged: " << params.optimizer_options.cg.line_converged
     << std::endl;
  os << "optCGLineSearchMaxIter: "
     << params.optimizer_options.cg.line_search_max_iter << std::endl;
  os << "optCGMaxIterBeforeReset: "
     << params.optimizer_options.cg.max_iter_before_reset << std::endl;
  os << "optSDAlpha: " << params.optimizer_options.sd.alpha << std::endl;
  os << "optSDTwoPoint: " << std::boolalpha
     << params.optimizer_options.sd.two_point << std::endl;

  os << "\n[Refine]" << std::endl;
  os << "refineOptMethod: "
     << magic_enum::enum_name(params.optimizer_options.refine.method)
     << std::endl;
  os << "refineThreshold: " << params.optimizer_options.refine.threshold
     << std::endl;

  os << "\n[Dimer]" << std::endl;
  os << "dimerRotationAngle: " << params.dimer_options.rotation_angle
     << std::endl;
  os << "dimerImproved: " << std::boolalpha << params.dimer_options.improved
     << std::endl;
  os << "dimerConvergedAngle: " << params.dimer_options.converged_angle
     << std::endl;
  os << "dimerMaxIterations: " << params.dimer_options.max_iterations
     << std::endl;
  os << "dimerOptMethod: " << params.dimer_options.opt_method << std::endl;
  os << "dimerRotationsMax: " << params.dimer_options.rotations_max
     << std::endl;
  os << "dimerRotationsMin: " << params.dimer_options.rotations_min
     << std::endl;
  os << "dimerTorqueMax: " << params.dimer_options.torque_max << std::endl;
  os << "dimerTorqueMin: " << params.dimer_options.torque_min << std::endl;
  os << "dimerRemoveRotation: " << std::boolalpha
     << params.dimer_options.remove_rotation << std::endl;

  os << "\n[GPR Dimer]" << std::endl;
  os << "gprDimerRotationAngle: " << params.gpr_dimer_options.rotation_angle
     << std::endl;
  os << "gprDimerConvergedAngle: " << params.gpr_dimer_options.converged_angle
     << std::endl;
  os << "gprDimerRelaxConvAngle: " << params.gpr_dimer_options.relax_conv_angle
     << std::endl;
  os << "gprDimerInitRotationsMax: "
     << params.gpr_dimer_options.init_rotations_max << std::endl;
  os << "gprDimerRelaxRotationsMax: "
     << params.gpr_dimer_options.relax_rotations_max << std::endl;
  os << "gprDimerDivisorTdimerGP: "
     << params.gpr_dimer_options.divisor_t_dimer_gp << std::endl;
  os << "gprDimerMaxOuterIterations: "
     << params.gpr_dimer_options.max_outer_iterations << std::endl;
  os << "gprDimerMaxInnerIterations: "
     << params.gpr_dimer_options.max_inner_iterations << std::endl;
  os << "gprDimerMidpointMaxDisp: "
     << params.gpr_dimer_options.midpoint_max_disp << std::endl;
  os << "gprDimerRotOptMethod: " << params.gpr_dimer_options.rot_opt_method
     << std::endl;
  os << "gprDimerTransOptMethod: " << params.gpr_dimer_options.trans_opt_method
     << std::endl;
  os << "gprActiveRadius: " << params.gpr_dimer_options.active_radius
     << std::endl;
  os << "gprDimerSep: " << params.gpr_dimer_options.dimer_sep << std::endl;
  os << "gprDimerConvStep: " << params.gpr_dimer_options.conv_step << std::endl;
  os << "gprDimerMaxStep: " << params.gpr_dimer_options.max_step << std::endl;
  os << "gprForceThreshold: " << params.saddle_search_options.converged_force
     << std::endl;
  os << "gprDimerRatioAtLimit: " << params.gpr_dimer_options.ratio_at_limit
     << std::endl;
  os << "gprDimerInitRotGP: " << std::boolalpha
     << params.gpr_dimer_options.init_rot_gp << std::endl;
  os << "gprDimerInitTransGP: " << std::boolalpha
     << params.gpr_dimer_options.init_trans_gp << std::endl;
  os << "gprDimerManyIterations: " << std::boolalpha
     << params.gpr_dimer_options.many_iterations << std::endl;
  os << "gprDimerHyperOptMethod: "
     << params.gpr_dimer_options.gpr_params.hyper_opt_method << std::endl;
  os << "gprDimerSigma2: " << params.gpr_dimer_options.gpr_params.sigma2
     << std::endl;
  os << "gprDimerJitterSigma2: "
     << params.gpr_dimer_options.gpr_params.jitter_sigma2 << std::endl;
  os << "gprDimerNoiseSigma2: "
     << params.gpr_dimer_options.gpr_params.noise_sigma2 << std::endl;
  os << "gprDimerPriorMu: " << params.gpr_dimer_options.gpr_params.prior_mu
     << std::endl;
  os << "gprDimerPriorSigma2: "
     << params.gpr_dimer_options.gpr_params.prior_sigma2 << std::endl;
  os << "gprDimerPriorNu: " << params.gpr_dimer_options.gpr_params.prior_nu
     << std::endl;
  os << "gprOptCheckDerivatives: " << std::boolalpha
     << params.gpr_dimer_options.opt_params.check_derivatives << std::endl;
  os << "gprOptMaxIterations: "
     << params.gpr_dimer_options.opt_params.max_iterations << std::endl;
  os << "gprOptTolFunc: " << params.gpr_dimer_options.opt_params.tol_func
     << std::endl;
  os << "gprOptTolSol: " << params.gpr_dimer_options.opt_params.tol_sol
     << std::endl;
  os << "gprOptLambdaLimit: "
     << params.gpr_dimer_options.opt_params.lambda_limit << std::endl;
  os << "gprOptLambdaInit: " << params.gpr_dimer_options.opt_params.lambda_init
     << std::endl;
  os << "gprUsePrune: " << std::boolalpha
     << params.gpr_dimer_options.prune_params.use_prune << std::endl;
  os << "gprPruneBegin: " << params.gpr_dimer_options.prune_params.begin
     << std::endl;
  os << "gprPruneNVals: " << params.gpr_dimer_options.prune_params.n_vals
     << std::endl;
  os << "gprPruneThreshold: " << params.gpr_dimer_options.prune_params.threshold
     << std::endl;
  os << "gprReportLevel: " << params.gpr_dimer_options.debug_params.report_level
     << std::endl;
  os << "gprDebugLevel: " << params.gpr_dimer_options.debug_params.debug_level
     << std::endl;
  os << "gprDebugOutDir: " << params.gpr_dimer_options.debug_params.out_dir
     << std::endl;
  os << "gprDebugPosFile: " << params.gpr_dimer_options.debug_params.pos_file
     << std::endl;
  os << "gprDebugEnergyFile: "
     << params.gpr_dimer_options.debug_params.energy_file << std::endl;
  os << "gprDebugGradFile: " << params.gpr_dimer_options.debug_params.grad_file
     << std::endl;
  os << "gprDebugOutExt: " << params.gpr_dimer_options.debug_params.out_ext
     << std::endl;
  os << "gprDebugOffsetMidPoint: "
     << params.gpr_dimer_options.debug_params.offset_mid_point << std::endl;
  os << "gprDebugDy: " << params.gpr_dimer_options.debug_params.dy << std::endl;
  os << "gprDebugDz: " << params.gpr_dimer_options.debug_params.dz << std::endl;

  os << "\n[Surrogate]" << std::endl;
  os << "use_surrogate: " << params.gp_surrogate_options.enabled << std::endl;
  os << "sub_job: "
     << magic_enum::enum_name(params.gp_surrogate_options.sub_job) << std::endl;
  os << "gp_uncertainty: " << params.gp_surrogate_options.uncertainty
     << std::endl;
  os << "gp_linear_path_always: " << std::boolalpha
     << params.gp_surrogate_options.linear_path_always << std::endl;
  os << "surrogatePotential: "
     << magic_enum::enum_name(params.gp_surrogate_options.potential)
     << std::endl;

  os << "\n[CatLearn]" << std::endl;
  os << "catl_path: " << params.catlearn_options.path << std::endl;
  os << "catl_model: " << params.catlearn_options.model << std::endl;
  os << "catl_prior: " << params.catlearn_options.prior << std::endl;
  os << "catl_use_deriv: " << std::boolalpha
     << params.catlearn_options.use_deriv << std::endl;
  os << "catl_use_fingerprint: " << std::boolalpha
     << params.catlearn_options.use_fingerprint << std::endl;
  os << "catl_parallel: " << std::boolalpha << params.catlearn_options.parallel
     << std::endl;

  os << "\n[ASE ORCA]" << std::endl;
  os << "orca_path: " << params.ase_orca_options.path << std::endl;
  os << "orca_nproc: " << params.ase_orca_options.nproc << std::endl;
  os << "orca_sline: " << params.ase_orca_options.simpleinput << std::endl;

  os << "\n[Lanczos]" << std::endl;
  os << "lanczosTolerance: " << params.lanczos_options.tolerance << std::endl;
  os << "lanczosMaxIterations: " << params.lanczos_options.max_iterations
     << std::endl;
  os << "lanczosQuitEarly: " << std::boolalpha
     << params.lanczos_options.quit_early << std::endl;

  os << "\n[Prefactor]" << std::endl;
  os << "prefactorDefaultValue: " << params.prefactor_options.default_value
     << std::endl;
  os << "prefactorMaxValue: " << params.prefactor_options.max_value
     << std::endl;
  os << "prefactorMinValue: " << params.prefactor_options.min_value
     << std::endl;
  os << "prefactorWithinRadius: " << params.prefactor_options.within_radius
     << std::endl;
  os << "prefactorMinDisplacement: "
     << params.prefactor_options.min_displacement << std::endl;
  os << "prefactorRate: " << params.prefactor_options.rate << std::endl;
  os << "prefactorConfiguration: " << params.prefactor_options.configuration
     << std::endl;
  os << "prefactorAllFreeAtoms: " << std::boolalpha
     << params.prefactor_options.all_free_atoms << std::endl;
  os << "prefactorFilterScheme: " << params.prefactor_options.filter_scheme
     << std::endl;
  os << "prefactorFilterFraction: " << params.prefactor_options.filter_fraction
     << std::endl;

  os << "\n[Hessian]" << std::endl;
  os << "hessianAtomList: " << params.hessian_options.atom_list << std::endl;
  os << "hessianZeroFreqValue: " << params.hessian_options.zero_freq_value
     << std::endl;

  os << "\n[Nudged Elastic Band]" << std::endl;
  os << "nebImages: " << params.nebImages << std::endl;
  os << "nebMaxIterations: " << params.nebMaxIterations << std::endl;
  os << "nebSpring: " << params.nebSpring << std::endl;
  os << "nebClimbingImageMethod: " << std::boolalpha
     << params.nebClimbingImageMethod << std::endl;
  os << "nebClimbingImageConvergedOnly: " << std::boolalpha
     << params.nebClimbingImageConvergedOnly << std::endl;
  os << "nebOldTangent: " << std::boolalpha << params.nebOldTangent
     << std::endl;
  os << "nebDoublyNudged: " << std::boolalpha << params.nebDoublyNudged
     << std::endl;
  os << "nebDoublyNudgedSwitching: " << std::boolalpha
     << params.nebDoublyNudgedSwitching << std::endl;
  os << "nebOptMethod: " << params.nebOptMethod << std::endl;
  os << "nebElasticBand: " << std::boolalpha << params.nebElasticBand
     << std::endl;
  os << "nebConvergedForce: " << params.nebConvergedForce << std::endl;
  os << "nebKSPMin: " << params.nebKSPMin << std::endl;
  os << "nebKSPMax: " << params.nebKSPMax << std::endl;
  os << "nebEnergyWeighted: " << std::boolalpha << params.nebEnergyWeighted
     << std::endl;

  os << "\n[Dynamics]" << std::endl;
  os << "mdTimeStepInput: " << params.dynamics_options.time_step_input
     << std::endl;
  os << "mdTimeStep: " << params.dynamics_options.time_step << std::endl;
  os << "mdTimeInput: " << params.dynamics_options.time_input << std::endl;
  os << "mdTime: " << params.dynamics_options.time << std::endl;
  os << "mdSteps: " << params.dynamics_options.steps << std::endl;

  os << "\n[Parallel Replica]" << std::endl;
  os << "parrepRefineTransition: " << std::boolalpha
     << params.parallel_replica_options.refine_transition << std::endl;
  os << "parrepAutoStop: " << std::boolalpha
     << params.parallel_replica_options.auto_stop << std::endl;
  os << "parrepDephaseLoopStop: " << std::boolalpha
     << params.parallel_replica_options.dephase_loop_stop << std::endl;
  os << "parrepDephaseTimeInput: "
     << params.parallel_replica_options.dephase_time_input << std::endl;
  os << "parrepDephaseTime: " << params.parallel_replica_options.dephase_time
     << std::endl;
  os << "parrepDephaseLoopMax: "
     << params.parallel_replica_options.dephase_loop_max << std::endl;
  os << "parrepStateCheckIntervalInput: "
     << params.parallel_replica_options.state_check_interval_input << std::endl;
  os << "parrepStateCheckInterval: "
     << params.parallel_replica_options.state_check_interval << std::endl;
  os << "parrepRecordIntervalInput: "
     << params.parallel_replica_options.record_interval_input << std::endl;
  os << "parrepRecordInterval: "
     << params.parallel_replica_options.record_interval << std::endl;
  os << "parrepCorrTimeInput: "
     << params.parallel_replica_options.corr_time_input << std::endl;
  os << "parrepCorrTime: " << params.parallel_replica_options.corr_time
     << std::endl;

  os << "\n[TAD]" << std::endl;
  os << "tadLowT: " << params.tadLowT << std::endl;
  os << "tadMinPrefactor: " << params.tadMinPrefactor << std::endl;
  os << "tadConfidence: " << params.tadConfidence << std::endl;

  os << "\n[Thermostat]" << std::endl;
  os << "thermostat: " << params.thermostat_options.kind << std::endl;
  os << "thermoAndersenAlpha: " << params.thermostat_options.andersen_alpha
     << std::endl;
  os << "thermoAndersenTcolInput: "
     << params.thermostat_options.andersen_tcol_input << std::endl;
  os << "thermoAndersenTcol: " << params.thermostat_options.andersen_tcol
     << std::endl;
  os << "thermoNoseMass: " << params.thermostat_options.nose_mass << std::endl;
  os << "thermoLangevinFrictionInput: "
     << params.thermostat_options.langevin_friction_input << std::endl;
  os << "thermoLangevinFriction: "
     << params.thermostat_options.langevin_friction << std::endl;

  os << "\n[Replica Exchange]" << std::endl;
  os << "repexcTemperatureDistribution: "
     << params.replica_exchange_options.temperature_distribution << std::endl;
  os << "repexcReplicas: " << params.replica_exchange_options.replicas
     << std::endl;
  os << "repexcExchangeTrials: "
     << params.replica_exchange_options.exchange_trials << std::endl;
  os << "repexcSamplingTimeInput: "
     << params.replica_exchange_options.sampling_time_input << std::endl;
  os << "repexcSamplingTime: " << params.replica_exchange_options.sampling_time
     << std::endl;
  os << "repexcTemperatureHigh: "
     << params.replica_exchange_options.temperature_high << std::endl;
  os << "repexcTemperatureLow: "
     << params.replica_exchange_options.temperature_low << std::endl;
  os << "repexcExchangePeriodInput: "
     << params.replica_exchange_options.exchange_period_input << std::endl;
  os << "repexcExchangePeriod: "
     << params.replica_exchange_options.exchange_period << std::endl;

  os << "\n[Hyperdynamics]" << std::endl;
  os << "biasPotential: " << params.hyperdynamics_options.bias_potential
     << std::endl;
  os << "bondBoostBALS: " << params.hyperdynamics_options.boost_atom_list
     << std::endl;
  os << "bondBoostRMDTimeInput: " << params.hyperdynamics_options.rmd_time_input
     << std::endl;
  os << "bondBoostRMDTime: " << params.hyperdynamics_options.rmd_time
     << std::endl;
  os << "bondBoostDVMAX: " << params.hyperdynamics_options.dvmax << std::endl;
  os << "bondBoostQRR: " << params.hyperdynamics_options.qrr << std::endl;
  os << "bondBoostPRR: " << params.hyperdynamics_options.prr << std::endl;
  os << "bondBoostQcut: " << params.hyperdynamics_options.qcut << std::endl;
  os << "basinHoppingDisplacement: "
     << params.basin_hopping_options.displacement << std::endl;
  os << "basinHoppingInitialRandomStructureProbability: "
     << params.basin_hopping_options.initial_random_structure_probability
     << std::endl;
  os << "basinHoppingPushApartDistance: "
     << params.basin_hopping_options.push_apart_distance << std::endl;
  os << "basinHoppingSteps: " << params.basin_hopping_options.steps
     << std::endl;
  os << "basinHoppingQuenchingSteps: "
     << params.basin_hopping_options.quenching_steps << std::endl;
  os << "basinHoppingSignificantStructure: " << std::boolalpha
     << params.basin_hopping_options.significant_structure << std::endl;
  os << "basinHoppingSingleAtomDisplace: " << std::boolalpha
     << params.basin_hopping_options.single_atom_displace << std::endl;
  os << "basinHoppingDisplacementAlgorithm: "
     << params.basin_hopping_options.displacement_algorithm << std::endl;
  os << "basinHoppingDisplacementDistribution: "
     << params.basin_hopping_options.displacement_distribution << std::endl;
  os << "basinHoppingSwapProbability: "
     << params.basin_hopping_options.swap_probability << std::endl;
  os << "basinHoppingJumpMax: " << params.basin_hopping_options.jump_max
     << std::endl;
  os << "basinHoppingJumpSteps: " << params.basin_hopping_options.jump_steps
     << std::endl;
  os << "basinHoppingAdjustDisplacement: " << std::boolalpha
     << params.basin_hopping_options.adjust_displacement << std::endl;
  os << "basinHoppingAdjustPeriod: "
     << params.basin_hopping_options.adjust_period << std::endl;
  os << "basinHoppingAdjustFraction: "
     << params.basin_hopping_options.adjust_fraction << std::endl;
  os << "basinHoppingTargetRatio: " << params.basin_hopping_options.target_ratio
     << std::endl;
  os << "basinHoppingWriteUnique: " << std::boolalpha
     << params.basin_hopping_options.write_unique << std::endl;
  os << "basinHoppingStopEnergy (-DBL_MAX): "
     << params.basin_hopping_options.stop_energy << std::endl;

  os << "\n[Global Optimization]" << std::endl;
  os << "globalOptimizationMoveMethod: "
     << params.global_optimization_options.move_method << std::endl;
  os << "globalOptimizationDecisionMethod: "
     << params.global_optimization_options.decision_method << std::endl;
  os << "globalOptimizationSteps: " << params.global_optimization_options.steps
     << std::endl;
  os << "globalOptimizationBeta: " << params.global_optimization_options.beta
     << std::endl;
  os << "globalOptimizationAlpha: " << params.global_optimization_options.alpha
     << std::endl;
  os << "globalOptimizationMdmin: " << params.global_optimization_options.mdmin
     << std::endl;
  os << "globalOptimizationTargetEnergy: "
     << params.global_optimization_options.target_energy << std::endl;

  os << "\n[Monte Carlo]" << std::endl;
  os << "monteCarloStepSize: " << params.monte_carlo_options.step_size
     << std::endl;
  os << "monteCarloSteps: " << params.monte_carlo_options.steps << std::endl;

  os << "\n[BGSD]" << std::endl;
  os << "alpha: " << params.bgsd_options.alpha << std::endl;
  os << "beta: " << params.bgsd_options.beta << std::endl;
  os << "gradientfinitedifference: "
     << params.bgsd_options.gradient_finite_difference << std::endl;
  os << "Hforceconvergence: " << params.bgsd_options.h_force_convergence
     << std::endl;
  os << "grad2energyconvergence: "
     << params.bgsd_options.grad2energy_convergence << std::endl;
  os << "grad2forceconvergence: " << params.bgsd_options.grad2force_convergence
     << std::endl;

  os << "\n[Debug]" << std::endl;
  os << "writeMovies: " << std::boolalpha << params.debug_options.write_movies
     << std::endl;
  os << "writeMoviesInterval: " << params.debug_options.write_movies_interval
     << std::endl;
  return os;
}

std::vector<Parameters> getParameters() {
  // Return test data for Parameters
  return {Parameters()};
}

TEST_CASE("VerifyParameters") {
  ApprovalTests::Approvals::verifyAll("parameters", getParameters());
}

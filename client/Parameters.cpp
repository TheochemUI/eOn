#include "Parameters.h"

#include "BondBoost.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "Hessian.h"
#include "INIFile.h"
#include "Log.h"
#include "NudgedElasticBand.h"
#include "ReplicaExchangeJob.h"

#include <errno.h>
#include <float.h>
#include <string>
#include <time.h>

std::string Parameters::toLowerCase(std::string s) {
    for (std::string::size_type i = 0; i < s.length(); ++i) {
        s[i] = tolower(s[i]);
    }
    return s;
}

int Parameters::load(string filename)
{
    FILE *fh;
    fh = fopen(filename.c_str(), "rb");
    if (fh == NULL) {
        std::cerr<<"error: "<<strerror(errno)<<"\n";
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

        job = toLowerCase(ini.GetValue("Main", "job"));
        temperature = ini.GetValueF("Main", "temperature", temperature);
        randomSeed = ini.GetValueL("Main", "random_seed", randomSeed);
        checkpoint = ini.GetValueB("Main", "checkpoint", checkpoint);
        quiet = ini.GetValueB("Main", "quiet", quiet);
        writeLog = ini.GetValueB("Main", "write_log", writeLog);
        finiteDifference = ini.GetValueF("Main", "finite_difference", finiteDifference);
        // Initialize random generator
        if (randomSeed < 0) {
            unsigned i = time(NULL);
            randomSeed = i;
            helper_functions::random(i);
        } else {
            helper_functions::random(randomSeed);
        }
        maxForceCalls = ini.GetValueL("Main", "max_force_calls", maxForceCalls);
        removeNetForce = ini.GetValueB("Main", "remove_net_force", removeNetForce);

        // [Potential] //

        potential = toLowerCase(ini.GetValue("Potential", "potential"));
        MPIPollPeriod = ini.GetValueF("Potential", "mpi_poll_period", MPIPollPeriod);
        LAMMPSLogging = ini.GetValueB("Potential", "lammps_logging", LAMMPSLogging);
        LAMMPSThreads = (int) ini.GetValueL("Potential", "lammps_threads", LAMMPSThreads);
        EMTRasmussen = ini.GetValueB("Potential", "emt_rasmussen", EMTRasmussen);
        extPotPath = ini.GetValue("Potential", "ext_pot_path", extPotPath);
        if (potential == "mpi" || potential == "vasp" || potential == "bopfox"
            || potential == "bop") {
            LogPotential = true;
        } else {
            LogPotential = false;
        }
        LogPotential = ini.GetValueB("Potential", "log_potential", LogPotential);

        // [AMS]
        if (potential == "ams") {
            engine = ini.GetValue("AMS", "engine", engine);
            forcefield = ini.GetValue("AMS", "forcefield", forcefield);
            resources = ini.GetValue("AMS", "resources", resources);
            model = ini.GetValue("AMS", "model", model);
            xc = ini.GetValue("AMS", "xc", xc);
            basis = ini.GetValue("AMS", "basis", basis);
        }
        // [AMS_IO]
        if (potential == "ams_io") {
            engine = ini.GetValue("AMS_IO", "engine", engine);
            forcefield = ini.GetValue("AMS_IO", "forcefield", forcefield);
            model = ini.GetValue("AMS_IO", "model", model);
            xc = ini.GetValue("AMS_IO", "xc", xc);
        }
        // [AMS_ENV]
        // This is only needed if the regular calls do not work
        // e.g. on a MacOS machine
        if (potential == "ams_io" || potential == "ams") {
            amshome = ini.GetValue("AMS_ENV", "amshome", amshome);
            scm_tmpdir = ini.GetValue("AMS_ENV", "scm_tmpdir", scm_tmpdir);
            scmlicense = ini.GetValue("AMS_ENV", "scmlicense", scmlicense);
            scm_pythondir = ini.GetValue("AMS_ENV", "scm_pythondir", scm_pythondir);
            amsbin = ini.GetValue("AMS_ENV", "amsbin", amsbin);
            amsresources = ini.GetValue("AMS_ENV", "amsresources", amsresources);
        }

        // [Debug] //

        writeMovies = ini.GetValueB("Debug", "write_movies", writeMovies);
        writeMoviesInterval = ini.GetValueL("Debug", "write_movies_interval", writeMoviesInterval);

        // [Structure Comparison] //

        distanceDifference
            = ini.GetValueF("Structure Comparison", "distance_difference", distanceDifference);
        neighborCutoff = ini.GetValueF("Structure Comparison", "neighbor_cutoff", neighborCutoff);
        checkRotation = ini.GetValueB("Structure Comparison", "check_rotation", checkRotation);
        energyDifference
            = ini.GetValueF("Structure Comparison", "energy_difference", energyDifference);
        indistinguishableAtoms = ini.GetValueB(
            "Structure Comparison", "indistinguishable_atoms", indistinguishableAtoms);
        removeTranslation
            = ini.GetValueB("Structure Comparison", "remove_translation", removeTranslation);

        // [Process Search] //

        processSearchMinimizeFirst
            = ini.GetValueB("Process Search", "minimize_first", processSearchMinimizeFirst);
        processSearchMinimizationOffset = ini.GetValueF(
            "Process Search", "minimization_offset", processSearchMinimizationOffset);

        // [Optimizers] //

        optMethod = toLowerCase(ini.GetValue("Optimizer", "opt_method", optMethod));
        optConvergenceMetric
            = toLowerCase(ini.GetValue("Optimizer", "convergence_metric", optConvergenceMetric));

        optConvergedForce = ini.GetValueF("Optimizer", "converged_force", optConvergedForce);
        optMaxIterations = ini.GetValueL("Optimizer", "max_iterations", optMaxIterations);
        optMaxMove = ini.GetValueF("Optimizer", "max_move", optMaxMove);
        processSearchMinimizationOffset = optMaxMove;
        optTimeStepInput = ini.GetValueF("Optimizer", "time_step", optTimeStepInput);
        optTimeStep = optTimeStepInput / timeUnit;
        optMaxTimeStepInput = ini.GetValueF("Optimizer", "time_step_max", optMaxTimeStepInput);
        optMaxTimeStep = optMaxTimeStepInput / timeUnit;
        optLBFGSMemory = ini.GetValueL("Optimizer", "lbfgs_memory", optLBFGSMemory);
        optLBFGSInverseCurvature
            = ini.GetValueF("Optimizer", "lbfgs_inverse_curvature", optLBFGSInverseCurvature);
        optLBFGSMaxInverseCurvature = ini.GetValueF(
            "Optimizer", "lbfgs_max_inverse_curvature", optLBFGSMaxInverseCurvature);
        optLBFGSAutoScale = ini.GetValueB("Optimizer", "lbfgs_auto_scale", optLBFGSAutoScale);
        optLBFGSAngleReset = ini.GetValueB("Optimizer", "lbfgs_angle_reset", optLBFGSAngleReset);
        optLBFGSDistanceReset
            = ini.GetValueB("Optimizer", "lbfgs_distance_reset", optLBFGSDistanceReset);
        optQMSteepestDecent
            = ini.GetValueB("Optimizer", "qm_steepest_descent", optQMSteepestDecent);
        optCGNoOvershooting
            = ini.GetValueB("Optimizer", "cg_no_overshooting", optCGNoOvershooting);
        optCGKnockOutMaxMove
            = ini.GetValueB("Optimizer", "cg_knock_out_max_move", optCGKnockOutMaxMove);
        optCGLineSearch = ini.GetValueB("Optimizer", "cg_line_search", optCGLineSearch);
        optCGLineConverged = ini.GetValueF("Optimizer", "cg_line_converged", optCGLineConverged);
        optCGMaxIterBeforeReset
            = ini.GetValueL("Optimizer", "cg_max_iter_before_reset", optCGMaxIterBeforeReset);
        optCGLineSearchMaxIter
            = ini.GetValueL("Optimizer", "cg_max_iter_line_search", optCGLineSearchMaxIter);
        optSDAlpha = ini.GetValueF("Optimizer", "sd_alpha", optSDAlpha);
        optSDTwoPoint = ini.GetValueB("Optimizer", "sd_twopoint", optSDTwoPoint);

        // [Dimer] //

        dimerRotationAngle = ini.GetValueF("Dimer", "finite_angle", dimerRotationAngle);
        dimerImproved = ini.GetValueB("Dimer", "improved", dimerImproved);
        dimerConvergedAngle = ini.GetValueF("Dimer", "converged_angle", dimerConvergedAngle);
        dimerMaxIterations = ini.GetValueL("Dimer", "max_iterations", dimerMaxIterations);
        dimerOptMethod = toLowerCase(ini.GetValue("Dimer", "opt_method", dimerOptMethod));
        dimerRotationsMin = ini.GetValueL("Dimer", "rotations_min", dimerRotationsMin); // old
        dimerRotationsMax
            = ini.GetValueL("Dimer", "rotations_max", dimerRotationsMax);      // old & new
        dimerTorqueMin = ini.GetValueF("Dimer", "torque_min", dimerTorqueMin); // old
        dimerTorqueMax = ini.GetValueF("Dimer", "torque_max", dimerTorqueMax); // old
        dimerRemoveRotation = ini.GetValueB("Dimer", "remove_rotation", dimerRemoveRotation);

        // [Lanczos] //

        lanczosTolerance = ini.GetValueF("Lanczos", "tolerance", lanczosTolerance);
        lanczosMaxIterations = ini.GetValueL("Lanczos", "max_iterations", lanczosMaxIterations);
        lanczosQuitEarly = ini.GetValueB("Lanczos", "quit_early", lanczosQuitEarly);

        // [GPR Dimer] //
        gprDimerRotationAngle = ini.GetValueF("GPR Dimer", "finite_angle", gprDimerRotationAngle);
        gprDimerConvergedAngle
            = ini.GetValueF("GPR Dimer", "converged_angle", gprDimerConvergedAngle);
        gprDimerRelaxConvAngle
            = ini.GetValueF("GPR Dimer", "relaxation_converged_angle", gprDimerRelaxConvAngle);
        gprDimerInitRotationsMax = (int) ini.GetValueL(
            "GPR Dimer", "max_initial_rotation_iterations", gprDimerInitRotationsMax);
        gprDimerRelaxRotationsMax = (int) ini.GetValueL(
            "GPR Dimer", "max_relaxation_rotation_iterations", gprDimerRelaxRotationsMax);
        gprDimerDivisorTdimerGP
            = (int) ini.GetValueL("GPR Dimer", "divisor_t_dimer", gprDimerDivisorTdimerGP);
        gprDimerMaxOuterIterations
            = (int) ini.GetValueL("GPR Dimer", "max_outer_iterations", gprDimerMaxOuterIterations);
        gprDimerMaxInnerIterations
            = (int) ini.GetValueL("GPR Dimer", "max_inner_iterations", gprDimerMaxInnerIterations);
        gprDimerMidpointMaxDisp
            = ini.GetValueF("GPR Dimer", "max_midpoint_displacement", gprDimerMidpointMaxDisp);
        gprDimerRotOptMethod
            = ini.GetValue("GPR Dimer", "rotation_opt_method", gprDimerRotOptMethod);
        gprDimerTransOptMethod
            = ini.GetValue("GPR Dimer", "translation_opt_method", gprDimerTransOptMethod);
        gprActiveRadius = ini.GetValueF("GPR Dimer", "active_radius", gprActiveRadius);
        gprDimerSep = ini.GetValueF("GPR Dimer", "dimer_separation", gprDimerSep);
        gprDimerConvStep = ini.GetValueF("GPR Dimer", "convex_region_step_size", gprDimerConvStep);
        gprDimerMaxStep = ini.GetValueF("GPR Dimer", "max_step_size", gprDimerMaxStep);
        gprForceThreshold = ini.GetValueF("GPR Dimer", "force_threshold", gprForceThreshold);
        gprDimerRatioAtLimit = ini.GetValueF("GPR Dimer", "ratio_at_limit", gprDimerRatioAtLimit);
        gprDimerInitRotGP
            = ini.GetValueB("GPR Dimer", "nogp_initial_rotations", gprDimerInitRotGP);
        gprDimerInitTransGP
            = ini.GetValueB("GPR Dimer", "nogp_init_translations", gprDimerInitTransGP);
        gprDimerManyIterations
            = ini.GetValueB("GPR Dimer", "has_many_iterations", gprDimerManyIterations);
        // GPR Params
        gprDimerHyperOptMethod
            = ini.GetValue("GPR Dimer", "hyperparameter_opt_method", gprDimerHyperOptMethod);
        gprDimerSigma2 = ini.GetValueF("GPR Dimer", "gpr_variance", gprDimerSigma2);
        gprDimerJitterSigma2
            = ini.GetValueF("GPR Dimer", "gpr_jitter_variance", gprDimerJitterSigma2);
        gprDimerNoiseSigma2
            = ini.GetValueF("GPR Dimer", "gpr_noise_variance", gprDimerNoiseSigma2);
        gprDimerPriorMu = ini.GetValueF("GPR Dimer", "prior_mean", gprDimerPriorMu);
        gprDimerPriorSigma2 = ini.GetValueF("GPR Dimer", "prior_variance", gprDimerPriorSigma2);
        gprDimerPriorNu = ini.GetValueF("GPR Dimer", "prior_degrees_of_freedom", gprDimerPriorNu);
        // GPR Optimization Parameters
        gprOptCheckDerivatives
            = ini.GetValueB("GPR Dimer", "check_derivatives", gprOptCheckDerivatives);
        gprOptMaxIterations
            = (int) ini.GetValueL("GPR Dimer", "opt_max_iterations", gprOptMaxIterations);
        gprOptTolFunc = ini.GetValueF("GPR Dimer", "opt_tol_func", gprOptTolFunc);
        gprOptTolSol = ini.GetValueF("GPR Dimer", "opt_tol_sol", gprOptTolSol);
        gprOptLambdaLimit = ini.GetValueF("GPR Dimer", "opt_lambda_limit", gprOptLambdaLimit);
        gprOptLambdaInit = ini.GetValueF("GPR Dimer", "opt_lambda_init", gprOptLambdaInit);
        // GPR Debugging Parameters
        gprReportLevel = (int) ini.GetValueL("GPR Dimer", "report_level", gprReportLevel);
        gprDebugLevel = (int) ini.GetValueL("GPR Dimer", "debug_level", gprDebugLevel);
        gprDebugOutDir = ini.GetValue("GPR Dimer", "debug_output_directory", gprDebugOutDir);
        ;
        gprDebugPosFile = ini.GetValue("GPR Dimer", "debug_position_basename", gprDebugPosFile);
        ;
        gprDebugEnergyFile
            = ini.GetValue("GPR Dimer", "debug_energy_basename", gprDebugEnergyFile);
        ;
        gprDebugGradFile = ini.GetValue("GPR Dimer", "debug_gradient_basename", gprDebugGradFile);
        ;
        gprDebugOffsetMidPoint
            = ini.GetValueF("GPR Dimer", "debug_midpoint_offset", gprDebugOffsetMidPoint);
        gprDebugDy = ini.GetValueF("GPR Dimer", "debug_y_step", gprDebugDy);
        gprDebugDz = ini.GetValueF("GPR Dimer", "debug_z_step", gprDebugDz);
        // GPR Prune
        gprUsePrune = ini.GetValueB("GPR Dimer", "use_prune", gprUsePrune);
        gprPruneBegin = (int) ini.GetValueL("GPR Dimer", "start_prune_at", gprPruneBegin);
        gprPruneNVals = (int) ini.GetValueL("GPR Dimer", "nprune_vals", gprPruneNVals);
        gprPruneThreshold = ini.GetValueF("GPR Dimer", "prune_threshold", gprPruneThreshold);

        // [Prefactor] //

        prefactorDefaultValue = ini.GetValueF("Prefactor", "default_value", prefactorDefaultValue);
        prefactorMaxValue = ini.GetValueF("Prefactor", "max_value", prefactorMaxValue);
        prefactorMinValue = ini.GetValueF("Prefactor", "min_value", prefactorMinValue);
        prefactorWithinRadius = ini.GetValueF("Prefactor", "within_radius", prefactorWithinRadius);
        prefactorMinDisplacement
            = ini.GetValueF("Prefactor", "min_displacement", prefactorMinDisplacement);
        prefactorRate = toLowerCase(ini.GetValue("Prefactor", "rate_estimation", prefactorRate));
        prefactorConfiguration
            = toLowerCase(ini.GetValue("Prefactor", "configuration", prefactorConfiguration));
        prefactorAllFreeAtoms
            = ini.GetValueB("Prefactor", "all_free_atoms", prefactorAllFreeAtoms);
        prefactorFilterScheme
            = toLowerCase(ini.GetValue("Prefactor", "filter_scheme", prefactorFilterScheme));
        prefactorFilterFraction
            = ini.GetValueF("Prefactor", "filter_fraction", prefactorFilterFraction);

        // [Hessian] //

        hessianAtomList = toLowerCase(ini.GetValue("Hessian", "atom_list", hessianAtomList));
        hessianZeroFreqValue = ini.GetValueF("Hessian", "zero_freq_value", hessianZeroFreqValue);

        // [Nudged Elastic Band] //

        nebImages = ini.GetValueL("Nudged Elastic Band", "images", nebImages);
        nebSpring = ini.GetValueF("Nudged Elastic Band", "spring", nebSpring);
        nebClimbingImageMethod = ini.GetValueB(
            "Nudged Elastic Band", "climbing_image_method", nebClimbingImageMethod);
        nebClimbingImageConvergedOnly = ini.GetValueB(
            "Nudged Elastic Band", "climbing_image_converged_only", nebClimbingImageConvergedOnly);
        nebOldTangent = ini.GetValueB("Nudged Elastic Band", "old_tangent", nebOldTangent);
        nebMaxIterations
            = ini.GetValueL("Nudged Elastic Band", "max_iterations", optMaxIterations);
        nebDoublyNudged = ini.GetValueB("Nudged Elastic Band", "doubly_nudged", nebDoublyNudged);
        nebDoublyNudgedSwitching = ini.GetValueB(
            "Nudged Elastic Band", "doubly_nudged_switching", nebDoublyNudgedSwitching);
        nebElasticBand = ini.GetValueB("Nudged Elastic Band", "elastic_band", nebElasticBand);
        nebConvergedForce
            = ini.GetValueF("Nudged Elastic Band", "converged_force", optConvergedForce);

        // [Dynamics] //

        mdTimeStepInput = ini.GetValueF("Dynamics", "time_step", mdTimeStepInput);
        mdTimeStep = mdTimeStepInput / timeUnit;
        mdTimeInput = ini.GetValueF("Dynamics", "time", mdTimeInput);
        mdTime = mdTimeInput / timeUnit;
        mdSteps = long(floor(mdTime / mdTimeStep + 0.5));
        thermostat = toLowerCase(ini.GetValue("Dynamics", "thermostat", "andersen"));
        thermoAndersenAlpha = ini.GetValueF("Dynamics", "andersen_alpha", thermoAndersenAlpha);
        thermoAndersenTcolInput
            = ini.GetValueF("Dynamics", "andersen_collision_period", thermoAndersenTcolInput);
        thermoAndersenTcol = thermoAndersenTcolInput / timeUnit;
        thermoNoseMass = ini.GetValueF("Dynamics", "nose_mass", thermoNoseMass);
        thermoLangevinFrictionInput
            = ini.GetValueF("Dynamics", "langevin_friction", thermoLangevinFrictionInput);
        thermoLangevinFriction = thermoLangevinFrictionInput * timeUnit;
        // thermoAtoms = helper_functions::split_std::string_int(ini.GetValue("Dynamics",
        // "thermo_atoms", ""), ",");

        // [Parallel Replica]

        parrepAutoStop
            = ini.GetValueB("Parallel Replica", "stop_after_transition", parrepAutoStop);
        parrepRefineTransition
            = ini.GetValueB("Parallel Replica", "refine_transition", parrepRefineTransition);
        parrepDephaseLoopStop
            = ini.GetValueB("Parallel Replica", "dephase_loop_stop", parrepDephaseLoopStop);
        parrepDephaseTimeInput
            = ini.GetValueF("Parallel Replica", "dephase_time", parrepDephaseTimeInput);
        parrepDephaseTime = parrepDephaseTimeInput / timeUnit;
        parrepDephaseLoopMax
            = ini.GetValueL("Parallel Replica", "dephase_loop_max", parrepDephaseLoopMax);
        parrepStateCheckIntervalInput = ini.GetValueF(
            "Parallel Replica", "state_check_interval", parrepStateCheckIntervalInput);
        parrepStateCheckInterval = parrepStateCheckIntervalInput / timeUnit;
        parrepRecordIntervalInput = ini.GetValueF(
            "Parallel Replica", "state_save_interval", 0.1 * parrepStateCheckIntervalInput);
        parrepRecordInterval = parrepRecordIntervalInput / timeUnit;
        parrepCorrTimeInput
            = ini.GetValueF("Parallel Replica", "post_transition_time", parrepCorrTimeInput);
        parrepCorrTime = parrepCorrTimeInput / timeUnit;

        //[Temperature Accelerated Dynamics] //

        tadLowT = ini.GetValueF("TAD", "low_temperature", tadLowT);
        tadMinPrefactor = ini.GetValueF("TAD", "min_prefactor", tadMinPrefactor);
        tadConfidence = ini.GetValueF("TAD", "confidence", tadConfidence);

        // [Replica Exchange] //

        repexcTemperatureDistribution = toLowerCase(ini.GetValue(
            "Replica Exchange", "temperature_distribution", repexcTemperatureDistribution));
        repexcReplicas = ini.GetValueL("Replica Exchange", "replicas", repexcReplicas);
        repexcExchangeTrials
            = ini.GetValueL("Replica Exchange", "exchange_trials", repexcExchangeTrials);
        repexcSamplingTimeInput
            = ini.GetValueF("Replica Exchange", "sampling_time", repexcSamplingTimeInput);
        repexcSamplingTime = repexcSamplingTimeInput / timeUnit;
        repexcTemperatureLow = ini.GetValueF("Replica Exchange", "temperature_low", temperature);
        repexcTemperatureHigh
            = ini.GetValueF("Replica Exchange", "temperature_high", repexcTemperatureHigh);
        repexcExchangePeriodInput
            = ini.GetValueF("Replica Exchange", "exchange_period", repexcExchangePeriodInput);
        repexcExchangePeriod = repexcExchangePeriodInput / timeUnit;

        // [Hyperdynamics] //

        bondBoostRMDTimeInput
            = ini.GetValueF("Hyperdynamics", "bb_rmd_time", bondBoostRMDTimeInput);
        bondBoostRMDTime = bondBoostRMDTimeInput / timeUnit;
        bondBoostBALS
            = toLowerCase(ini.GetValue("Hyperdynamics", "bb_boost_atomlist", bondBoostBALS));
        bondBoostDVMAX = ini.GetValueF("Hyperdynamics", "bb_dvmax", bondBoostDVMAX);
        bondBoostQRR = ini.GetValueF("Hyperdynamics", "bb_stretch_threshold", bondBoostQRR);
        bondBoostPRR = ini.GetValueF("Hyperdynamics", "bb_ds_curvature", bondBoostPRR);
        bondBoostQcut = ini.GetValueF("Hyperdynamics", "bb_rcut", bondBoostQcut);
        biasPotential
            = toLowerCase(ini.GetValue("Hyperdynamics", "bias_potential", biasPotential));

        // [Saddle Search] //

        saddleMethod = toLowerCase(ini.GetValue("Saddle Search", "method", saddleMethod));
        saddleMinmodeMethod
            = toLowerCase(ini.GetValue("Saddle Search", "min_mode_method", saddleMinmodeMethod));
        saddleDisplaceMagnitude
            = ini.GetValueF("Saddle Search", "displace_magnitude", saddleDisplaceMagnitude);
        saddleDisplaceRadius
            = ini.GetValueF("Saddle Search", "displace_radius", saddleDisplaceRadius);
        saddleMaxEnergy = ini.GetValueF("Saddle Search", "max_energy", saddleMaxEnergy);
        saddleMaxIterations = ini.GetValueL("Saddle Search", "max_iterations", optMaxIterations);
        saddleNonnegativeDisplacementAbort = ini.GetValueB(
            "Saddle Search", "nonnegative_displacement_abort", saddleNonnegativeDisplacementAbort);
        saddleMaxSingleDisplace
            = ini.GetValueF("Saddle Search", "max_single_displace", saddleMaxSingleDisplace);
        // must be loaded after optConvergedForce
        saddleConvergedForce
            = ini.GetValueF("Saddle Search", "converged_force", optConvergedForce);
        saddlePerpForceRatio = ini.GetValueF(
            "Saddle Search", "perp_force_ratio", saddlePerpForceRatio); // undocumented
        saddleDisplaceType = toLowerCase(
            ini.GetValue("Saddle Search", "client_displace_type", EpiCentersStrings::DISP_LOAD));
        saddleNonlocalCountAbort = ini.GetValueL(
            "Saddle Search", "nonlocal_count_abort", saddleNonlocalCountAbort); // undocumented
        saddleNonlocalDistanceAbort = ini.GetValueF("Saddle Search",
                                                    "nonlocal_distance_abort",
                                                    saddleNonlocalDistanceAbort); // undocumented
        if (saddleDisplaceType != EpiCentersStrings::DISP_NOT_FCC_OR_HCP
            && saddleDisplaceType != EpiCentersStrings::DISP_MIN_COORDINATED
            && saddleDisplaceType != EpiCentersStrings::DISP_LAST_ATOM
            && saddleDisplaceType != EpiCentersStrings::DISP_RANDOM) {
            saddleDisplaceType = EpiCentersStrings::DISP_LOAD;
        }
        saddleConfinePositive
            = ini.GetValueB("Saddle Search", "confine_positive", saddleConfinePositive);
        if (saddleConfinePositive) {
            saddleBowlBreakout
                = ini.GetValueB("Saddle Search", "bowl_breakout", saddleBowlBreakout);
            saddleBowlActive
                = ini.GetValueL("Saddle Search", "bowl_active_atoms", saddleBowlActive);
            saddleConfinePositiveMinForce = ini.GetValueF(
                "Saddle Search", "confine_positive_min_move", saddleConfinePositiveMinForce);
            saddleConfinePositiveScaleRatio = ini.GetValueF(
                "Saddle Search", "confine_positive_scale_ratio", saddleConfinePositiveScaleRatio);
            saddleConfinePositiveBoost = ini.GetValueF(
                "Saddle Search", "confine_positive_boost", saddleConfinePositiveBoost);
            saddleConfinePositiveMinActive = ini.GetValueL(
                "Saddle Search", "confine_positive_min_active", saddleConfinePositiveMinActive);
        }
        saddleDynamicsTemperature = temperature;
        saddleDynamicsTemperature
            = ini.GetValueF("Saddle Search", "dynamics_temperature", saddleDynamicsTemperature);
        saddleDynamicsStateCheckIntervalInput
            = ini.GetValueF("Saddle Search",
                            "dynamics_state_check_interval",
                            saddleDynamicsStateCheckIntervalInput);
        saddleDynamicsStateCheckInterval = saddleDynamicsStateCheckIntervalInput / timeUnit;
        saddleDynamicsRecordIntervalInput = ini.GetValueF(
            "Saddle Search", "dynamics_record_interval", saddleDynamicsRecordIntervalInput);
        saddleDynamicsRecordInterval = saddleDynamicsRecordIntervalInput / timeUnit;
        saddleDynamicsLinearInterpolation = ini.GetValueB(
            "Saddle Search", "dynamics_linear_interpolation", saddleDynamicsLinearInterpolation);
        saddleRemoveRotation
            = ini.GetValueB("Saddle Search", "remove_rotation", saddleRemoveRotation);
        saddleDynamicsMaxInitCurvature = ini.GetValueF(
            "Saddle Search", "dynamics_max_init_curvature", saddleDynamicsMaxInitCurvature);
        saddleZeroModeAbortCurvature = ini.GetValueF(
            "Saddle Search", "zero_mode_abort_curvature", saddleZeroModeAbortCurvature);

        // [Basin Hopping] //

        basinHoppingDisplacement
            = ini.GetValueF("Basin Hopping", "displacement", basinHoppingDisplacement);
        basinHoppingPushApartDistance
            = ini.GetValueF("Basin Hopping", "push_apart_distance", basinHoppingPushApartDistance);
        basinHoppingInitialRandomStructureProbability
            = ini.GetValueF("Basin Hopping",
                            "initial_random_structure_probability",
                            basinHoppingInitialRandomStructureProbability);
        basinHoppingSteps = ini.GetValueL("Basin Hopping", "steps", basinHoppingSteps);
        basinHoppingQuenchingSteps
            = ini.GetValueL("Basin Hopping", "quenching_steps", basinHoppingQuenchingSteps);
        basinHoppingSingleAtomDisplace = ini.GetValueB(
            "Basin Hopping", "single_atom_displace", basinHoppingSingleAtomDisplace);
        basinHoppingSignificantStructure = ini.GetValueB(
            "Basin Hopping", "significant_structure", basinHoppingSignificantStructure);
        basinHoppingDisplacementAlgorithm = toLowerCase(ini.GetValue(
            "Basin Hopping", "displacement_algorithm", basinHoppingDisplacementAlgorithm));
        basinHoppingDisplacementDistribution = toLowerCase(ini.GetValue(
            "Basin Hopping", "displacement_distribution", basinHoppingDisplacementDistribution));
        basinHoppingSwapProbability
            = ini.GetValueF("Basin Hopping", "swap_probability", basinHoppingSwapProbability);
        basinHoppingJumpMax = ini.GetValueL("Basin Hopping", "jump_max", basinHoppingJumpMax);
        basinHoppingJumpSteps
            = ini.GetValueL("Basin Hopping", "jump_steps", basinHoppingJumpSteps);
        basinHoppingAdjustDisplacement = ini.GetValueB(
            "Basin Hopping", "adjust_displacement", basinHoppingAdjustDisplacement);
        basinHoppingAdjustPeriod
            = ini.GetValueL("Basin Hopping", "adjust_period", basinHoppingAdjustPeriod);
        basinHoppingAdjustFraction
            = ini.GetValueF("Basin Hopping", "adjust_fraction", basinHoppingAdjustFraction);
        basinHoppingTargetRatio
            = ini.GetValueF("Basin Hopping", "target_ratio", basinHoppingTargetRatio);
        basinHoppingWriteUnique
            = ini.GetValueB("Basin Hopping", "write_unique", basinHoppingWriteUnique);
        basinHoppingStopEnergy
            = ini.GetValueF("Basin Hopping", "stop_energy", basinHoppingStopEnergy);

        // [Global Optimization] //

        globalOptimizationMoveMethod = toLowerCase(
            ini.GetValue("Global Optimization", "move_method", globalOptimizationMoveMethod));
        globalOptimizationDecisionMethod = toLowerCase(ini.GetValue(
            "Global Optimization", "decision_method", globalOptimizationDecisionMethod));
        globalOptimizationSteps
            = ini.GetValueL("Global Optimization", "steps", globalOptimizationSteps);
        globalOptimizationBeta
            = ini.GetValueF("Global Optimization", "beta", globalOptimizationBeta);
        globalOptimizationAlpha
            = ini.GetValueF("Global Optimization", "alpha", globalOptimizationAlpha);
        globalOptimizationMdmin
            = ini.GetValueL("Global Optimization", "mdmin", globalOptimizationMdmin);
        globalOptimizationTargetEnergy = ini.GetValueF(
            "Global Optimization", "target_energy", globalOptimizationTargetEnergy);

        // [BGSD] //

        alpha = ini.GetValueF("BGSD", "alpha", alpha);
        beta = ini.GetValueF("BGSD", "beta", beta);
        gradientfinitedifference
            = ini.GetValueF("BGSD", "gradientfinitedifference", gradientfinitedifference);
        grad2energyconvergence
            = ini.GetValueF("BGSD", "grad2energyconvergence", grad2energyconvergence);
        grad2forceconvergence
            = ini.GetValueF("BGSD", "grad2forceconvergence", grad2forceconvergence);

        // [Monte Carlo] //

        monteCarloStepSize = ini.GetValueF("Monte Carlo", "step_size", monteCarloStepSize);
        monteCarloSteps = ini.GetValueL("Monte Carlo", "steps", monteCarloSteps);

        log_close();
        log_init(this, (char *) "client.log");

        // Sanity Checks
        if (parrepStateCheckInterval > mdTime && job == "parallel_replica") {
            char msg[] = "error: [Parallel Replica] state_check_interval must be <= time\n";
            fprintf(stderr, "%s", msg);
            log(msg);
            exit(1);
        }

        if (saddleDynamicsRecordIntervalInput > saddleDynamicsStateCheckIntervalInput) {
            char msg[] = "error:  [Saddle Search] dynamics_record_interval must be <= "
                         "dynamics_state_check_interval\n";
            fprintf(stderr, "%s", msg);
            log(msg);
            exit(1);
        }

        if (potential == "ams"s || potential == "ams_io"s) {
            if (forcefield.empty() && model.empty() && xc.empty()) {
                char msg[] = "error:  [AMS] Must provide atleast forcefield or model or xc\n";
                fprintf(stderr, "%s", msg);
                log(msg);
                exit(1);
            }

            if (!forcefield.empty() && !model.empty() && !xc.empty()) {
                char msg[] = "error:  [AMS] Must provide either forcefield or model\n";
                fprintf(stderr, "%s", msg);
                log(msg);
                exit(1);
            }
        }

    } else {
        fprintf(stderr, "Couldn't parse the ini file.\n");
        error = 1;
    }

    return error;
}

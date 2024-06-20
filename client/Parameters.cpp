#include "Parameters.h"
#include "BaseStructures.h"
#include "EpiCenters.h"
#include "ImprovedDimer.h"
#include "ParallelReplicaJob.h"
#include "Prefactor.h"
#include "PrefactorJob.h"
#include "ReplicaExchangeJob.h"
#include "magic_enum/magic_enum.hpp"
#include <errno.h>
#include <float.h>
#include <string>
#include <thirdparty/toml.hpp>

using namespace std::string_literals;

Parameters::Parameters() {

  kB = 8.6173324e-5;     // eV/K
  timeUnit = 10.1805055; // fs
  MPIPotentialRank = -1;

  // [Main] //
  main.job = JobType::Process_Search;
  main.randomSeed = -1;
  main.temperature = 300.0;
  main.checkpoint = false;
  main.quiet = false;
  main.writeLog = true;
  main.inpFilename = "config.toml"s;
  main.conFilename = "pos.con"s;
  main.finiteDifference = 0.01;
  main.maxForceCalls = 0;
  main.removeNetForce = true;

  // [Prefactor] //
  prefactor.defaultValue = 0.0;
  prefactor.maxValue = 1e+21;
  prefactor.minValue = 1e+9;
  prefactor.withinRadius = 3.3;
  prefactor.minDisplacement = 0.25;
  prefactor.rate = ::Prefactor::RATE_HTST;
  prefactor.configuration = PrefactorType::REACTANT;
  prefactor.allFreeAtoms = false;
  prefactor.filterScheme = ::Prefactor::FILTER_FRACTION;
  prefactor.filterFraction = 0.90;

  // [Potential] //
  pot.potential = PotType::LJ;
  pot.MPIPollPeriod = 0.25; // seconds
  pot.LogPotential = false;
  pot.LAMMPSLogging = false;
  pot.LAMMPSThreads = 0;
  pot.EMTRasmussen = false;
  pot.extPotPath = "./ext_pot"s;

  // [AMS] //
  ams.engine = ""s;     // One of REAXFF MOPAC
  ams.forcefield = ""s; // OPt.ff or something else
  ams.model = ""s;      // PM7 PM3 or something
  ams.xc = ""s;         // exchange-correlation functional
  ams.basis = ""s;      // with xc
  ams.resources = ""s;  // For DFTB

  // [AMS_ENV] //
  // Horrid little section to mimic amsrc.sh
  // Assumes the entire thing is going to be set
  amsenv.amshome = ""s;       // "/some/path/to/amshome/"s;
  amsenv.scm_tmpdir = ""s;    // "/tmp"s;
  amsenv.scm_pythondir = ""s; // "/.scm/python"s;
  amsenv.amsbin = ""s;        // amshome.append("/bin"s);
  amsenv.scmlicense = ""s;    // amshome.append("license.txt"s);
  amsenv.amsresources = ""s;  // amshome.append("/atomicdata"s);

  // [XTBPot] //
  xtbpot.paramset = "GFNFF"s;
  xtbpot.acc = 1.0;
  xtbpot.elec_temperature = 0.0;
  xtbpot.maxiter = 250;

  // [Structure Comparison] //
  structcomp.distanceDifference = 0.1;
  structcomp.neighborCutoff = 3.3;
  structcomp.checkRotation = false;
  structcomp.indistinguishableAtoms = true;
  structcomp.energyDifference = 0.01;
  structcomp.removeTranslation = true;

  // [Debug] //
  debug.writeMovies = false;
  debug.writeMoviesInterval = 1;

  // [Saddle Search] //
  saddle.displaceType = EpiCenters::DISP_LOAD;
  saddle.method = "min_mode"s;
  saddle.minmodeMethod = LowestEigenmode::MINMODE_DIMER;
  saddle.maxEnergy = 20.0;
  saddle.maxIterations = 1000;
  saddle.displaceRadius = 4.0;
  saddle.displaceMagnitude = 0.1;
  saddle.maxSingleDisplace = 10.;
  //    saddle.ConvergedForce = optConvergedForce; default value is set after
  //    the value of optConvergedForce is loaded
  saddle.nonnegativeDisplacementAbort = false;
  saddle.nonlocalCountAbort = 0;
  saddle.nonlocalDistanceAbort = 0.0;
  saddle.removeRotation = false;
  saddle.perpForceRatio = 0.0;    // undocumented
  saddle.confinePositive = false; // undocumented
  saddle.bowlBreakout = false;    // undocumented
  saddle.bowlActive = 20;
  saddle.confinePositiveMinForce = 0.5;           // undocumented
  saddle.confinePositiveScaleRatio = 0.9;         // undocumented
  saddle.confinePositiveBoost = 10.;              // undocumented
  saddle.confinePositiveMinActive = 30;           // undocumented
  saddle.dynamicsTemperature = 0.0;               // defaults to temperature
  saddle.dynamicsStateCheckIntervalInput = 100.0; // fs
  saddle.dynamicsStateCheckInterval =
      saddle.dynamicsStateCheckIntervalInput / timeUnit;
  saddle.dynamicsRecordIntervalInput = 10.0;      // fs
  saddle.dynamicsRecordInterval =
    saddle.dynamicsRecordIntervalInput / timeUnit;
  saddle.dynamicsLinearInterpolation = true;
  saddle.dynamicsMaxInitCurvature = 0.0; // eV/Ang^2
  saddle.zeroModeAbortCurvature = 0.0;   // eV/Ang^2

  // [Optimizers] //
  optim.method = OptType::CG;
  optim.convergenceMetric = "norm"s;
  optim.refineOptMethod = OptType::None;
  optim.refineThreshold = 0.5;
  optim.maxIterations = 1000;
  optim.convergedForce = 0.01;
  optim.maxMove = 0.2;
  optim.timeStepInput = 1.0;
  optim.maxTimeStepInput = 2.5;
  optim.maxTimeStep = optim.maxTimeStepInput / timeUnit;

  optim.LBFGSMemory = 20;
  // assumes stiffest curvature at minimum is 100 eV/A^2
  optim.LBFGSInverseCurvature = 0.01;
  optim.LBFGSAutoScale = true;
  optim.LBFGSAngleReset = true;
  optim.LBFGSDistanceReset = true;

  optim.QMSteepestDecent = false;
  optim.CGNoOvershooting = false;
  optim.CGKnockOutMaxMove = false;
  optim.CGLineConverged = 0.1;
  optim.CGLineSearch = false;
  optim.CGMaxIterBeforeReset = 0;
  optim.CGLineSearchMaxIter = 10;
  optim.SDAlpha = 0.1;
  optim.SDTwoPoint = false;

  // [Process Search] //
  procsearch.minimizeFirst = true;
  procsearch.minimizationOffset = optim.maxMove;

  // [Dimer] //
  dimer.rotationAngle = 0.005;
  dimer.improved = true;
  dimer.convergedAngle = 5.0; // degrees
  dimer.maxIterations = 1000;
  dimer.optMethod = ImprovedDimer::OPT_CG;
  dimer.torqueMin = 0.1;   // old dimer
  dimer.torqueMax = 1.0;   // old dimer
  dimer.rotationsMin = 1;  // old dimer
  dimer.rotationsMax = 10; // old dimer and new dimer
  dimer.removeRotation = false;

  // [ASE_ORCA] //
  aseorca.orca_path = ""s;
  aseorca.orca_nproc = "1"s;

  // [Lanczos] //
  lanczos.tolerance = 0.01;
  lanczos.maxIterations = 20;
  lanczos.quitEarly = true;

  // [GPR Dimer] //
  gprd.rotationAngle = 0.005;
  gprd.convergedAngle = 0.08;     // T_anglerot_init
  gprd.relaxConvAngle = 0.001;    // T_anglerot_gp
  gprd.initRotationsMax = 6;      // num_iter_initrot; should be DoF
  gprd.relaxRotationsMax = 10;    // num_iter_rot_gp
  gprd.divisorTdimerGP = 10;      // divisor_T_dimer_gp
  gprd.maxOuterIterations = 300;  // num_bigiter
  gprd.maxInnerIterations = 1000; // num_iter
  gprd.midpointMaxDisp = 0.5;     // disp_max
  gprd.rotOptMethod = "lbfgs"s;   // method_rot
  gprd.transOptMethod = "lbfgs"s; // method_trans
  gprd.activeRadius = 5.0;        // actidst_fro
  gprd.dimerSep = 0.01;           // dimer_sep
  gprd.convStep = 0.1;            // param_trans[0]
  gprd.maxStep = 0.1;             // param_trans[1]
  gprd.forceThreshold = 0.01;     // T_dimer
  gprd.ratioAtLimit = 0.66667;    // ratio_at_limit
  gprd.initRotGP = 0;             // initrot_nogp
  gprd.initTransGP = 0;           // inittrans_nogp
  gprd.manyIterations = true;     // islarge_num_iter
  // GPR Params
  gprd.hyperOptMethod = "scg"s; // optimization_alg
  gprd.sigma2 = 1e-8;           // gp_sigma2
  gprd.jitterSigma2 = 0;        // jitter_sigma2
  gprd.noiseSigma2 = 1e-8;      // sigma2
  gprd.priorMu = 0;             // prior_mu
  gprd.priorSigma2 = 1;         // prior_s2
  gprd.priorNu = 20;            // prior_nu
  // GPR Optimization Parameters
  gprd.optCheckDerivatives = false; // check_derivative
  gprd.optMaxIterations = 400;      // max_iter
  gprd.optTolFunc = 1e-4;           // tolerance_func
  gprd.optTolSol = 1e-4;            // tolerance_sol
  gprd.optLambdaLimit = 1e17;       // lambda_limit
  gprd.optLambdaInit = 10.0;        // lambda
  // GPR Prune Parameters
  gprd.usePrune = false;     // use_prune
  gprd.pruneBegin = 8;       // start_prune_at
  gprd.pruneNVals = 3;       // nprune_vals
  gprd.pruneThreshold = 0.5; // prune_threshold
  // GPR Debugging Parameters
  gprd.reportLevel = 1;             // report_level
  gprd.debugLevel = 2;              // debug_level
  gprd.debugOutDir = "output"s;     // debug_output_dir
  gprd.debugPosFile = "position"s;  // debug_output_file_R
  gprd.debugEnergyFile = "energy"s; // debug_output_file_E
  gprd.debugGradFile = "gradient"s; // debug_output_file_G
  gprd.debugOutExt = "dat"s;        // debug_output_file_extension
  gprd.debugOffsetMidPoint = 3.;    // debug_offset_from_mid_point
  gprd.debugDy = 0.1;               // debug_dy
  gprd.debugDz = 0.1;               // debug_dz

  // GP Surrogate Parameters
  surrogate.use = false;
  surrogate.sub_job = JobType::Unknown;
  surrogate.gp_uncertainity = 0.05;
  surrogate.gp_linear_path_always = false;
  surrogate.potential = PotType::CatLearn;

  // [Hessian] //
  hessian.atomList = std::string("All"s);
  hessian.zeroFreqValue = 1e-6;

  // [Nudged Elastic Band] //
  neb.images = 5;
  neb.spring = 5.0;
  neb.climbingImageMethod = true;
  neb.climbingImageConvergedOnly = true;
  neb.oldTangent = false;
  neb.maxIterations = 1000;
  neb.doublyNudged = false;
  neb.doublyNudgedSwitching = false;
  neb.elasticBand = false;
  neb.convergedForce = optim.convergedForce;
  neb.energyWeighted = false;
  neb.KSPMin = 0.97;
  neb.KSPMax = 9.7;

  // [Dynamics] //
  md.timeStepInput = 1.0;
  md.timeInput = 1000.0;
  md.timeStep = md.timeStepInput / timeUnit;
  md.time = md.timeInput / timeUnit;
  md.steps = long(floor(md.time / md.timeStep + 0.5));

  // [Thermostat] //
  thermostat.kind = Dynamics::NONE;
  thermostat.andersenAlpha = 1.0;       // collision strength
  thermostat.andersenTcolInput = 100.0; // collision frequency in unit of fs
  thermostat.noseMass = 1.0;
  thermostat.langevinFrictionInput = 0.01;
  thermostat.langevinFriction = thermostat.langevinFrictionInput * timeUnit;

  // [Parallel Replica] //
  parrep.refineTransition = true;
  parrep.autoStop = false;
  parrep.dephaseLoopStop = false;
  parrep.dephaseTimeInput = 1000.0;
  parrep.dephaseLoopMax = 5;
  parrep.stateCheckIntervalInput = 1000.0;
  parrep.recordIntervalInput = 50.0;
  parrep.corrTimeInput = 1000.0;

  parrep.dephaseTime = parrep.dephaseTimeInput / timeUnit;
  parrep.stateCheckInterval = parrep.stateCheckIntervalInput / timeUnit;
  parrep.recordInterval = parrep.recordIntervalInput / timeUnit;
  parrep.corrTime = parrep.corrTimeInput / timeUnit;

  // [Temperature Accelerated Dynamics] //
  tad.lowT = 300.0;
  tad.minPrefactor = 0.001; // in unit of fs-1
  tad.confidence = 0.001;

  // [Replica Exchange] //
  repexc.temperatureDistribution = "exponential"s;
  repexc.replicas = 10;
  repexc.exchangeTrials = repexc.replicas;
  repexc.samplingTimeInput = 1000.0;
  repexc.samplingTime = repexc.samplingTimeInput / timeUnit;
  repexc.temperatureLow = 0.0;
  repexc.temperatureHigh = 0.0;
  repexc.exchangePeriod = 100.0;

  // [Hyperdynamics] //
  bondBoost.biasPotential = "none"s;    // Hyperdynamics::NONE
  bondBoost.BALS = "ALL"s; // boosted atom list string
  bondBoost.DVMAX = 0.0;
  bondBoost.QRR = 0.2; // can not be set to 0
  bondBoost.PRR = 0.95;
  bondBoost.Qcut = 3.0;
  bondBoost.RMDTimeInput = 100.0;
  bondBoost.RMDTime = bondBoost.RMDTimeInput / timeUnit;

  // [Basin Hopping] //
  bhop.displacement = 0.5;
  bhop.pushApartDistance = 0.4;
  bhop.initialRandomStructureProbability = 0.0;
  bhop.steps = 10000;
  bhop.quenchingSteps = 0;
  bhop.singleAtomDisplace = false;
  bhop.significantStructure = true;
  bhop.displacementAlgorithm = "standard"s;
  bhop.displacementDistribution = "uniform"s;
  bhop.swapProbability = 0.0;
  bhop.jumpMax = 10;
  bhop.jumpSteps = 0;
  bhop.adjustDisplacement = true;
  bhop.adjustPeriod = 10;
  bhop.adjustFraction = 0.05;
  bhop.targetRatio = 0.5;
  bhop.writeUnique = false;
  bhop.stopEnergy = -DBL_MAX;

  // [Global Optimization] //
  globopt.moveMethod = "md"s;
  globopt.decisionMethod = "npew"s;
  globopt.steps = 10000;
  globopt.beta = 1.05;
  globopt.alpha = 1.02;
  globopt.mdmin = 3;
  globopt.targetEnergy = -1.E50;

  // [Monte Carlo] //
  monte_carlo.stepSize = 0.005;
  monte_carlo.steps = 1000;

  // [BGSD] //
  bgsd.alpha = 10.0;
  bgsd.beta = 0.2;
  bgsd.gradientfinitedifference = 0.000001;
  bgsd.Hforceconvergence = 0.01;
  bgsd.grad2energyconvergence = 0.000001;
  bgsd.grad2forceconvergence = 0.0001;

  // [CatLearn] //
  // No reasonable default for catl_path
  catl.path = ""s;
  catl.model = "gp"s;
  catl.prior = "median"s;
  catl.use_deriv = true;
  catl.use_fingerprint = false;
  catl.parallel = false;
}

std::string Parameters::toLowerCase(std::string s) {
  for (std::string::size_type i = 0; i < s.length(); ++i) {
    s[i] = tolower(s[i]);
  }
  return s;
}

int Parameters::load(const std::string &filename) {
  try {
    auto config = toml::parse_file(filename);

    // Main section
    main.job = magic_enum::enum_cast<JobType>(
                   config["Main"]["job"].value_or("process_search"s),
                   magic_enum::case_insensitive)
                   .value_or(JobType::Process_Search);
    main.temperature = config["Main"]["temperature"].value_or(300.0);
    main.randomSeed = config["Main"]["random_seed"].value_or(-1L);
    main.checkpoint = config["Main"]["checkpoint"].value_or(false);
    main.quiet = config["Main"]["quiet"].value_or(false);
    main.writeLog = config["Main"]["write_log"].value_or(true);
    main.inpFilename = config["Main"]["ini_filename"].value_or("config.toml"s);
    main.conFilename = config["Main"]["con_filename"].value_or("pos.con"s);
    main.finiteDifference = config["Main"]["finite_difference"].value_or(0.01);
    main.maxForceCalls = config["Main"]["max_force_calls"].value_or(0L);
    main.removeNetForce = config["Main"]["remove_net_force"].value_or(true);

    // Initialize random generator
    if (main.randomSeed < 0) {
      unsigned i = time(NULL);
      main.randomSeed = i;
      helper_functions::random(i);
    } else {
      helper_functions::random(main.randomSeed);
    }

    // Potential section
    pot.potential = magic_enum::enum_cast<PotType>(
                        config["Potential"]["potential"].value_or("LJ"s),
                        magic_enum::case_insensitive)
                        .value_or(PotType::LJ);
    pot.MPIPollPeriod = config["Potential"]["mpi_poll_period"].value_or(0.25);
    pot.LAMMPSLogging = config["Potential"]["lammps_logging"].value_or(false);
    pot.LAMMPSThreads = config["Potential"]["lammps_threads"].value_or(0);
    pot.EMTRasmussen = config["Potential"]["emt_rasmussen"].value_or(false);
    pot.extPotPath = config["Potential"]["ext_pot_path"].value_or("./ext_pot"s);
    pot.LogPotential = config["Potential"]["log_potential"].value_or(
        pot.potential == PotType::MPI || pot.potential == PotType::VASP ||
        pot.potential == PotType::BOPFOX || pot.potential == PotType::BOP);

    // AMS section
    ams.engine = config["AMS"]["engine"].value_or(""s); // One of REAXFF MOPA
    ams.forcefield =
        config["AMS"]["forcefield"].value_or(""s); // OPt.ff or something else
    ams.resources =
        config["AMS"]["resources"].value_or(""s); // PM7 PM3 or something
    ams.model =
        config["AMS"]["model"].value_or(""s); // exchange-correlation functional
    ams.xc = config["AMS"]["xc"].value_or(""s);       // with xc
    ams.basis = config["AMS"]["basis"].value_or(""s); // For DFTB

    // AMS_ENV section
    // Horrid little section to mimic amsrc.sh
    // Assumes the entire thing is going to be set
    amsenv.amshome = config["AMS_ENV"]["amshome"].value_or(
        ""s); // "/some/path/to/amshome/"s;
    amsenv.scm_tmpdir =
        config["AMS_ENV"]["scm_tmpdir"].value_or(""s); // "/tmp"s;
    amsenv.scmlicense =
        config["AMS_ENV"]["scmlicense"].value_or(""s); // "/.scm/python"s;
    amsenv.scm_pythondir = config["AMS_ENV"]["scm_pythondir"].value_or(
        ""s); // amshome.append("/bin"s);
    amsenv.amsbin = config["AMS_ENV"]["amsbin"].value_or(
        ""s); // amshome.append("license.txt"s);
    amsenv.amsresources = config["AMS_ENV"]["amsresources"].value_or(
        ""s); // amshome.append("/atomicdata"s);

    // XTBPot section
    xtbpot.paramset = config["XTBPot"]["paramset"].value_or("GFNFF"s);
    xtbpot.acc = config["XTBPot"]["accuracy"].value_or(1.0);
    xtbpot.elec_temperature =
        config["XTBPot"]["electronic_temperature"].value_or(0.0);
    xtbpot.maxiter = config["XTBPot"]["max_iterations"].value_or(250L);

    // Structure Comparison section
    structcomp.distanceDifference =
        config["Structure_Comparison"]["distance_difference"].value_or(0.1);
    structcomp.neighborCutoff =
        config["Structure_Comparison"]["neighbor_cutoff"].value_or(3.3);
    structcomp.checkRotation =
        config["Structure_Comparison"]["check_rotation"].value_or(false);
    structcomp.indistinguishableAtoms =
        config["Structure_Comparison"]["indistinguishable_atoms"].value_or(
            true);
    structcomp.energyDifference =
        config["Structure_Comparison"]["energy_difference"].value_or(0.01);
    structcomp.removeTranslation =
        config["Structure_Comparison"]["remove_translation"].value_or(true);

    // Process Search section
    procsearch.minimizeFirst =
        config["Process_Search"]["minimize_first"].value_or(true);
    procsearch.minimizationOffset =
        config["Process_Search"]["minimization_offset"].value_or(optim.maxMove);

    // Saddle Search section
    saddle.maxJumpAttempts =
        config["Saddle_Search"]["max_jump_attempts"].value_or(0L);
    saddle.maxIterations =
        config["Saddle_Search"]["max_iterations"].value_or(1000L);
    saddle.method = config["Saddle_Search"]["method"].value_or("min_mode"s);
    saddle.minmodeMethod = config["Saddle_Search"]["min_mode_method"].value_or(
        "LowestEigenmode::MINMODE_DIMER"s);
    saddle.displaceType =
        config["Saddle_Search"]["displace_type"].value_or("load"s);
    saddle.maxEnergy = config["Saddle_Search"]["max_energy"].value_or(20.0);
    saddle.displaceMagnitude =
        config["Saddle_Search"]["displace_magnitude"].value_or(0.1);
    saddle.maxSingleDisplace =
        config["Saddle_Search"]["max_single_displace"].value_or(10.0);
    saddle.displaceRadius =
        config["Saddle_Search"]["displace_radius"].value_or(4.0);
    saddle.convergedForce = config["Saddle_Search"]["converged_force"].value_or(
        optim.convergedForce);
    saddle.perpForceRatio =
        config["Saddle_Search"]["perp_force_ratio"].value_or(0.0);
    saddle.nonnegativeDisplacementAbort =
        config["Saddle_Search"]["nonnegative_displacement_abort"].value_or(
            false);
    saddle.nonlocalCountAbort =
        config["Saddle_Search"]["nonlocal_count_abort"].value_or(0L);
    saddle.nonlocalDistanceAbort =
        config["Saddle_Search"]["nonlocal_distance_abort"].value_or(0.0);
    saddle.removeRotation =
        config["Saddle_Search"]["remove_rotation"].value_or(false);
    saddle.dynamicsTemperature =
        config["Saddle_Search"]["dynamics_temperature"].value_or(0.0);
    saddle.dynamicsStateCheckIntervalInput =
        config["Saddle_Search"]["dynamics_state_check_interval"].value_or(
            100.0);
    saddle.dynamicsStateCheckInterval =
        saddle.dynamicsStateCheckIntervalInput / timeUnit;
    saddle.dynamicsRecordIntervalInput =
        config["Saddle_Search"]["dynamics_record_interval"].value_or(10.0);
    saddle.dynamicsRecordInterval =
        saddle.dynamicsRecordIntervalInput / timeUnit;
    saddle.dynamicsLinearInterpolation =
        config["Saddle_Search"]["dynamics_linear_interpolation"].value_or(true);
    saddle.dynamicsMaxInitCurvature =
        config["Saddle_Search"]["dynamics_max_init_curvature"].value_or(0.0);
    saddle.confinePositive =
        config["Saddle_Search"]["confine_positive"].value_or(false);
    saddle.bowlBreakout =
        config["Saddle_Search"]["bowl_breakout"].value_or(false);
    saddle.bowlActive = config["Saddle_Search"]["bowl_active"].value_or(20L);
    saddle.confinePositiveMinForce =
        config["Saddle_Search"]["confine_positive_min_force"].value_or(0.5);
    saddle.confinePositiveScaleRatio =
        config["Saddle_Search"]["confine_positive_scale_ratio"].value_or(0.9);
    saddle.confinePositiveBoost =
        config["Saddle_Search"]["confine_positive_boost"].value_or(10.0);
    saddle.confinePositiveMinActive =
        config["Saddle_Search"]["confine_positive_min_active"].value_or(30L);
    saddle.zeroModeAbortCurvature =
        config["Saddle_Search"]["zero_mode_abort_curvature"].value_or(0.0);

    // Optimizer section
    optim.method = magic_enum::enum_cast<OptType>(
                       config["Optimizer"]["opt_method"].value_or("none"s),
                       magic_enum::case_insensitive)
                       .value_or(OptType::CG);
    optim.refineOptMethod =
        magic_enum::enum_cast<OptType>(
            config["Optimizer"]["refine_opt"].value_or("none"s),
            magic_enum::case_insensitive)
            .value_or(OptType::CG);
    optim.refineOptMethod =
        magic_enum::enum_cast<OptType>(
            config["Refine"]["opt_method"].value_or("none"s),
            magic_enum::case_insensitive)
            .value_or(OptType::None);
    optim.refineThreshold = config["Refine"]["threshold"].value_or(0.5);
    optim.convergenceMetric = toLowerCase(
        config["Optimizer"]["convergence_metric"].value_or("norm"s));
    optim.convergenceMetricLabel =
        config["Optimizer"]["convergence_metric_label"].value_or("||Force||"s);
    optim.maxIterations = config["Optimizer"]["max_iterations"].value_or(1000L);
    optim.maxMove = config["Optimizer"]["max_move"].value_or(0.2);
    optim.convergedForce =
        config["Optimizer"]["converged_force"].value_or(0.01);
    optim.timeStepInput = config["Optimizer"]["time_step"].value_or(1.0);
    optim.timeStep = optim.timeStepInput / timeUnit;
    optim.maxTimeStepInput = config["Optimizer"]["max_time_step"].value_or(2.5);
    optim.maxTimeStep = optim.maxTimeStepInput / timeUnit;
    optim.LBFGSMemory = config["Optimizer"]["lbfgs_memory"].value_or(20L);
    optim.LBFGSInverseCurvature =
        config["Optimizer"]["lbfgs_inverse_curvature"].value_or(0.01);
    optim.LBFGSMaxInverseCurvature =
        config["Optimizer"]["lbfgs_max_inverse_curvature"].value_or(0.01);
    optim.LBFGSAutoScale =
        config["Optimizer"]["lbfgs_auto_scale"].value_or(true);
    optim.LBFGSAngleReset =
        config["Optimizer"]["lbfgs_angle_reset"].value_or(true);
    optim.LBFGSDistanceReset =
        config["Optimizer"]["lbfgs_distance_reset"].value_or(true);
    optim.QMSteepestDecent =
        config["Optimizer"]["qm_steepest_descent"].value_or(false);
    optim.CGNoOvershooting =
        config["Optimizer"]["cg_no_overshooting"].value_or(false);
    optim.CGKnockOutMaxMove =
        config["Optimizer"]["cg_knock_out_max_move"].value_or(false);
    optim.CGLineSearch = config["Optimizer"]["cg_line_search"].value_or(false);
    optim.CGLineConverged =
        config["Optimizer"]["cg_line_converged"].value_or(0.1);
    optim.CGLineSearchMaxIter =
        config["Optimizer"]["cg_max_iter_line_search"].value_or(10L);
    optim.CGMaxIterBeforeReset =
        config["Optimizer"]["cg_max_iter_before_reset"].value_or(0L);
    optim.SDAlpha = config["Optimizer"]["sd_alpha"].value_or(0.1);
    optim.SDTwoPoint = config["Optimizer"]["sd_twopoint"].value_or(false);

    // Dimer section
    dimer.rotationAngle = config["Dimer"]["finite_angle"].value_or(0.005);
    dimer.improved = config["Dimer"]["improved"].value_or(true);
    dimer.convergedAngle = config["Dimer"]["converged_angle"].value_or(5.0);
    dimer.maxIterations = config["Dimer"]["max_iterations"].value_or(1000L);
    dimer.optMethod = config["Dimer"]["opt_method"].value_or("OPT_CG"s);
    dimer.rotationsMin = config["Dimer"]["rotations_min"].value_or(1L);
    dimer.rotationsMax = config["Dimer"]["rotations_max"].value_or(10L);
    dimer.torqueMin = config["Dimer"]["torque_min"].value_or(0.1);
    dimer.torqueMax = config["Dimer"]["torque_max"].value_or(1.0);
    dimer.removeRotation = config["Dimer"]["remove_rotation"].value_or(false);

    // GPR Dimer section
    gprd.rotationAngle = config["GPR_Dimer"]["finite_angle"].value_or(0.005);
    gprd.convergedAngle = config["GPR_Dimer"]["converged_angle"].value_or(0.08);
    gprd.relaxConvAngle =
        config["GPR_Dimer"]["relaxation_converged_angle"].value_or(0.001);
    gprd.initRotationsMax =
        config["GPR_Dimer"]["max_initial_rotation_iterations"].value_or(6L);
    gprd.relaxRotationsMax =
        config["GPR_Dimer"]["max_relaxation_rotation_iterations"].value_or(10L);
    gprd.divisorTdimerGP = config["GPR_Dimer"]["divisor_t_dimer"].value_or(10L);
    gprd.maxOuterIterations =
        config["GPR_Dimer"]["max_outer_iterations"].value_or(300L);
    gprd.maxInnerIterations =
        config["GPR_Dimer"]["max_inner_iterations"].value_or(1000L);
    gprd.midpointMaxDisp =
        config["GPR_Dimer"]["max_midpoint_displacement"].value_or(0.5);
    gprd.rotOptMethod =
        config["GPR_Dimer"]["rotation_opt_method"].value_or("lbfgs"s);
    gprd.transOptMethod =
        config["GPR_Dimer"]["translation_opt_method"].value_or("lbfgs"s);
    gprd.activeRadius = config["GPR_Dimer"]["active_radius"].value_or(5.0);
    gprd.dimerSep = config["GPR_Dimer"]["dimer_separation"].value_or(0.01);
    gprd.convStep =
        config["GPR_Dimer"]["convex_region_step_size"].value_or(0.1);
    gprd.maxStep = config["GPR_Dimer"]["max_step_size"].value_or(0.1);
    gprd.forceThreshold = config["GPR_Dimer"]["force_threshold"].value_or(0.01);
    gprd.ratioAtLimit = config["GPR_Dimer"]["ratio_at_limit"].value_or(0.66667);
    gprd.initRotGP =
        config["GPR_Dimer"]["nogp_initial_rotations"].value_or(false);
    gprd.initTransGP =
        config["GPR_Dimer"]["nogp_initial_translations"].value_or(false);
    gprd.manyIterations =
        config["GPR_Dimer"]["has_many_iterations"].value_or(true);
    gprd.hyperOptMethod =
        config["GPR_Dimer"]["hyperparameter_opt_method"].value_or("scg"s);
    gprd.sigma2 = config["GPR_Dimer"]["gpr_variance"].value_or(1e-8);
    gprd.jitterSigma2 =
        config["GPR_Dimer"]["gpr_jitter_variance"].value_or(0.0);
    gprd.noiseSigma2 = config["GPR_Dimer"]["gpr_noise_variance"].value_or(1e-8);
    gprd.priorMu = config["GPR_Dimer"]["prior_mean"].value_or(0.0);
    gprd.priorSigma2 = config["GPR_Dimer"]["prior_variance"].value_or(1.0);
    gprd.priorNu =
        config["GPR_Dimer"]["prior_degrees_of_freedom"].value_or(20L);
    gprd.optCheckDerivatives =
        config["GPR_Dimer"]["check_derivative"].value_or(false);
    gprd.optMaxIterations =
        config["GPR_Dimer"]["opt_max_iterations"].value_or(400L);
    gprd.optTolFunc = config["GPR_Dimer"]["opt_tol_func"].value_or(1e-4);
    gprd.optTolSol = config["GPR_Dimer"]["opt_tol_sol"].value_or(1e-4);
    gprd.optLambdaLimit = config["GPR_Dimer"]["lambda_limit"].value_or(1e17);
    gprd.optLambdaInit = config["GPR_Dimer"]["lambda"].value_or(10.0);
    gprd.usePrune = config["GPR_Dimer"]["use_prune"].value_or(false);
    gprd.pruneBegin = config["GPR_Dimer"]["start_prune_at"].value_or(8L);
    gprd.pruneNVals = config["GPR_Dimer"]["nprune_vals"].value_or(3L);
    gprd.pruneThreshold = config["GPR_Dimer"]["prune_threshold"].value_or(0.5);
    gprd.reportLevel = config["GPR_Dimer"]["report_level"].value_or(1);
    gprd.debugLevel = config["GPR_Dimer"]["debug_level"].value_or(2);
    gprd.debugOutDir =
        config["GPR_Dimer"]["debug_output_directory"].value_or("output"s);
    gprd.debugPosFile =
        config["GPR_Dimer"]["debug_position_basename"].value_or("position"s);
    gprd.debugEnergyFile =
        config["GPR_Dimer"]["debug_energy_basename"].value_or("energy"s);
    gprd.debugGradFile =
        config["GPR_Dimer"]["debug_gradient_basename"].value_or("gradient"s);
    gprd.debugOutExt =
        config["GPR_Dimer"]["debug_output_file_extension"].value_or("dat"s);
    gprd.debugOffsetMidPoint =
        config["GPR_Dimer"]["debug_midpoint_offset"].value_or(3.0);
    gprd.debugDy = config["GPR_Dimer"]["debug_y_step"].value_or(0.1);
    gprd.debugDz = config["GPR_Dimer"]["debug_z_step"].value_or(0.1);

    // Lanczos section
    lanczos.tolerance = config["Lanczos"]["tolerance"].value_or(0.01);
    lanczos.maxIterations = config["Lanczos"]["max_iterations"].value_or(20L);
    lanczos.quitEarly = config["Lanczos"]["quit_early"].value_or(true);

    // Prefactor section
    prefactor.defaultValue = config["Prefactor"]["default_value"].value_or(0.0);
    prefactor.maxValue = config["Prefactor"]["max_value"].value_or(1e21);
    prefactor.minValue = config["Prefactor"]["min_value"].value_or(1e9);
    prefactor.withinRadius = config["Prefactor"]["within_radius"].value_or(3.3);
    prefactor.minDisplacement =
        config["Prefactor"]["min_displacement"].value_or(0.25);
    prefactor.rate =
        config["Prefactor"]["rate_estimation"].value_or("RATE_HTST"s);
    prefactor.configuration =
        magic_enum::enum_cast<PrefactorType>(
            config["Prefactor"]["configuration"].value_or("reactant"s),
            magic_enum::case_insensitive)
            .value_or(PrefactorType::REACTANT);
    prefactor.allFreeAtoms =
        config["Prefactor"]["all_free_atoms"].value_or(false);
    prefactor.filterScheme =
        config["Prefactor"]["filter_scheme"].value_or("fraction"s);
    prefactor.filterFraction =
        config["Prefactor"]["filter_fraction"].value_or(0.90);

    // Hessian section
    hessian.atomList = config["Hessian"]["atom_list"].value_or("All"s);
    hessian.zeroFreqValue = config["Hessian"]["zero_freq_value"].value_or(1e-6);

    // Nudged Elastic Band section
    neb.images = config["NEB"]["images"].value_or(5L);
    neb.maxIterations = config["NEB"]["max_iterations"].value_or(1000L);
    neb.spring = config["NEB"]["spring"].value_or(5.0);
    neb.climbingImageMethod =
        config["NEB"]["climbing_image_method"].value_or(true);
    neb.climbingImageConvergedOnly =
        config["NEB"]["climbing_image_converged_only"].value_or(true);
    neb.oldTangent = config["NEB"]["old_tangent"].value_or(false);
    neb.doublyNudged = config["NEB"]["doubly_nudged"].value_or(false);
    neb.doublyNudgedSwitching =
        config["NEB"]["doubly_nudged_switching"].value_or(false);
    neb.elasticBand = config["NEB"]["elastic_band"].value_or(false);
    neb.convergedForce =
        config["NEB"]["converged_force"].value_or(optim.convergedForce);
    neb.KSPMin = config["NEB"]["ew_ksp_min"].value_or(0.97);
    neb.KSPMax = config["NEB"]["ew_ksp_max"].value_or(9.7);
    neb.energyWeighted = config["NEB"]["energy_weighted"].value_or(false);

    // Dynamics section
    md.timeStepInput = config["Dynamics"]["time_step"].value_or(1.0);
    md.timeStep = md.timeStepInput / timeUnit;
    md.timeInput = config["Dynamics"]["time"].value_or(1000.0);
    md.time = md.timeInput / timeUnit;
    md.steps = static_cast<long>(std::floor(md.time / md.timeStep + 0.5));
    thermostat.kind = config["Dynamics"]["thermostat"].value_or("andersen"s);
    thermostat.andersenAlpha =
        config["Dynamics"]["andersen_alpha"].value_or(1.0);
    thermostat.andersenTcolInput =
        config["Dynamics"]["andersen_collision_period"].value_or(100.0);
    thermostat.andersenTcol = thermostat.andersenTcolInput / timeUnit;
    thermostat.noseMass = config["Dynamics"]["nose_mass"].value_or(1.0);
    thermostat.langevinFrictionInput =
        config["Dynamics"]["langevin_friction"].value_or(0.01);
    thermostat.langevinFriction = thermostat.langevinFrictionInput * timeUnit;

    // Parallel Replica section
    parrep.autoStop =
        config["Parallel_Replica"]["stop_after_transition"].value_or(false);
    parrep.refineTransition =
        config["Parallel_Replica"]["refine_transition"].value_or(true);
    parrep.dephaseLoopStop =
        config["Parallel_Replica"]["dephase_loop_stop"].value_or(false);
    parrep.dephaseTimeInput =
        config["Parallel_Replica"]["dephase_time"].value_or(1000.0);
    parrep.dephaseTime = parrep.dephaseTimeInput / timeUnit;
    parrep.dephaseLoopMax =
        config["Parallel_Replica"]["dephase_loop_max"].value_or(5L);
    parrep.stateCheckIntervalInput =
        config["Parallel_Replica"]["state_check_interval"].value_or(1000.0);
    parrep.stateCheckInterval = parrep.stateCheckIntervalInput / timeUnit;
    parrep.recordIntervalInput =
        config["Parallel_Replica"]["state_save_interval"].value_or(
            0.1 * parrep.stateCheckIntervalInput);
    parrep.recordInterval = parrep.recordIntervalInput / timeUnit;
    parrep.corrTimeInput =
        config["Parallel_Replica"]["post_transition_time"].value_or(1000.0);
    parrep.corrTime = parrep.corrTimeInput / timeUnit;

    // TAD section
    tad.lowT = config["TAD"]["low_temperature"].value_or(300.0);
    tad.minPrefactor = config["TAD"]["min_prefactor"].value_or(0.001);
    tad.confidence = config["TAD"]["confidence"].value_or(0.001);

    // GP Surrogate Parameters
    surrogate.use = config["Surrogate"]["use"].value_or(false);
    // If use_surrogate is true, job -> gp_surrogate and sub_job->job
    if (surrogate.use) {
      // TODO: What about other jobs
      surrogate.sub_job = main.job;
      main.job = JobType::GP_Surrogate;
    }
    surrogate.gp_uncertainity =
        config["Surrogate"]["gp_uncertainity"].value_or(0.05);
    surrogate.gp_linear_path_always =
        config["Surrogate"]["gp_linear_path_always"].value_or(false);
    if (config.contains("Surrogate"s)) {
      surrogate.potential =
          magic_enum::enum_cast<PotType>(
              config["Surrogate"]["potential"].value_or("catlearn"s),
              magic_enum::case_insensitive)
              .value_or(PotType::UNKNOWN);
      if (surrogate.potential != PotType::CatLearn) {
        throw std::runtime_error("We only support catlearn for GP right now"s);
      }
    }
    // [CatLearn]
    if (config.contains("CatLearn"s)) {
      // Case sensitive!!
      catl.path = config["CatLearn"]["path"].value_or(""s);
      catl.model = config["CatLearn"]["model"].value_or("gp"s);
      catl.prior = config["CatLearn"]["prior"].value_or("median"s);
      catl.use_deriv = config["CatLearn"]["use_derivative"].value_or(true);
      catl.parallel = config["CatLearn"]["parallel"].value_or(true);
      catl.use_deriv = config["CatLearn"]["use_fingerprint"].value_or(false);
    }
    // [ASE_ORCA]
    if (config.contains("ASE_ORCA"s)) {
      // Case sensitive!!
      // Can be used with environment variables
      aseorca.orca_path = config["ASE_ORCA"]["orca_path"].value_or(""s);
      aseorca.orca_nproc = config["ASE_ORCA"]["nproc"].value_or("1"s);
      aseorca.simpleinput = config["ASE_ORCA"]["simpleinput"].value_or(""s);
    }
    // Replica Exchange section
    repexc.temperatureDistribution =
        config["Replica_Exchange"]["temperature_distribution"].value_or(
            "exponential"s);
    repexc.replicas = config["Replica_Exchange"]["replicas"].value_or(10L);
    repexc.exchangeTrials =
        config["Replica_Exchange"]["exchange_trials"].value_or(repexc.replicas);
    repexc.samplingTimeInput =
        config["Replica_Exchange"]["sampling_time"].value_or(1000.0);
    repexc.samplingTime = repexc.samplingTimeInput / timeUnit;
    repexc.temperatureLow =
        config["Replica_Exchange"]["temperature_low"].value_or(0.0);
    repexc.temperatureHigh =
        config["Replica_Exchange"]["temperature_high"].value_or(0.0);
    repexc.exchangePeriodInput =
        config["Replica_Exchange"]["exchange_period"].value_or(100.0);
    repexc.exchangePeriod = repexc.exchangePeriodInput / timeUnit;

    // Hyperdynamics section
    bondBoost.biasPotential =
        config["Hyperdynamics"]["bias_potential"].value_or("NONE"s);
    bondBoost.BALS =
        config["Hyperdynamics"]["bb_boost_atomlist"].value_or("ALL"s);
    bondBoost.RMDTimeInput =
        config["Hyperdynamics"]["bb_rmd_time"].value_or(100.0);
    bondBoost.RMDTime = bondBoost.RMDTimeInput / timeUnit;
    bondBoost.DVMAX = config["Hyperdynamics"]["bb_dvmax"].value_or(0.0);
    bondBoost.QRR =
        config["Hyperdynamics"]["bb_stretch_threshold"].value_or(0.2);
    bondBoost.PRR = config["Hyperdynamics"]["bb_ds_curvature"].value_or(0.95);
    bondBoost.Qcut = config["Hyperdynamics"]["bb_rcut"].value_or(3.0);

    // Basin Hopping section
    bhop.displacement = config["Basin_Hopping"]["displacement"].value_or(0.5);
    bhop.initialRandomStructureProbability =
        config["Basin_Hopping"]["initial_random_structure_probability"]
            .value_or(0.0);
    bhop.pushApartDistance =
        config["Basin_Hopping"]["push_apart_distance"].value_or(0.4);
    bhop.steps = config["Basin_Hopping"]["steps"].value_or(10000L);
    bhop.quenchingSteps =
        config["Basin_Hopping"]["quenching_steps"].value_or(0L);
    bhop.significantStructure =
        config["Basin_Hopping"]["significant_structure"].value_or(true);
    bhop.singleAtomDisplace =
        config["Basin_Hopping"]["single_atom_displace"].value_or(false);
    bhop.displacementAlgorithm =
        config["Basin_Hopping"]["displacement_algorithm"].value_or("standard"s);
    bhop.displacementDistribution =
        config["Basin_Hopping"]["displacement_distribution"].value_or(
            "uniform"s);
    bhop.swapProbability =
        config["Basin_Hopping"]["swap_probability"].value_or(0.0);
    bhop.jumpMax = config["Basin_Hopping"]["jump_max"].value_or(10L);
    bhop.jumpSteps = config["Basin_Hopping"]["jump_steps"].value_or(0L);
    bhop.adjustDisplacement =
        config["Basin_Hopping"]["adjust_displacement"].value_or(true);
    bhop.adjustPeriod = config["Basin_Hopping"]["adjust_period"].value_or(10L);
    bhop.adjustFraction =
        config["Basin_Hopping"]["adjust_fraction"].value_or(0.05);
    bhop.targetRatio = config["Basin_Hopping"]["target_ratio"].value_or(0.5);
    bhop.writeUnique = config["Basin_Hopping"]["write_unique"].value_or(false);
    bhop.stopEnergy = config["Basin_Hopping"]["stop_energy"].value_or(-DBL_MAX);

    // Global Optimization section
    globopt.moveMethod =
        config["Global_Optimization"]["move_method"].value_or("md"s);
    globopt.decisionMethod =
        config["Global_Optimization"]["decision_method"].value_or("npew"s);
    globopt.steps = config["Global_Optimization"]["steps"].value_or(10000L);
    globopt.beta = config["Global_Optimization"]["beta"].value_or(1.05);
    globopt.alpha = config["Global_Optimization"]["alpha"].value_or(1.02);
    globopt.mdmin = config["Global_Optimization"]["mdmin"].value_or(3L);
    globopt.targetEnergy =
        config["Global_Optimization"]["target_energy"].value_or(-1.E50);

    // Monte Carlo section
    monte_carlo.stepSize = config["Monte_Carlo"]["step_size"].value_or(0.005);
    monte_carlo.steps = config["Monte_Carlo"]["steps"].value_or(1000);

    // BGSD section
    bgsd.alpha = config["BGSD"]["alpha"].value_or(10.0);
    bgsd.beta = config["BGSD"]["beta"].value_or(0.2);
    bgsd.gradientfinitedifference =
        config["BGSD"]["gradientfinitedifference"].value_or(0.000001);
    bgsd.Hforceconvergence = config["BGSD"]["Hforceconvergence"].value_or(0.01);
    bgsd.grad2energyconvergence =
        config["BGSD"]["grad2energyconvergence"].value_or(0.000001);
    bgsd.grad2forceconvergence =
        config["BGSD"]["grad2forceconvergence"].value_or(0.0001);

    // Debug section
    debug.writeMovies = config["Debug"]["write_movies"].value_or(false);
    debug.writeMoviesInterval =
        config["Debug"]["write_movies_interval"].value_or(1L);

    // Sanity Checks
    if (parrep.stateCheckInterval > md.time &&
        magic_enum::enum_name<JobType>(main.job) == "parallel_replica"s) {
      spdlog::error("[Parallel Replica] state_check_interval must be <= time"s);
      return 1;
    }

    if (saddle.dynamicsRecordIntervalInput >
        saddle.dynamicsStateCheckIntervalInput) {
      spdlog::error("[Saddle Search] dynamics_record_interval must be <= "
                    "dynamics_state_check_interval"s);
      return 1;
    }

    if (pot.potential == PotType::AMS || pot.potential == PotType::AMS_IO) {
      if (ams.forcefield.empty() && ams.model.empty() && ams.xc.empty()) {
        spdlog::error("[AMS] Must provide at least forcefield or model or xc"s);
        return 1;
      }

      if (!ams.forcefield.empty() && !ams.model.empty() && !ams.xc.empty()) {
        spdlog::error("[AMS] Must provide either forcefield or model"s);
        return 1;
      }
    }

  } catch (const std::exception &e) {
    spdlog::error("Error parsing the configuration file: {}", e.what());
    return 1;
  }

  return 0;
}

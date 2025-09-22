#include "Parameters.h"
#include "BaseStructures.h"
#include "BondBoost.h"
#include "Dynamics.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "INIFile.h"
#include "ImprovedDimer.h"
#include "Prefactor.h"
#include "PrefactorJob.h"
#include "magic_enum/magic_enum.hpp"
#include <errno.h>
#include <float.h>
#include <stdexcept>
#include <time.h>

Parameters::Parameters() {

  kB = 8.6173324e-5;     // eV/K
  timeUnit = 10.1805055; // fs

  // [Main] //
  job = JobType::Process_Search;
  randomSeed = -1;
  temperature = 300.0;
  checkpoint = false;
  quiet = false;
  writeLog = true;
  iniFilename = "config.ini";
  conFilename = "pos.con";
  finiteDifference = 0.01;
  maxForceCalls = 0;
  removeNetForce = true;

  // [Prefactor] //
  prefactorDefaultValue = 0.0;
  prefactorMaxValue = 1e+21;
  prefactorMinValue = 1e+9;
  prefactorWithinRadius = 3.3;
  prefactorMinDisplacement = 0.25;
  prefactorRate = Prefactor::RATE_HTST;
  prefactorConfiguration = PrefactorJob::PREFACTOR_REACTANT;
  prefactorAllFreeAtoms = false;
  prefactorFilterScheme = Prefactor::FILTER_FRACTION;
  prefactorFilterFraction = 0.90;

  // [Potential] //
  potential = PotType::LJ;
  MPIPollPeriod = 0.25; // seconds
  MPIPotentialRank = -1;
  LogPotential = false;
  LAMMPSLogging = false;
  LAMMPSThreads = 0;
  EMTRasmussen = false;
  extPotPath = "./ext_pot";

  // [AMS] //
  engine = "";     // One of REAXFF MOPAC
  forcefield = ""; // OPt.ff or something else
  model = "";      // PM7 PM3 or something
  xc = "";         // exchange-correlation functional
  basis = "";      // with xc
  resources = "";  // For DFTB

  // [AMS_ENV] //
  // Horrid little section to mimic amsrc.sh
  // Assumes the entire thing is going to be set
  amshome = "";       // "/some/path/to/amshome/";
  scm_tmpdir = "";    // "/tmp";
  scm_pythondir = ""; // "/.scm/python";
  amsbin = "";        // amshome.append("/bin");
  scmlicense = "";    // amshome.append("license.txt");
  amsresources = "";  // amshome.append("/atomicdata");

  // [XTBPot] //
  xtb_paramset = "GFNFF";
  xtb_acc = 1.0;
  xtb_elec_temperature = 0.0;
  xtb_maxiter = 250;

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
  distanceDifference = 0.1;
  neighborCutoff = 3.3;
  checkRotation = false;
  indistinguishableAtoms = true;
  energyDifference = 0.01;
  removeTranslation = true;

  // [Debug] //
  writeMovies = false;
  writeMoviesInterval = 1;
  estNEBeig = false;
  nebMMF = LowestEigenmode::MINMODE_DIMER;

  // [Saddle Search] //
  saddleDisplaceType = EpiCenters::DISP_LOAD;
  saddleMethod = "min_mode";
  saddleMinmodeMethod = LowestEigenmode::MINMODE_DIMER;
  saddleMaxEnergy = 20.0;
  saddleMaxIterations = 1000;
  saddleDisplaceRadius = 4.0;
  saddleDisplaceMagnitude = 0.1;
  saddleMaxSingleDisplace = 10.;
  //    saddleConvergedForce = optConvergedForce; default value is set after the
  //    value of optConvergedForce is loaded
  saddleNonnegativeDisplacementAbort = false;
  saddleNonlocalCountAbort = 0;
  saddleNonlocalDistanceAbort = 0.0;
  saddleRemoveRotation = false;
  saddlePerpForceRatio = 0.0;    // undocumented
  saddleConfinePositive = false; // undocumented
  saddleBowlBreakout = false;    // undocumented
  saddleBowlActive = 20;
  saddleConfinePositiveMinForce = 0.5;           // undocumented
  saddleConfinePositiveScaleRatio = 0.9;         // undocumented
  saddleConfinePositiveBoost = 10.;              // undocumented
  saddleConfinePositiveMinActive = 30;           // undocumented
  saddleDynamicsTemperature = 0.0;               // defaults to temperature
  saddleDynamicsStateCheckIntervalInput = 100.0; // fs
  saddleDynamicsStateCheckInterval =
      saddleDynamicsStateCheckIntervalInput / timeUnit;
  saddleDynamicsRecordIntervalInput = 10.0; // fs
  saddleDynamicsLinearInterpolation = true;
  saddleDynamicsMaxInitCurvature = 0.0; // eV/Ang^2
  saddleZeroModeAbortCurvature = 0.0;   // eV/Ang^2

  // [Optimizers] //
  optMethod = OptType::CG;
  optConvergenceMetric = "norm"s;
  refineOptMethod = OptType::None;
  refineThreshold = 0.5;
  optMaxIterations = 1000;
  optConvergedForce = 0.01;
  optMaxMove = 0.2;
  optTimeStepInput = 1.0;
  optMaxTimeStepInput = 2.5;
  optMaxTimeStep = optMaxTimeStepInput / timeUnit;

  optLBFGSMemory = 20;
  optLBFGSInverseCurvature =
      0.01; // assumes stiffest curvature at minimum is 100 eV/A^2
  optLBFGSAutoScale = true;
  optLBFGSAngleReset = true;
  optLBFGSDistanceReset = true;

  optQMSteepestDecent = false;
  optCGNoOvershooting = false;
  optCGKnockOutMaxMove = false;
  optCGLineConverged = 0.1;
  optCGLineSearch = false;
  optCGMaxIterBeforeReset = 0;
  optCGLineSearchMaxIter = 10;
  optSDAlpha = 0.1;
  optSDTwoPoint = false;

  // [Process Search] //
  processSearchMinimizeFirst = true;
  processSearchMinimizationOffset = optMaxMove;

  // [Dimer] //
  dimerRotationAngle = 0.005;
  dimerImproved = true;
  dimerConvergedAngle = 5.0; // degrees
  dimerMaxIterations = 1000;
  dimerOptMethod = ImprovedDimer::OPT_CG;
  dimerTorqueMin = 0.1;   // old dimer
  dimerTorqueMax = 1.0;   // old dimer
  dimerRotationsMin = 1;  // old dimer
  dimerRotationsMax = 10; // old dimer and new dimer
  dimerRemoveRotation = false;

  // [ASE_ORCA] //
  orca_path = ""s;
  orca_nproc = "1"s;

  // [ASE_NWCHEM] //
  nwchem_path = ""s;
  nwchem_nproc = "1"s;
  nwchem_multiplicity = ""s;
  nwchem_scf_thresh = 1e-5;
  nwchem_scf_maxiter = 200;

  // [Metatomic] //
  metatomic_options.model_path = ""s;
  metatomic_options.device = "cpu"s;
  metatomic_options.length_unit = "angstrom"s;
  metatomic_options.extensions_directory = ""s;
  metatomic_options.check_consistency = false;

  // [Lanczos] //
  lanczosTolerance = 0.01;
  lanczosMaxIterations = 20;
  lanczosQuitEarly = true;

  // [GPR Dimer] //
  gprDimerRotationAngle = 0.005;
  gprDimerConvergedAngle = 0.08;            // T_anglerot_init
  gprDimerRelaxConvAngle = 0.001;           // T_anglerot_gp
  gprDimerInitRotationsMax = 6;             // num_iter_initrot; should be DoF
  gprDimerRelaxRotationsMax = 10;           // num_iter_rot_gp
  gprDimerDivisorTdimerGP = 10;             // divisor_T_dimer_gp
  gprDimerMaxOuterIterations = 300;         // num_bigiter
  gprDimerMaxInnerIterations = 1000;        // num_iter
  gprDimerMidpointMaxDisp = 0.5;            // disp_max
  gprDimerRotOptMethod = "lbfgs";           // method_rot
  gprDimerTransOptMethod = "lbfgs";         // method_trans
  gprActiveRadius = 5.0;                    // actidst_fro
  gprDimerSep = 0.01;                       // dimer_sep
  gprd_trans_options.step_length = 0.1;     // param_trans[0]
  gprd_trans_options.max_step_length = 0.1; // param_trans[1]
  gprd_trans_options.rotrem_thresh = std::numeric_limits<
      double>::infinity();        // rotation_removal_projection_threshold
  gprDimerRatioAtLimit = 0.66667; // ratio_at_limit
  gprDimerInitRotGP = 0;          // initrot_nogp
  gprDimerManyIterations = true;  // islarge_num_iter
  // GPR Params
  gprDimerSigma2 = 1e-8;      // gp_sigma2
  gprDimerJitterSigma2 = 0;   // jitter_sigma2
  gprDimerNoiseSigma2 = 1e-8; // sigma2
  gprDimerPriorMu = 0;        // prior_mu
  gprDimerPriorSigma2 = 1;    // prior_s2
  gprDimerPriorNu = 20;       // prior_nu
  // GPR Optimization Parameters
  gpr_hypopt_options.hopt_method = "Ceres_opt"s; // optimization_alg
  gpr_hypopt_options.check_derivative = false; // check_derivative
  gpr_hypopt_options.tol_func = 1e-4;
  gpr_hypopt_options.tol_sol = 1e-4;
  gpr_hypopt_options.max_iter = 400;
  gpr_hypopt_options.scg.lambda_limit = 1e17;
  gpr_hypopt_options.scg.lambda = 10.0;
  gpr_hypopt_options.adam.lr = 1e-3;
  gpr_hypopt_options.adam.lrd = 0.999;
  gpr_hypopt_options.adam.b1 = 0.9;
  gpr_hypopt_options.adam.b2 = 0.99;
  gpr_hypopt_options.adam.eps = 1.0e-8;
  gpr_hypopt_options.adam.weight_decay = 0.0;
  gpr_hypopt_options.adam.amsgrad = true;
  early_stopping_options.dist_metrics = "EMD"s;
  early_stopping_options.threshold = 1.2;
  fps_options.history = 5;
  fps_options.metric = "EMD"s;
  // GPR Prune Parameters
  gprUsePrune = false;     // use_prune
  gprPruneBegin = 8;       // start_prune_at
  gprPruneNVals = 3;       // nprune_vals
  gprPruneThreshold = 0.5; // prune_threshold
  // GPR Debugging Parameters
  gprReportLevel = 1;            // report_level
  gprDebugLevel = 2;             // debug_level
  gprDebugOutDir = "output";     // debug_output_dir
  gprDebugPosFile = "position";  // debug_output_file_R
  gprDebugEnergyFile = "energy"; // debug_output_file_E
  gprDebugGradFile = "gradient"; // debug_output_file_G
  gprDebugOutExt = "dat";        // debug_output_file_extension
  gprDebugOffsetMidPoint = 3.;   // debug_offset_from_mid_point
  gprDebugDy = 0.1;              // debug_dy
  gprDebugDz = 0.1;              // debug_dz

  // GP Surrogate Parameters
  use_surrogate = false;
  sub_job = JobType::Unknown;
  gp_uncertainity = 0.05;
  gp_linear_path_always = false;
  surrogatePotential = PotType::CatLearn;

  // [Hessian] //
  hessianAtomList = string("All");
  hessianZeroFreqValue = 1e-6;

  // [Nudged Elastic Band] //
  nebImages = 5;
  nebSpring = 5.0;
  nebClimbingImageMethod = true;
  nebClimbingImageConvergedOnly = true;
  nebOldTangent = false;
  nebMaxIterations = 1000;
  nebDoublyNudged = false;
  nebDoublyNudgedSwitching = false;
  nebElasticBand = false;
  nebConvergedForce = optConvergedForce;
  nebEnergyWeighted = false;
  nebKSPMin = 0.97;
  nebKSPMax = 9.7;
  nebIpath = ""s;
  nebMinimEP = true;
  nebciAfter = std::numeric_limits<double>::infinity();
  nebciWithMMF = false;
  nebciMMFAfter = 0.5;
  nebciMMFnSteps = 10;

  // [Dynamics] //
  mdTimeStepInput = 1.0;
  mdTimeInput = 1000.0;
  mdTimeStep = mdTimeStepInput / timeUnit;
  mdTime = mdTimeInput / timeUnit;
  mdSteps = long(floor(mdTime / mdTimeStep + 0.5));

  // [Thermostat] //
  thermostat = Dynamics::NONE;
  thermoAndersenAlpha = 1.0;       // collision strength
  thermoAndersenTcolInput = 100.0; // collision frequency in unit of fs
  thermoNoseMass = 1.0;
  thermoLangevinFrictionInput = 0.01;
  thermoLangevinFriction = thermoLangevinFrictionInput * timeUnit;

  // [Parallel Replica] //
  parrepRefineTransition = true;
  parrepAutoStop = false;
  parrepDephaseLoopStop = false;
  parrepDephaseTimeInput = 1000.0;
  parrepDephaseLoopMax = 5;
  parrepStateCheckIntervalInput = 1000.0;
  parrepRecordIntervalInput = 50.0;
  parrepCorrTimeInput = 1000.0;
  parrepDephaseTime = parrepDephaseTimeInput / timeUnit;
  parrepStateCheckInterval = parrepStateCheckIntervalInput / timeUnit;
  parrepRecordInterval = parrepRecordIntervalInput / timeUnit;
  parrepCorrTime = parrepCorrTimeInput / timeUnit;

  // [Temperature Accelerated Dynamics] //
  tadLowT = 300.0;
  tadMinPrefactor = 0.001; // in unit of fs-1
  tadConfidence = 0.001;

  // [Replica Exchange] //
  repexcTemperatureDistribution = "exponential";
  repexcReplicas = 10;
  repexcExchangeTrials = repexcReplicas;
  repexcSamplingTimeInput = 1000.0;
  repexcSamplingTime = repexcSamplingTimeInput / timeUnit;
  repexcTemperatureLow = 0.0;
  repexcTemperatureHigh = 0.0;
  repexcExchangePeriod = 100.0;

  // [Hyperdynamics] //
  biasPotential = Hyperdynamics::NONE;
  bondBoostBALS = string("All"); // boosted atom list string
  bondBoostDVMAX = 0.0;
  bondBoostQRR = 0.2; // can not be set to 0
  bondBoostPRR = 0.95;
  bondBoostQcut = 3.0;
  bondBoostRMDTimeInput = 100.0;
  bondBoostRMDTime = bondBoostRMDTimeInput / timeUnit;

  // [Basin Hopping] //
  basinHoppingDisplacement = 0.5;
  basinHoppingPushApartDistance = 0.4;
  basinHoppingInitialRandomStructureProbability = 0.0;
  basinHoppingSteps = 10000;
  basinHoppingQuenchingSteps = 0;
  basinHoppingSingleAtomDisplace = false;
  basinHoppingSignificantStructure = true;
  basinHoppingDisplacementAlgorithm = "standard";
  basinHoppingDisplacementDistribution = "uniform";
  basinHoppingSwapProbability = 0.0;
  basinHoppingJumpMax = 10;
  basinHoppingJumpSteps = 0;
  basinHoppingAdjustDisplacement = true;
  basinHoppingAdjustPeriod = 10;
  basinHoppingAdjustFraction = 0.05;
  basinHoppingTargetRatio = 0.5;
  basinHoppingWriteUnique = false;
  basinHoppingStopEnergy = -DBL_MAX;

  // [Global Optimization] //
  globalOptimizationMoveMethod = "md";
  globalOptimizationDecisionMethod = "npew";
  globalOptimizationSteps = 10000;
  globalOptimizationBeta = 1.05;
  globalOptimizationAlpha = 1.02;
  globalOptimizationMdmin = 3;
  globalOptimizationTargetEnergy = -1.E50;

  // [Monte Carlo] //
  monteCarloStepSize = 0.005;
  monteCarloSteps = 1000;

  // [BGSD] //
  alpha = 10.0;
  beta = 0.2;
  gradientfinitedifference = 0.000001;
  Hforceconvergence = 0.01;
  grad2energyconvergence = 0.000001;
  grad2forceconvergence = 0.0001;

  // [CatLearn] //
  // No reasonable default for catl_path
  catl_model = "gp";
  catl_prior = "median";
  catl_use_deriv = true;
  catl_use_fingerprint = false;
  catl_parallel = false;
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

    job = magic_enum::enum_cast<JobType>(ini.GetValue("Main", "job"),
                                         magic_enum::case_insensitive)
              .value_or(JobType::Unknown);
    temperature = ini.GetValueF("Main", "temperature", temperature);
    randomSeed = ini.GetValueL("Main", "random_seed", randomSeed);
    checkpoint = ini.GetValueB("Main", "checkpoint", checkpoint);
    quiet = ini.GetValueB("Main", "quiet", quiet);
    writeLog = ini.GetValueB("Main", "write_log", writeLog);
    finiteDifference =
        ini.GetValueF("Main", "finite_difference", finiteDifference);
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

    potential =
        magic_enum::enum_cast<PotType>(ini.GetValue("Potential", "potential"),
                                       magic_enum::case_insensitive)
            .value_or(PotType::UNKNOWN);
    MPIPollPeriod =
        ini.GetValueF("Potential", "mpi_poll_period", MPIPollPeriod);
    LAMMPSLogging = ini.GetValueB("Potential", "lammps_logging", LAMMPSLogging);
    LAMMPSThreads =
        (int)ini.GetValueL("Potential", "lammps_threads", LAMMPSThreads);
    EMTRasmussen = ini.GetValueB("Potential", "emt_rasmussen", EMTRasmussen);
    extPotPath = ini.GetValue("Potential", "ext_pot_path", extPotPath);

    if (potential == PotType::MPI || potential == PotType::VASP ||
        potential == PotType::BOPFOX || potential == PotType::BOP) {
      LogPotential = true;
    } else {
      LogPotential = false;
    }
    LogPotential = ini.GetValueB("Potential", "log_potential", LogPotential);

    // [AMS]
    if (potential == PotType::AMS) {
      engine = ini.GetValue("AMS", "engine", engine);
      forcefield = ini.GetValue("AMS", "forcefield", forcefield);
      resources = ini.GetValue("AMS", "resources", resources);
      model = ini.GetValue("AMS", "model", model);
      xc = ini.GetValue("AMS", "xc", xc);
      basis = ini.GetValue("AMS", "basis", basis);
    }
    // [AMS_IO]
    if (potential == PotType::AMS_IO) {
      engine = ini.GetValue("AMS_IO", "engine", engine);
      forcefield = ini.GetValue("AMS_IO", "forcefield", forcefield);
      model = ini.GetValue("AMS_IO", "model", model);
      xc = ini.GetValue("AMS_IO", "xc", xc);
    }
    // [AMS_ENV]
    // This is only needed if the regular calls do not work
    // e.g. on a MacOS machine
    if (potential == PotType::AMS_IO || potential == PotType::AMS) {
      amshome = ini.GetValue("AMS_ENV", "amshome", amshome);
      scm_tmpdir = ini.GetValue("AMS_ENV", "scm_tmpdir", scm_tmpdir);
      scmlicense = ini.GetValue("AMS_ENV", "scmlicense", scmlicense);
      scm_pythondir = ini.GetValue("AMS_ENV", "scm_pythondir", scm_pythondir);
      amsbin = ini.GetValue("AMS_ENV", "amsbin", amsbin);
      amsresources = ini.GetValue("AMS_ENV", "amsresources", amsresources);
    }
    // [XTBPot]
    if (potential == PotType::XTB) {
      xtb_paramset = ini.GetValue("XTBPot", "paramset", xtb_paramset);
      xtb_acc = ini.GetValueF("XTBPot", "accuracy", xtb_acc);
      xtb_elec_temperature = ini.GetValueF("XTBPot", "electronic_temperature",
                                           xtb_elec_temperature);
      xtb_maxiter = ini.GetValueL("XTBPot", "max_iterations", xtb_maxiter);
    }
    // [ZBLPot]
    if (potential == PotType::ZBL) {
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
    if (potential == PotType::SocketNWChem) {
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

    writeMovies = ini.GetValueB("Debug", "write_movies", writeMovies);
    writeMoviesInterval =
        ini.GetValueL("Debug", "write_movies_interval", writeMoviesInterval);
    estNEBeig = ini.GetValueB("Debug", "estimate_neb_eigenvalues", estNEBeig);
    nebMMF = toLowerCase(ini.GetValue("Debug", "neb_mmf_estimator", nebMMF));

    // [Structure Comparison] //

    distanceDifference = ini.GetValueF(
        "Structure Comparison", "distance_difference", distanceDifference);
    neighborCutoff = ini.GetValueF("Structure Comparison", "neighbor_cutoff",
                                   neighborCutoff);
    checkRotation =
        ini.GetValueB("Structure Comparison", "check_rotation", checkRotation);
    energyDifference = ini.GetValueF("Structure Comparison",
                                     "energy_difference", energyDifference);
    indistinguishableAtoms =
        ini.GetValueB("Structure Comparison", "indistinguishable_atoms",
                      indistinguishableAtoms);
    removeTranslation = ini.GetValueB("Structure Comparison",
                                      "remove_translation", removeTranslation);

    // [Process Search] //

    processSearchMinimizeFirst = ini.GetValueB(
        "Process Search", "minimize_first", processSearchMinimizeFirst);
    processSearchMinimizationOffset =
        ini.GetValueF("Process Search", "minimization_offset",
                      processSearchMinimizationOffset);

    // [Optimizers] //
    auto inp_optMethod = magic_enum::enum_cast<OptType>(
                             ini.GetValue("Optimizer", "opt_method", "none"),
                             magic_enum::case_insensitive)
                             .value_or(OptType::Unknown);
    if (inp_optMethod != OptType::None) {
      optMethod = inp_optMethod;
    }

    optConvergenceMetric = toLowerCase(
        ini.GetValue("Optimizer", "convergence_metric", optConvergenceMetric));
    if (optConvergenceMetric == "max_atom") {
      optConvergenceMetricLabel = "Max atom force";
    } else if (optConvergenceMetric == "max_component") {
      optConvergenceMetricLabel = "Max force comp";
    } else if (optConvergenceMetric == "norm") {
      // optConvergenceMetricLabel = "\u2016Force\u2016";
      optConvergenceMetricLabel = "||Force||";
    } else {
      fprintf(stderr, "unknown convergence_metric %s\n",
              optConvergenceMetric.c_str());
      exit(1);
    }

    if (ini.FindKey("Refine") != -1) {
      refineOptMethod =
          magic_enum::enum_cast<OptType>(ini.GetValue("Refine", "opt_method"),
                                         magic_enum::case_insensitive)
              .value_or(OptType::None);
      refineThreshold = ini.GetValueF("Refine", "threshold", refineThreshold);
    }

    optConvergedForce =
        ini.GetValueF("Optimizer", "converged_force", optConvergedForce);
    optMaxIterations = static_cast<size_t>(
        ini.GetValueL("Optimizer", "max_iterations", optMaxIterations));
    optMaxMove = ini.GetValueF("Optimizer", "max_move", optMaxMove);
    processSearchMinimizationOffset = optMaxMove;
    // Handle each optimizer separately
    if (ini.FindKey("QuickMin") != -1) {
      optTimeStepInput =
          ini.GetValueF("QuickMin", "time_step", optTimeStepInput);
      optTimeStep = optTimeStepInput / timeUnit;
      optQMSteepestDecent = ini.GetValueB("Optimizer", "qm_steepest_descent",
                                          optQMSteepestDecent);
    }
    if (ini.FindKey("FIRE") != -1) {
      SPDLOG_WARN("Overwriting QuickMin timestep with Fire timestep!!");
      optTimeStepInput = ini.GetValueF("FIRE", "time_step", optTimeStepInput);
      optTimeStep = optTimeStepInput / timeUnit;
      optMaxTimeStepInput =
          ini.GetValueF("FIRE", "time_step_max", optMaxTimeStepInput);
      optMaxTimeStep = optMaxTimeStepInput / timeUnit;
    }
    if (ini.FindKey("LBFGS") != -1) {
      optLBFGSMemory = ini.GetValueL("LBFGS", "lbfgs_memory", optLBFGSMemory);
      optLBFGSInverseCurvature = ini.GetValueF(
          "LBFGS", "lbfgs_inverse_curvature", optLBFGSInverseCurvature);
      optLBFGSMaxInverseCurvature = ini.GetValueF(
          "LBFGS", "lbfgs_max_inverse_curvature", optLBFGSMaxInverseCurvature);
      optLBFGSAutoScale =
          ini.GetValueB("LBFGS", "lbfgs_auto_scale", optLBFGSAutoScale);
      optLBFGSAngleReset =
          ini.GetValueB("LBFGS", "lbfgs_angle_reset", optLBFGSAngleReset);
      optLBFGSDistanceReset =
          ini.GetValueB("LBFGS", "lbfgs_distance_reset", optLBFGSDistanceReset);
    }
    if (ini.FindKey("CG") != -1) {
      optCGNoOvershooting =
          ini.GetValueB("CG", "cg_no_overshooting", optCGNoOvershooting);
      optCGKnockOutMaxMove =
          ini.GetValueB("CG", "cg_knock_out_max_move", optCGKnockOutMaxMove);
      optCGLineSearch = ini.GetValueB("CG", "cg_line_search", optCGLineSearch);
      optCGLineConverged =
          ini.GetValueF("CG", "cg_line_converged", optCGLineConverged);
      optCGMaxIterBeforeReset = ini.GetValueL("CG", "cg_max_iter_before_reset",
                                              optCGMaxIterBeforeReset);
      optCGLineSearchMaxIter = ini.GetValueL("CG", "cg_max_iter_line_search",
                                             optCGLineSearchMaxIter);
    }
    if (ini.FindKey("SD") != -1) {
      optSDAlpha = ini.GetValueF("SD", "sd_alpha", optSDAlpha);
      optSDTwoPoint = ini.GetValueB("SD", "sd_twopoint", optSDTwoPoint);
    }

    // [Dimer] //

    dimerRotationAngle =
        ini.GetValueF("Dimer", "finite_angle", dimerRotationAngle);
    dimerImproved = ini.GetValueB("Dimer", "improved", dimerImproved);
    dimerConvergedAngle =
        ini.GetValueF("Dimer", "converged_angle", dimerConvergedAngle);
    dimerMaxIterations =
        ini.GetValueL("Dimer", "max_iterations", dimerMaxIterations);
    dimerOptMethod =
        toLowerCase(ini.GetValue("Dimer", "opt_method", dimerOptMethod));
    dimerRotationsMin =
        ini.GetValueL("Dimer", "rotations_min", dimerRotationsMin); // old
    dimerRotationsMax = ini.GetValueL("Dimer", "rotations_max",
                                      dimerRotationsMax); // old & new
    dimerTorqueMin =
        ini.GetValueF("Dimer", "torque_min", dimerTorqueMin); // old
    dimerTorqueMax =
        ini.GetValueF("Dimer", "torque_max", dimerTorqueMax); // old
    dimerRemoveRotation =
        ini.GetValueB("Dimer", "remove_rotation", dimerRemoveRotation);

    // GP Surrogate Parameters
    use_surrogate = ini.GetValueB("Surrogate", "use_surrogate", false);
    // If use_surrogate is true, job -> gp_surrogate and sub_job->job
    if (use_surrogate) {
      // TODO: What about other jobs
      sub_job = job;
      job = JobType::GP_Surrogate;
    }
    gp_uncertainity =
        ini.GetValueF("Surrogate", "gp_uncertainity", gp_uncertainity);
    if (ini.FindKey("Surrogate") != -1) {
      surrogatePotential =
          magic_enum::enum_cast<PotType>(ini.GetValue("Surrogate", "potential"),
                                         magic_enum::case_insensitive)
              .value_or(PotType::UNKNOWN);
      if (surrogatePotential != PotType::CatLearn) {
        throw std::runtime_error("We only support catlearn for GP right now");
      }
    }
    // [CatLearn]
    if (ini.FindKey("CatLearn") != -1) {
      // Case sensitive!!
      catl_path = ini.GetValue("CatLearn", "catl_path");
      catl_model = ini.GetValue("CatLearn", "model", "catl_model");
      catl_prior = ini.GetValue("CatLearn", "prior", "catl_prior");
      catl_use_deriv =
          ini.GetValueB("CatLearn", "use_derivatives", "catl_deriv");
      catl_use_fingerprint =
          ini.GetValueB("CatLearn", "use_fingerprint", "catl_fingerprint");
      catl_parallel = ini.GetValueB("CatLearn", "parallel_hyperparameter_opt",
                                    "catl_parallel");
    }
    // [ASE_ORCA]
    if (ini.FindKey("ASE_ORCA") != -1) {
      // Case sensitive!!
      // TODO: This should be handled in clienteon so you can still call
      // eonclient for single point calculations easily
      orca_path = ini.GetValue("ASE_ORCA", "orca_path", "");
      orca_nproc = ini.GetValue("ASE_ORCA", "nproc");
      orca_sline = ini.GetValue("ASE_ORCA", "simpleinput", "");
    }
    // [ASE_NWCHEM]
    if (ini.FindKey("ASE_NWCHEM") != -1) {
      // Case sensitive!!
      // TODO: This should be handled in clienteon so you can still call
      // eonclient for single point calculations easily
      nwchem_path = ini.GetValue("ASE_NWCHEM", "nwchem_path", "");
      nwchem_nproc = ini.GetValue("ASE_NWCHEM", "nproc");
      nwchem_multiplicity = ini.GetValue("ASE_NWCHEM", "multiplicity");
      nwchem_scf_thresh = ini.GetValueF("ASE_NWCHEM", "scf_thresh", 1e-5);
      nwchem_scf_maxiter = ini.GetValueL("ASE_NWCHEM", "scf_maxiter", 200);
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
    }
    // GP_NEB only
    gp_linear_path_always = ini.GetValueB("Surrogate", "gp_linear_path_always",
                                          gp_linear_path_always);
    // [Lanczos] //

    lanczosTolerance = ini.GetValueF("Lanczos", "tolerance", lanczosTolerance);
    lanczosMaxIterations =
        ini.GetValueL("Lanczos", "max_iterations", lanczosMaxIterations);
    lanczosQuitEarly = ini.GetValueB("Lanczos", "quit_early", lanczosQuitEarly);

    // [GPR Dimer] //
    gprDimerRotationAngle =
        ini.GetValueF("GPR Dimer", "finite_angle", gprDimerRotationAngle);
    gprDimerConvergedAngle =
        ini.GetValueF("GPR Dimer", "converged_angle", gprDimerConvergedAngle);
    gprDimerRelaxConvAngle = ini.GetValueF(
        "GPR Dimer", "relaxation_converged_angle", gprDimerRelaxConvAngle);
    gprDimerInitRotationsMax =
        (int)ini.GetValueL("GPR Dimer", "max_initial_rotation_iterations",
                           gprDimerInitRotationsMax);
    gprDimerRelaxRotationsMax =
        (int)ini.GetValueL("GPR Dimer", "max_relaxation_rotation_iterations",
                           gprDimerRelaxRotationsMax);
    gprDimerDivisorTdimerGP = (int)ini.GetValueL("GPR Dimer", "divisor_t_dimer",
                                                 gprDimerDivisorTdimerGP);
    gprDimerMaxOuterIterations = (int)ini.GetValueL(
        "GPR Dimer", "max_outer_iterations", gprDimerMaxOuterIterations);
    gprDimerMaxInnerIterations = (int)ini.GetValueL(
        "GPR Dimer", "max_inner_iterations", gprDimerMaxInnerIterations);
    gprDimerMidpointMaxDisp = ini.GetValueF(
        "GPR Dimer", "max_midpoint_displacement", gprDimerMidpointMaxDisp);
    gprDimerRotOptMethod =
        ini.GetValue("GPR Dimer", "rotation_opt_method", gprDimerRotOptMethod);
    gprDimerTransOptMethod = ini.GetValue("GPR Dimer", "translation_opt_method",
                                          gprDimerTransOptMethod);
    gprActiveRadius =
        ini.GetValueF("GPR Dimer", "active_radius", gprActiveRadius);
    gprDimerSep = ini.GetValueF("GPR Dimer", "dimer_separation", gprDimerSep);
    gprd_trans_options.step_length = ini.GetValueF(
        "GPR Dimer", "convex_region_step_size", gprd_trans_options.step_length);
    gprd_trans_options.max_step_length = ini.GetValueF(
        "GPR Dimer", "max_step_size", gprd_trans_options.max_step_length);
    gprd_trans_options.rotrem_thresh =
        ini.GetValueF("GPR Dimer", "rotation_removal_projection_threshold",
                      gprd_trans_options.rotrem_thresh);
    gprDimerRatioAtLimit =
        ini.GetValueF("GPR Dimer", "ratio_at_limit", gprDimerRatioAtLimit);
    gprDimerInitRotGP =
        ini.GetValueB("GPR Dimer", "nogp_initial_rotations", gprDimerInitRotGP);
    gprDimerManyIterations = ini.GetValueB("GPR Dimer", "has_many_iterations",
                                           gprDimerManyIterations);
    // GPR Params
    gprDimerSigma2 = ini.GetValueF("GPR Dimer", "gpr_variance", gprDimerSigma2);
    gprDimerJitterSigma2 =
        ini.GetValueF("GPR Dimer", "gpr_jitter_variance", gprDimerJitterSigma2);
    gprDimerNoiseSigma2 =
        ini.GetValueF("GPR Dimer", "gpr_noise_variance", gprDimerNoiseSigma2);
    gprDimerPriorMu = ini.GetValueF("GPR Dimer", "prior_mean", gprDimerPriorMu);
    gprDimerPriorSigma2 =
        ini.GetValueF("GPR Dimer", "prior_variance", gprDimerPriorSigma2);
    gprDimerPriorNu =
        ini.GetValueF("GPR Dimer", "prior_degrees_of_freedom", gprDimerPriorNu);
    // GPR Optimization Parameters
    gpr_hypopt_options.hopt_method =
        ini.GetValue("GPR Dimer", "hyperparameter_opt_method",
                     gpr_hypopt_options.hopt_method);
    gpr_hypopt_options.check_derivative = ini.GetValueB(
        "GPR Dimer", "check_derivatives", gpr_hypopt_options.check_derivative);
    gpr_hypopt_options.tol_func =
        ini.GetValueF("GPR Dimer", "opt_tol_func", gpr_hypopt_options.tol_func);
    gpr_hypopt_options.tol_sol =
        ini.GetValueF("GPR Dimer", "opt_tol_sol", gpr_hypopt_options.tol_sol);
    gpr_hypopt_options.max_iter = ini.GetValueI(
        "GPR Dimer", "opt_max_iterations", gpr_hypopt_options.max_iter);
    gpr_hypopt_options.scg.lambda_limit = ini.GetValueF(
        "GPR Dimer", "scg_lambda_limit", gpr_hypopt_options.scg.lambda_limit);
    gpr_hypopt_options.scg.lambda = ini.GetValueF(
        "GPR Dimer", "scg_lambda_init", gpr_hypopt_options.scg.lambda);
    gpr_hypopt_options.adam.lr = ini.GetValueF(
        "GPR Dimer", "adam_learning_rate", gpr_hypopt_options.adam.lr);
    gpr_hypopt_options.adam.lrd = ini.GetValueF(
        "GPR Dimer", "adam_learning_rate_decay", gpr_hypopt_options.adam.lrd);
    gpr_hypopt_options.adam.b1 =
        ini.GetValueF("GPR Dimer", "adam_beta1", gpr_hypopt_options.adam.b1);
    gpr_hypopt_options.adam.b2 =
        ini.GetValueF("GPR Dimer", "adam_beta2", gpr_hypopt_options.adam.b2);
    gpr_hypopt_options.adam.eps =
        ini.GetValueF("GPR Dimer", "adam_epsilon", gpr_hypopt_options.adam.eps);
    gpr_hypopt_options.adam.weight_decay = ini.GetValueF(
        "GPR Dimer", "adam_weight_decay", gpr_hypopt_options.adam.weight_decay);
    gpr_hypopt_options.adam.amsgrad = ini.GetValueB(
        "GPR Dimer", "adam_amsgrad", gpr_hypopt_options.adam.amsgrad);
    early_stopping_options.dist_metrics = ini.GetValue(
        "GPR Dimer", "es_dist_metric", early_stopping_options.dist_metrics);
    early_stopping_options.threshold = ini.GetValueF(
        "GPR Dimer", "es_threshold", early_stopping_options.threshold);
    fps_options.metric =
        ini.GetValue("GPR Dimer", "fps_metric", fps_options.metric);
    fps_options.history =
        ini.GetValueF("GPR Dimer", "fps_history", fps_options.history);
    // GPR Debugging Parameters
    gprReportLevel =
        (int)ini.GetValueL("GPR Dimer", "report_level", gprReportLevel);
    gprDebugLevel =
        (int)ini.GetValueL("GPR Dimer", "debug_level", gprDebugLevel);
    gprDebugOutDir =
        ini.GetValue("GPR Dimer", "debug_output_directory", gprDebugOutDir);
    ;
    gprDebugPosFile =
        ini.GetValue("GPR Dimer", "debug_position_basename", gprDebugPosFile);
    ;
    gprDebugEnergyFile =
        ini.GetValue("GPR Dimer", "debug_energy_basename", gprDebugEnergyFile);
    ;
    gprDebugGradFile =
        ini.GetValue("GPR Dimer", "debug_gradient_basename", gprDebugGradFile);
    ;
    gprDebugOutExt = ini.GetValue("GPR Dimer", "debug_out_ext", gprDebugOutExt);
    ;
    gprDebugOffsetMidPoint = ini.GetValueF("GPR Dimer", "debug_midpoint_offset",
                                           gprDebugOffsetMidPoint);
    gprDebugDy = ini.GetValueF("GPR Dimer", "debug_y_step", gprDebugDy);
    gprDebugDz = ini.GetValueF("GPR Dimer", "debug_z_step", gprDebugDz);
    // GPR Prune
    gprUsePrune = ini.GetValueB("GPR Dimer", "use_prune", gprUsePrune);
    gprPruneBegin =
        (int)ini.GetValueL("GPR Dimer", "start_prune_at", gprPruneBegin);
    gprPruneNVals =
        (int)ini.GetValueL("GPR Dimer", "nprune_vals", gprPruneNVals);
    gprPruneThreshold =
        ini.GetValueF("GPR Dimer", "prune_threshold", gprPruneThreshold);

    // [Prefactor] //

    prefactorDefaultValue =
        ini.GetValueF("Prefactor", "default_value", prefactorDefaultValue);
    prefactorMaxValue =
        ini.GetValueF("Prefactor", "max_value", prefactorMaxValue);
    prefactorMinValue =
        ini.GetValueF("Prefactor", "min_value", prefactorMinValue);
    prefactorWithinRadius =
        ini.GetValueF("Prefactor", "within_radius", prefactorWithinRadius);
    prefactorMinDisplacement = ini.GetValueF("Prefactor", "min_displacement",
                                             prefactorMinDisplacement);
    prefactorRate = toLowerCase(
        ini.GetValue("Prefactor", "rate_estimation", prefactorRate));
    prefactorConfiguration = toLowerCase(
        ini.GetValue("Prefactor", "configuration", prefactorConfiguration));
    prefactorAllFreeAtoms =
        ini.GetValueB("Prefactor", "all_free_atoms", prefactorAllFreeAtoms);
    prefactorFilterScheme = toLowerCase(
        ini.GetValue("Prefactor", "filter_scheme", prefactorFilterScheme));
    prefactorFilterFraction =
        ini.GetValueF("Prefactor", "filter_fraction", prefactorFilterFraction);

    // [Hessian] //

    hessianAtomList =
        toLowerCase(ini.GetValue("Hessian", "atom_list", hessianAtomList));
    hessianZeroFreqValue =
        ini.GetValueF("Hessian", "zero_freq_value", hessianZeroFreqValue);

    // [Nudged Elastic Band] //

    nebImages = ini.GetValueL("Nudged Elastic Band", "images", nebImages);
    nebSpring = ini.GetValueF("Nudged Elastic Band", "spring", nebSpring);
    nebClimbingImageMethod = ini.GetValueB(
        "Nudged Elastic Band", "climbing_image_method", nebClimbingImageMethod);
    nebClimbingImageConvergedOnly =
        ini.GetValueB("Nudged Elastic Band", "climbing_image_converged_only",
                      nebClimbingImageConvergedOnly);
    nebOldTangent =
        ini.GetValueB("Nudged Elastic Band", "old_tangent", nebOldTangent);
    nebMaxIterations = ini.GetValueL("Nudged Elastic Band", "max_iterations",
                                     optMaxIterations);
    nebDoublyNudged =
        ini.GetValueB("Nudged Elastic Band", "doubly_nudged", nebDoublyNudged);
    nebDoublyNudgedSwitching =
        ini.GetValueB("Nudged Elastic Band", "doubly_nudged_switching",
                      nebDoublyNudgedSwitching);
    nebElasticBand =
        ini.GetValueB("Nudged Elastic Band", "elastic_band", nebElasticBand);
    nebConvergedForce = ini.GetValueF("Nudged Elastic Band", "converged_force",
                                      optConvergedForce);
    nebEnergyWeighted = ini.GetValueB("Nudged Elastic Band", "energy_weighted",
                                      nebEnergyWeighted);
    nebKSPMin = ini.GetValueF("Nudged Elastic Band", "ew_ksp_min", nebKSPMin);
    nebKSPMax = ini.GetValueF("Nudged Elastic Band", "ew_ksp_max", nebKSPMax);
    nebIpath = ini.GetValue("Nudged Elastic Band", "initial_path_in", nebIpath);
    nebMinimEP =
        ini.GetValueB("Nudged Elastic Band", "minimize_endpoints", nebMinimEP);
    nebciAfter = ini.GetValueF("Nudged Elastic Band", "ci_after", nebciAfter);
    nebciWithMMF = ini.GetValueB("Nudged Elastic Band", "ci_mmf", nebciWithMMF);
    nebciMMFAfter =
        ini.GetValueF("Nudged Elastic Band", "ci_mmf_after", nebciMMFAfter);
    nebciMMFnSteps =
        ini.GetValueL("Nudged Elastic Band", "ci_mmf_nsteps", nebciMMFnSteps);

    // [Dynamics] //

    mdTimeStepInput = ini.GetValueF("Dynamics", "time_step", mdTimeStepInput);
    mdTimeStep = mdTimeStepInput / timeUnit;
    mdTimeInput = ini.GetValueF("Dynamics", "time", mdTimeInput);
    mdTime = mdTimeInput / timeUnit;
    mdSteps = long(floor(mdTime / mdTimeStep + 0.5));
    thermostat =
        toLowerCase(ini.GetValue("Dynamics", "thermostat", "andersen"));
    thermoAndersenAlpha =
        ini.GetValueF("Dynamics", "andersen_alpha", thermoAndersenAlpha);
    thermoAndersenTcolInput = ini.GetValueF(
        "Dynamics", "andersen_collision_period", thermoAndersenTcolInput);
    thermoAndersenTcol = thermoAndersenTcolInput / timeUnit;
    thermoNoseMass = ini.GetValueF("Dynamics", "nose_mass", thermoNoseMass);
    thermoLangevinFrictionInput = ini.GetValueF("Dynamics", "langevin_friction",
                                                thermoLangevinFrictionInput);
    thermoLangevinFriction = thermoLangevinFrictionInput * timeUnit;
    // thermoAtoms =
    // helper_functions::split_string_int(ini.GetValue("Dynamics",
    // "thermo_atoms", ""), ",");

    // [Parallel Replica]

    parrepAutoStop = ini.GetValueB("Parallel Replica", "stop_after_transition",
                                   parrepAutoStop);
    parrepRefineTransition = ini.GetValueB(
        "Parallel Replica", "refine_transition", parrepRefineTransition);
    parrepDephaseLoopStop = ini.GetValueB(
        "Parallel Replica", "dephase_loop_stop", parrepDephaseLoopStop);
    parrepDephaseTimeInput = ini.GetValueF("Parallel Replica", "dephase_time",
                                           parrepDephaseTimeInput);
    parrepDephaseTime = parrepDephaseTimeInput / timeUnit;
    parrepDephaseLoopMax = ini.GetValueL("Parallel Replica", "dephase_loop_max",
                                         parrepDephaseLoopMax);
    parrepStateCheckIntervalInput =
        ini.GetValueF("Parallel Replica", "state_check_interval",
                      parrepStateCheckIntervalInput);
    parrepStateCheckInterval = parrepStateCheckIntervalInput / timeUnit;
    parrepRecordIntervalInput =
        ini.GetValueF("Parallel Replica", "state_save_interval",
                      0.1 * parrepStateCheckIntervalInput);
    parrepRecordInterval = parrepRecordIntervalInput / timeUnit;
    parrepCorrTimeInput = ini.GetValueF(
        "Parallel Replica", "post_transition_time", parrepCorrTimeInput);
    parrepCorrTime = parrepCorrTimeInput / timeUnit;

    //[Temperature Accelerated Dynamics] //

    tadLowT = ini.GetValueF("TAD", "low_temperature", tadLowT);
    tadMinPrefactor = ini.GetValueF("TAD", "min_prefactor", tadMinPrefactor);
    tadConfidence = ini.GetValueF("TAD", "confidence", tadConfidence);

    // [Replica Exchange] //

    repexcTemperatureDistribution =
        toLowerCase(ini.GetValue("Replica Exchange", "temperature_distribution",
                                 repexcTemperatureDistribution));
    repexcReplicas =
        ini.GetValueL("Replica Exchange", "replicas", repexcReplicas);
    repexcExchangeTrials = ini.GetValueL("Replica Exchange", "exchange_trials",
                                         repexcExchangeTrials);
    repexcSamplingTimeInput = ini.GetValueF("Replica Exchange", "sampling_time",
                                            repexcSamplingTimeInput);
    repexcSamplingTime = repexcSamplingTimeInput / timeUnit;
    repexcTemperatureLow =
        ini.GetValueF("Replica Exchange", "temperature_low", temperature);
    repexcTemperatureHigh = ini.GetValueF(
        "Replica Exchange", "temperature_high", repexcTemperatureHigh);
    repexcExchangePeriodInput = ini.GetValueF(
        "Replica Exchange", "exchange_period", repexcExchangePeriodInput);
    repexcExchangePeriod = repexcExchangePeriodInput / timeUnit;

    // [Hyperdynamics] //

    bondBoostRMDTimeInput =
        ini.GetValueF("Hyperdynamics", "bb_rmd_time", bondBoostRMDTimeInput);
    bondBoostRMDTime = bondBoostRMDTimeInput / timeUnit;
    bondBoostBALS = toLowerCase(
        ini.GetValue("Hyperdynamics", "bb_boost_atomlist", bondBoostBALS));
    bondBoostDVMAX = ini.GetValueF("Hyperdynamics", "bb_dvmax", bondBoostDVMAX);
    bondBoostQRR =
        ini.GetValueF("Hyperdynamics", "bb_stretch_threshold", bondBoostQRR);
    bondBoostPRR =
        ini.GetValueF("Hyperdynamics", "bb_ds_curvature", bondBoostPRR);
    bondBoostQcut = ini.GetValueF("Hyperdynamics", "bb_rcut", bondBoostQcut);
    biasPotential = toLowerCase(
        ini.GetValue("Hyperdynamics", "bias_potential", biasPotential));

    // [Saddle Search] //

    saddleMethod =
        toLowerCase(ini.GetValue("Saddle Search", "method", saddleMethod));
    saddleMinmodeMethod = toLowerCase(
        ini.GetValue("Saddle Search", "min_mode_method", saddleMinmodeMethod));
    saddleDisplaceMagnitude = ini.GetValueF(
        "Saddle Search", "displace_magnitude", saddleDisplaceMagnitude);
    saddleDisplaceRadius =
        ini.GetValueF("Saddle Search", "displace_radius", saddleDisplaceRadius);
    saddleMaxEnergy =
        ini.GetValueF("Saddle Search", "max_energy", saddleMaxEnergy);
    saddleMaxIterations =
        ini.GetValueL("Saddle Search", "max_iterations", optMaxIterations);
    saddleNonnegativeDisplacementAbort =
        ini.GetValueB("Saddle Search", "nonnegative_displacement_abort",
                      saddleNonnegativeDisplacementAbort);
    saddleMaxSingleDisplace = ini.GetValueF(
        "Saddle Search", "max_single_displace", saddleMaxSingleDisplace);
    // must be loaded after optConvergedForce
    saddleConvergedForce =
        ini.GetValueF("Saddle Search", "converged_force", optConvergedForce);
    saddlePerpForceRatio = ini.GetValueF("Saddle Search", "perp_force_ratio",
                                         saddlePerpForceRatio); // undocumented
    saddleDisplaceType = toLowerCase(ini.GetValue(
        "Saddle Search", "client_displace_type", EpiCenters::DISP_LOAD));
    saddleNonlocalCountAbort =
        ini.GetValueL("Saddle Search", "nonlocal_count_abort",
                      saddleNonlocalCountAbort); // undocumented
    saddleNonlocalDistanceAbort =
        ini.GetValueF("Saddle Search", "nonlocal_distance_abort",
                      saddleNonlocalDistanceAbort); // undocumented
    if (saddleDisplaceType != EpiCenters::DISP_NOT_FCC_OR_HCP &&
        saddleDisplaceType != EpiCenters::DISP_MIN_COORDINATED &&
        saddleDisplaceType != EpiCenters::DISP_LAST_ATOM &&
        saddleDisplaceType != EpiCenters::DISP_RANDOM) {
      saddleDisplaceType = EpiCenters::DISP_LOAD;
    }
    saddleConfinePositive = ini.GetValueB("Saddle Search", "confine_positive",
                                          saddleConfinePositive);
    if (saddleConfinePositive) {
      saddleBowlBreakout =
          ini.GetValueB("Saddle Search", "bowl_breakout", saddleBowlBreakout);
      saddleBowlActive =
          ini.GetValueL("Saddle Search", "bowl_active_atoms", saddleBowlActive);
      saddleConfinePositiveMinForce =
          ini.GetValueF("Saddle Search", "confine_positive_min_move",
                        saddleConfinePositiveMinForce);
      saddleConfinePositiveScaleRatio =
          ini.GetValueF("Saddle Search", "confine_positive_scale_ratio",
                        saddleConfinePositiveScaleRatio);
      saddleConfinePositiveBoost =
          ini.GetValueF("Saddle Search", "confine_positive_boost",
                        saddleConfinePositiveBoost);
      saddleConfinePositiveMinActive =
          ini.GetValueL("Saddle Search", "confine_positive_min_active",
                        saddleConfinePositiveMinActive);
    }
    saddleDynamicsTemperature = temperature;
    saddleDynamicsTemperature = ini.GetValueF(
        "Saddle Search", "dynamics_temperature", saddleDynamicsTemperature);
    saddleDynamicsStateCheckIntervalInput =
        ini.GetValueF("Saddle Search", "dynamics_state_check_interval",
                      saddleDynamicsStateCheckIntervalInput);
    saddleDynamicsStateCheckInterval =
        saddleDynamicsStateCheckIntervalInput / timeUnit;
    saddleDynamicsRecordIntervalInput =
        ini.GetValueF("Saddle Search", "dynamics_record_interval",
                      saddleDynamicsRecordIntervalInput);
    saddleDynamicsRecordInterval = saddleDynamicsRecordIntervalInput / timeUnit;
    saddleDynamicsLinearInterpolation =
        ini.GetValueB("Saddle Search", "dynamics_linear_interpolation",
                      saddleDynamicsLinearInterpolation);
    saddleRemoveRotation =
        ini.GetValueB("Saddle Search", "remove_rotation", saddleRemoveRotation);
    saddleDynamicsMaxInitCurvature =
        ini.GetValueF("Saddle Search", "dynamics_max_init_curvature",
                      saddleDynamicsMaxInitCurvature);
    saddleZeroModeAbortCurvature =
        ini.GetValueF("Saddle Search", "zero_mode_abort_curvature",
                      saddleZeroModeAbortCurvature);

    // [Basin Hopping] //

    basinHoppingDisplacement = ini.GetValueF("Basin Hopping", "displacement",
                                             basinHoppingDisplacement);
    basinHoppingPushApartDistance = ini.GetValueF(
        "Basin Hopping", "push_apart_distance", basinHoppingPushApartDistance);
    basinHoppingInitialRandomStructureProbability =
        ini.GetValueF("Basin Hopping", "initial_random_structure_probability",
                      basinHoppingInitialRandomStructureProbability);
    basinHoppingSteps =
        ini.GetValueL("Basin Hopping", "steps", basinHoppingSteps);
    basinHoppingQuenchingSteps = ini.GetValueL(
        "Basin Hopping", "quenching_steps", basinHoppingQuenchingSteps);
    basinHoppingSingleAtomDisplace =
        ini.GetValueB("Basin Hopping", "single_atom_displace",
                      basinHoppingSingleAtomDisplace);
    basinHoppingSignificantStructure =
        ini.GetValueB("Basin Hopping", "significant_structure",
                      basinHoppingSignificantStructure);
    basinHoppingDisplacementAlgorithm =
        toLowerCase(ini.GetValue("Basin Hopping", "displacement_algorithm",
                                 basinHoppingDisplacementAlgorithm));
    basinHoppingDisplacementDistribution =
        toLowerCase(ini.GetValue("Basin Hopping", "displacement_distribution",
                                 basinHoppingDisplacementDistribution));
    basinHoppingSwapProbability = ini.GetValueF(
        "Basin Hopping", "swap_probability", basinHoppingSwapProbability);
    basinHoppingJumpMax =
        ini.GetValueL("Basin Hopping", "jump_max", basinHoppingJumpMax);
    basinHoppingJumpSteps =
        ini.GetValueL("Basin Hopping", "jump_steps", basinHoppingJumpSteps);
    basinHoppingAdjustDisplacement = ini.GetValueB(
        "Basin Hopping", "adjust_displacement", basinHoppingAdjustDisplacement);
    basinHoppingAdjustPeriod = ini.GetValueL("Basin Hopping", "adjust_period",
                                             basinHoppingAdjustPeriod);
    basinHoppingAdjustFraction = ini.GetValueF(
        "Basin Hopping", "adjust_fraction", basinHoppingAdjustFraction);
    basinHoppingTargetRatio =
        ini.GetValueF("Basin Hopping", "target_ratio", basinHoppingTargetRatio);
    basinHoppingWriteUnique =
        ini.GetValueB("Basin Hopping", "write_unique", basinHoppingWriteUnique);
    basinHoppingStopEnergy =
        ini.GetValueF("Basin Hopping", "stop_energy", basinHoppingStopEnergy);

    // [Global Optimization] //

    globalOptimizationMoveMethod = toLowerCase(ini.GetValue(
        "Global Optimization", "move_method", globalOptimizationMoveMethod));
    globalOptimizationDecisionMethod =
        toLowerCase(ini.GetValue("Global Optimization", "decision_method",
                                 globalOptimizationDecisionMethod));
    globalOptimizationSteps =
        ini.GetValueL("Global Optimization", "steps", globalOptimizationSteps);
    globalOptimizationBeta =
        ini.GetValueF("Global Optimization", "beta", globalOptimizationBeta);
    globalOptimizationAlpha =
        ini.GetValueF("Global Optimization", "alpha", globalOptimizationAlpha);
    globalOptimizationMdmin =
        ini.GetValueL("Global Optimization", "mdmin", globalOptimizationMdmin);
    globalOptimizationTargetEnergy = ini.GetValueF(
        "Global Optimization", "target_energy", globalOptimizationTargetEnergy);

    // [BGSD] //

    alpha = ini.GetValueF("BGSD", "alpha", alpha);
    beta = ini.GetValueF("BGSD", "beta", beta);
    gradientfinitedifference = ini.GetValueF("BGSD", "gradientfinitedifference",
                                             gradientfinitedifference);
    grad2energyconvergence =
        ini.GetValueF("BGSD", "grad2energyconvergence", grad2energyconvergence);
    grad2forceconvergence =
        ini.GetValueF("BGSD", "grad2forceconvergence", grad2forceconvergence);

    // [Monte Carlo] //

    monteCarloStepSize =
        ini.GetValueF("Monte Carlo", "step_size", monteCarloStepSize);
    monteCarloSteps = ini.GetValueI("Monte Carlo", "steps", monteCarloSteps);

    // Sanity Checks
    if (parrepStateCheckInterval > mdTime &&
        magic_enum::enum_name<JobType>(job) == "parallel_replica") {
      SPDLOG_ERROR("[Parallel Replica] state_check_interval must be <= time");
      error = 1;
    }

    if (saddleDynamicsRecordIntervalInput >
        saddleDynamicsStateCheckIntervalInput) {
      SPDLOG_ERROR("[Saddle Search] dynamics_record_interval must be <= "
                   "dynamics_state_check_interval");
      error = 1;
    }

    if (potential == PotType::AMS || potential == PotType::AMS_IO) {
      if (forcefield.empty() && model.empty() && xc.empty()) {
        SPDLOG_ERROR("[AMS] Must provide atleast forcefield or model or xc");
        error = 1;
      }

      if (!forcefield.empty() && !model.empty() && !xc.empty()) {
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

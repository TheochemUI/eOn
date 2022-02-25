#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstring>
#include <cfloat>
#include <ctype.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <string_view>
#include <fstream>

#include "StaticData.h"

using namespace std;
using namespace std::string_literals; // For ""s

#ifdef EONMPI
#    include "mpi.h"
#endif

/** Contains all runtime parameters and results. No functionality just bookkeeping.*/
class Parameters {

public:
    Parameters()
        : // clang-format off
        // [ Main ] //
        kB { 8.6173324e-5 /* eV/K */ },
        timeUnit { 10.1805055 /* fs */ },
        job { JobStrings::PROCESS_SEARCH /* default job */ },
        randomSeed { 1995 /* randomized default */ },
        temperature { 300.0 /* room temperature */ },
        checkpoint { false /* i/o is slow */ },
        quiet { false /* verbosity */ },
        writeLog { true /* logging is good */ },
        iniFilename { "config.ini"s /*  */ },
        conFilename { "pos.con"s /*  */ },
        finiteDifference { 0.01 /*  */ },
        maxForceCalls { 0 /*  */ },
        removeNetForce { true /*  */ },
        // [Prefactor] //
        prefactorDefaultValue { 0.0 /* default prefactor; calculate explicitly if zero */ },
        prefactorMaxValue { 1e+21 /* max prefactor allowed */ },
        prefactorMinValue { 1e+9 /* min prefactor allowed */ },
        prefactorWithinRadius { 3.3 /* atoms within this radius of the displaced atoms are put in the Hessian, unless filterMode is fraction */ },
        prefactorMinDisplacement { 0.25 /* atoms with displacement between min1 or min2 and the saddle point are put in the Hessian method to estimate prefactor */ },
        prefactorRate { PrefactorStrings::RATE_HTST /*  */ },
        prefactorConfiguration { PrefactorJobStrings::PREFACTOR_REACTANT /* configuration for which the frequencies should be determined */ },
        prefactorAllFreeAtoms { false /* use all free atom when determining the prefactor */ },
        prefactorFilterScheme { PrefactorStrings::FILTER_FRACTION /* "cutoff" or "fraction", which use prefactorMinDisplacement or prefactorFilterFraction, respectively. */ },
        prefactorFilterFraction { 0.90 /* Include atoms whose summed motion comprise more than prefactorFilterFraction in the prefactor calculation. Prioritizes atoms that move more. */ },
        // [Potential] //
        potential { PotentialStrings::POT_LJ /*  */},
        MPIPollPeriod { 0.25 /* second */},
        MPIPotentialRank { -1 /*  */ },
        LogPotential { false /*  */ },
        LAMMPSLogging { false /*  */ },
        LAMMPSThreads { 0 /*  */ },
        EMTRasmussen { false /*  */ },
        extPotPath { "./ext_pot"s /*  */ },
        // [AMS] //
        engine { ""s /* MOPAC, ADF, BAND, REAXFF, FORCEFIELD */ },
        forcefield { ""s /* OPt.ff etc. (REAXFF) */ },
        model { ""s /* PM7 PM3 ->  Model hamiltonian (MOPAC) */ },
        xc { ""s /* Exchange (BAND, ADF) */ },
        basis { ""s /* Basis (BAND, ADF) */ },
        resources { ""s /* DFTB */ },
        // [AMS_ENV] //
        // Horrid little section to mimic amsrc.sh
        // Assumes the entire thing is going to be set
        amshome { ""s /* "/some/path/to/amshome/"s */ },
        scm_tmpdir { ""s /* "/tmp"s */ },
        scm_pythondir { ""s /* "/.scm/python" */ },
        amsbin { ""s /* amshome.append("/bin") */ },
        scmlicense { ""s /* amshome.append("license.txt") */ },
        amsresources { ""s /* amshome.append("/atomicdata") */ },
        // [Structure Comparison] //
        distanceDifference { 0.1 /* The distance criterion for comparing geometries */ },
        neighborCutoff { 3.3 /* Radius used in the local atomic structure analysis */ },
        checkRotation { false /*  */ },
        indistinguishableAtoms { true /*  */ },
        energyDifference { 0.01 /*  */ },
        removeTranslation { true /*  */ },
        // [Debug] //
        writeMovies { false /*  */ },
        writeMoviesInterval { 1 /*  */ },
        // [Saddle Search] //
        saddleDisplaceType { EpiCentersStrings::DISP_LOAD /* displacement type to use */ },
        saddleMethod { "min_mode"s /*  */ },
        saddleMinmodeMethod { LowestEigenmodeStrings::MINMODE_DIMER /* algorithm to be used for lowest eigenmode determination */ },
        saddleMaxEnergy { 20.0 /* energy above product state that will cause termination of the saddle point search */ },
        saddleMaxIterations { 1000 /*  */ },
        saddleDisplaceRadius { 4.0 /* atoms within this radius of the displacement atoms are also displaced */ },
        saddleDisplaceMagnitude { 0.1 /* norm of the displacement vector */ },
        saddleMaxSingleDisplace { 10.0 /*  number of displacements to reach a convex region;  if 0, a search is started after the displacement max iterations for saddle point searches and minimization */ },
        /* ^ maximum value of displacement in x, y and z direction for atoms being displaced */
//    saddleConvergedForce = optConvergedForce; default value is set after the value of optConvergedForce is loaded
//    force convergence criterion required for a saddle point search ^
        saddleNonnegativeDisplacementAbort { false /* abort the saddle search if the displacement does not have a negative mode */ },
        saddleNonlocalCountAbort { 0 /* abort the search if this many atoms move more than NonlocalDistanceAbort */ },
        saddleNonlocalDistanceAbort { 0.0 /* abort the search if NonlocalCountAbort atoms move more than this distance */ },
        saddleRemoveRotation { false /* remove dominant rotational component when system is translated */ },
        saddlePerpForceRatio { 0.0 /* proportion to keep of the perpendicular force when the lowest eigenvalue is positive */ },
        saddleConfinePositive { false /* undocumented */ },
        saddleBowlBreakout { false /* undocumented */ },
        saddleBowlActive { 20 /* number of active atoms */ },
        saddleConfinePositiveMinForce { 0.5 /* undocumented */ },
        saddleConfinePositiveScaleRatio { 0.9 /* undocumented */ },
        saddleConfinePositiveBoost { 10. /* undocumented */ },
        saddleConfinePositiveMinActive { 30 /* undocumented */ },
        saddleDynamicsTemperature { temperature /*defaults to temperature, temperature for dynamics saddle search method */ },
        saddleDynamicsStateCheckIntervalInput { 100.0 /*fs */ },
        saddleDynamicsStateCheckInterval { 5.0 /* fs, how often to minimize */},
        saddleDynamicsRecordIntervalInput { 10.0 /*fs */ },
        saddleDynamicsLinearInterpolation { true /*  */ },
        saddleDynamicsMaxInitCurvature { 0.0 /* eV/Ang^2 */ },
        saddleZeroModeAbortCurvature { 0.0 /* eV/Ang^2 */ },
        // [Optimizers] //
        optMethod { "cg"s /*  */ },
        optConvergenceMetric { "norm"s /* norm, max_atom, max_component */ },
        optMaxIterations { 1000 /* maximum iterations for saddle point searches and minimization */ },
        optConvergedForce { 0.01 /* force convergence criterion required for an optimization */ },
        optMaxMove { 0.2 /* maximum displacement vector for a step during optimization */ },
        optTimeStepInput { 1.0 /*  */ },
        optTimeStep { optTimeStepInput /* time step size used in quickmin */},
        optMaxTimeStepInput { 2.5 /* maximum time step for FIRE. */ },
        optLBFGSMemory { 20 /* number of previous forces to keep in the bfgs memory */ },
        optLBFGSInverseCurvature { 0.01 /* assumes stiffest curvature at minimum is 100 eV/A^2 */ },
        optLBFGSAutoScale { true /*  */ },
        optLBFGSAngleReset { true /*  */ },
        optLBFGSDistanceReset { true /*  */ },
        optQMSteepestDecent { false /* if set the velocity will always be set to zero in quickmin */ },
        optCGNoOvershooting { false /* if set it is ensured that the approximate line search in conjugate gradients never overshoot the minimum along the search line*/ },
        optCGKnockOutMaxMove { false /* if set the old search direction is nullified when steps larger than the optMaxMove are conducted */ },
        optCGLineConverged { 0.1 /* convergence criteria for line search, ratio between force component along search line and the orthogonal part */ },
        optCGLineSearch { false /* if set full line search is conducted */ },
        optCGMaxIterBeforeReset { 0 /* max nr of cg steps before reset, if 0 no resetting is done */ },
        optCGLineSearchMaxIter { 10 /* maximal nr of iterations during line search */ },
        optSDAlpha { 0.1 /*  */ },
        optSDTwoPoint { false /*  */ },
        // [Process Search] //
        processSearchMinimizeFirst { true /*  */ },
        processSearchMinimizationOffset { optMaxMove /* how far from the saddle to displace the minimization images */ },
        // [Dimer] //
        dimerRotationAngle { 0.005 /* finite difference rotation angle */ },
        dimerImproved { true /* turn on the improved dimer method */ },
        dimerConvergedAngle { 5.0 /* degrees, stop rotating when angle drops below this value */ },
        dimerOptMethod { ImprovedDimerStrings::OPT_CG /* method to determine the next rotation direction */ },
        dimerTorqueMin { 0.1 /* old dimer */ },
        dimerTorqueMax { 1.0 /* old dimer */ },
        dimerRotationsMin { 1 /* old dimer */ },
        dimerRotationsMax { 10 /* old dimer and new dimer */ },
        dimerMaxIterations { 1000 /* maximum number of rotation iterations */ },
        dimerRemoveRotation { false /* remove dominant rotational component when estimating the eigenmode */ },
        // [Lanczos] //
        lanczosTolerance { 0.01 /* difference between the lowest eigenvalues of two successive iterations */ },
        lanczosMaxIterations { 20 /* maximum number of iterations */ },
        lanczosQuitEarly { true /*  */ },
        // [GPR Dimer] //
        gprDimerRotationAngle { 0.005 /* finite difference rotation angle */ },
        gprDimerConvergedAngle { 0.08 /* {T_anglerot_init} stop rotating when angle drops below this value */ },
        gprDimerRelaxConvAngle { 0.001 /* {T_anglerot_gp} stop rotating when angle drops below this value during relaxation */ },
        gprDimerInitRotationsMax { 6 /* {num_iter_initrot} should be DoF */ },
        gprDimerRelaxRotationsMax { 10 /* {num_iter_rot_gp} */ },
        gprDimerDivisorTdimerGP { 10 /* {divisor_T_dimer_gp} */ },
        gprDimerMaxOuterIterations { 300 /* {num_bigiter} maximum number of outer iterations or new sets of observations */ },
        gprDimerMaxInnerIterations { 1000 /* {num_iter} maximum number of steps during the relaxation phase */ },
        gprDimerMidpointMaxDisp { 0.5 /* {disp_max} */ },
        gprDimerRotOptMethod { "lbfgs"s /* {method_rot} method to determine the next rotation direction */ },
        gprDimerTransOptMethod { "lbfgs"s /* {method_trans} method to determine the next rotation direction */ },
        gprActiveRadius { 5.0 /* {actidst_fro} activation radius for inclusion in covariance matrix */ },
        gprDimerSep { 0.01 /* {dimer_sep} distance from the middle point of the dimer to the two images */ },
        gprDimerConvStep { 0.1 /* {param_trans[0]} step length for convex regions */ },
        gprDimerMaxStep { 0.1 /* {param_trans[1]} maximum step length */ },
        gprForceThreshold { 0.01 /* {T_dimer} maximum component of the force acting on the middle point of the dimer (stops when accurate force components are below this value) */ },
        gprDimerRatioAtLimit { 0.66667 /* {ratio_at_limit} defines the limit for the ratio of inter-atomic distances between the image and its "nearest observed data point" */ },
        gprDimerInitRotGP { false /* {initrot_nogp} initial rotations without GP (1) or with GP (0) */ },
        gprDimerInitTransGP { false /* {inittrans_nogp} initial translations without GP (1) or with GP (0) */ },
        gprDimerManyIterations { true /* {islarge_num_iter} indicates of the number of iterations is larger than required for dimer convergence on the accurate energy surface (1), otherwise (0) the relaxation phase is continued from the current dimer if the maximum iterations in the relaxation phase is reached */ },
        // GPR Params
        gprDimerHyperOptMethod { "scg"s /* {optimization_alg} method to optimize hyperparameters */ },
        gprDimerSigma2 { 1e-8 /* {gp_sigma2} GPR variance {gp_sigma2} */ },
        gprDimerJitterSigma2 { 0 /* {jitter_sigma2} GPR jitter variance */ },
        gprDimerNoiseSigma2 { 1e-8 /* {sigma2} noise Variance */ },
        gprDimerPriorMu { 0 /* {prior_mu} prior mean */ },
        gprDimerPriorSigma2 { 1 /* {prior_s2} prior variance */ },
        gprDimerPriorNu { 20 /* {prior_nu} prior degrees of freedom */ },
        // GPR Optimization Parameters
        gprOptCheckDerivatives { false /* {check_derivative} */ },
        gprOptMaxIterations { 400 /* {max_iter} */ },
        gprOptTolFunc { 1e-4 /* {tolerance_func} */ },
        gprOptTolSol { 1e-4 /* {tolerance_sol} */ },
        gprOptLambdaLimit { 1e17 /* {lambda_limit} */ },
        gprOptLambdaInit { 10.0 /* {lambda} */ },
        // GPR Prune Parameters
        gprUsePrune { false /* {use_prune} */ },
        gprPruneBegin { 8 /* {start_prune_at} */ },
        gprPruneNVals { 3 /* {nprune_vals} */ },
        gprPruneThreshold { 0.5 /* {prune_threshold} */ },
        // GPR Debugging Parameters
        gprReportLevel { 1 /* {report_level} */ },
        gprDebugLevel { 2 /* {debug_level} */ },
        gprDebugOutDir { "output"s /* {debug_output_dir} */ },
        gprDebugPosFile { "position"s /* {debug_output_file_R} */ },
        gprDebugEnergyFile { "energy"s /* {debug_output_file_E} */ },
        gprDebugGradFile { "gradient"s /* {debug_output_file_G} */ },
        gprDebugOutExt { "dat"s /* {debug_output_file_extension} */ },
        gprDebugOffsetMidPoint { 3. /* {debug_offset_from_mid_point} */ },
        gprDebugDy { 0.1 /* {debug_dy} */ },
        gprDebugDz { 0.1 /* {debug_dz} */ },
        // [Hessian] //
        hessianAtomList { std::string{"All"s} /*  */ },
        hessianZeroFreqValue { 1e-6 /*  */ },
        // [Nudged Elastic Band] //
        nebImages { 5 /*  */ },
        nebSpring { 5.0 /*  */ },
        nebClimbingImageMethod { true /*  */ },
        nebClimbingImageConvergedOnly { true /*  */ },
        nebOldTangent { false /*  */ },
        nebMaxIterations { 1000 /*  */ },
        nebDoublyNudged { false /*  */ },
        nebDoublyNudgedSwitching { false /*  */ },
        nebElasticBand { false /*  */ },
        nebConvergedForce { optConvergedForce /* force convergence criterion required for an optimization */ },
        // [Dynamics] //
        mdTimeStepInput { 1.0 /*  */ },
        mdTimeInput { 1000.0 /*  */ },
        // [Thermostat] //
        thermostat { DynamicsStrings::NONE /*  */ },
        thermoAndersenAlpha { 1.0 /* collision strength */ },
        thermoAndersenTcolInput { 100.0 /* collision frequency in unit of fs */ },
        thermoNoseMass { 1.0 /*  */},
        thermoLangevinFrictionInput { 0.01 /*  */},
        // [Parallel Replica] //
        parrepRefineTransition { true /*  */ },
        parrepAutoStop { false /*  */ },
        parrepDephaseLoopStop { false /*  */ },
        parrepDephaseTimeInput { 1000.0 /*  */ },
        parrepDephaseLoopMax { 5 /*  */ },
        parrepStateCheckIntervalInput { 1000.0 /*  */ },
        parrepRecordIntervalInput { 50.0 /*  */ },
        parrepCorrTimeInput { 1000.0 /*  */ },
        // [Temperature Accelerated Dynamics] //
        tadLowT { 300.0 /*  */ },
        tadMinPrefactor { 0.001 /* in unit of fs-1 */ },
        tadConfidence { 0.001 /*  */ },
        // [Replica Exchange] //
        repexcTemperatureDistribution { "exponential"s /*  */ },
        repexcReplicas { 10 /*  */ },
        repexcExchangeTrials { repexcReplicas /*  */ },
        repexcSamplingTimeInput { 1000.0 /*  */ },
        repexcTemperatureLow { 0.0 /*  */ },
        repexcTemperatureHigh { 0.0 /*  */ },
        repexcExchangePeriod { 100.0 /*  */ },
        // [Hyperdynamics] //
        biasPotential { Hyperdynamics::NONE /*  */ },
        bondBoostBALS { std::string{"ALL"s} /* boosted atom list string */ },
        bondBoostDVMAX { 0.0 /*  */ },
        bondBoostQRR { 0.2 /* can not be set to 0 */ },
        bondBoostPRR { 0.95 /*  */ },
        bondBoostQcut { 3.0 /*  */ },
        bondBoostRMDTimeInput { 100.0 /*  */ },
        // [Basin Hopping] //
        basinHoppingDisplacement { 0.5 /*  */ },
        basinHoppingPushApartDistance { 0.4 /*  */ },
        basinHoppingInitialRandomStructureProbability { 0.0 /*  */ },
        basinHoppingSteps { 10000 /*  */ },
        basinHoppingQuenchingSteps { 0 /*  */ },
        basinHoppingSingleAtomDisplace { false /*  */ },
        basinHoppingSignificantStructure { true /*  */ },
        basinHoppingDisplacementAlgorithm { "standard"s /*  */ },
        basinHoppingDisplacementDistribution { "uniform"s /*  */ },
        basinHoppingSwapProbability { 0.0 /*  */ },
        basinHoppingJumpMax { 10 /*  */ },
        basinHoppingJumpSteps { 0 /*  */ },
        basinHoppingAdjustDisplacement { true /*  */ },
        basinHoppingAdjustPeriod { 10 /*  */ },
        basinHoppingAdjustFraction { 0.05 /*  */ },
        basinHoppingTargetRatio { 0.5 /*  */ },
        basinHoppingWriteUnique { false /*  */ },
        basinHoppingStopEnergy { -DBL_MAX /*  */ },
        // [Global Optimization] //
        globalOptimizationMoveMethod { "md"s /*  */ },
        globalOptimizationDecisionMethod { "npew"s /*  */ },
        globalOptimizationSteps { 10000 /*  */ },
        globalOptimizationBeta { 1.05 /*  */ },
        globalOptimizationAlpha { 1.02 /*  */ },
        globalOptimizationMdmin { 3 /*  */ },
        globalOptimizationTargetEnergy { -1.E50 /*  */ },
        // [Monte Carlo] //
        monteCarloStepSize { 0.005 /*  */ },
        monteCarloSteps { 1000 /*  */ },
        // [BGSD] //
        alpha { 10.0 /*  */ },
        beta { 0.2 /*  */ },
        gradientfinitedifference { 0.000001 /*  */ },
        Hforceconvergence { 0.01 /*  */ },
        grad2energyconvergence { 0.000001 /*  */ },
        grad2forceconvergence { 0.0001 /*  */ }
        {}
          // clang-format on
    ~Parameters(){};
    // bool load(std::string filename);
    int load(std::string filename);
    int load(FILE *file);

    /** string constants: declared here, defined in Parameters.cpp. **/

    // potentials //
    // jobs //

    // Physical Constants
    double kB;
    double timeUnit;

    /** input parameters **/

    // [Main] //
    std::string job;
    double randomSeed;
    double temperature;
    bool checkpoint;
    bool quiet;
    bool writeLog;
    std::string iniFilename;
    std::string conFilename;
    double finiteDifference;
    size_t maxForceCalls;
    bool removeNetForce;

    // [Prefactor] //
    double prefactorDefaultValue;
    double prefactorMaxValue;
    double prefactorMinValue;
    double prefactorWithinRadius;
    double prefactorMinDisplacement;
    std::string prefactorRate;
    std::string prefactorConfiguration;
    bool prefactorAllFreeAtoms;
    std::string prefactorFilterScheme;
    double prefactorFilterFraction;

    // [Potential] //
    std::string potential;
    // MPI stuff, not actually specified in config file
    // it is used to pass information to the GPAW MPI potential.
    double MPIPollPeriod;
    int MPIPotentialRank;
    bool LogPotential;
    bool LAMMPSLogging;
    int LAMMPSThreads;
    bool EMTRasmussen;
    std::string extPotPath;

    // [AMS] and [AMS_IO] //
    std::string engine;
    std::string forcefield;
    std::string model;
    std::string xc;
    std::string basis;
    std::string resources;

    // [AMS_ENV] //
    std::string amshome;
    std::string scm_tmpdir;
    std::string scm_pythondir;
    std::string amsbin;
    std::string scmlicense;
    std::string amsresources;

    // [Structure Comparison] //
    double distanceDifference;
    double neighborCutoff;
    bool checkRotation;
    bool indistinguishableAtoms;
    double energyDifference;
    bool removeTranslation;

    // [Debug] //
    bool writeMovies;
    long int writeMoviesInterval;

    // [Saddle Search] //
    std::string saddleDisplaceType;
    std::string saddleMethod;
    std::string saddleMinmodeMethod;
    double saddleMaxEnergy;
    size_t saddleMaxIterations;
    double saddleDisplaceRadius;
    double saddleDisplaceMagnitude;
    double saddleMaxSingleDisplace;
    bool saddleNonnegativeDisplacementAbort;
    size_t saddleNonlocalCountAbort;
    double saddleNonlocalDistanceAbort;
    bool saddleRemoveRotation;
    float saddlePerpForceRatio;
    bool saddleConfinePositive;
    bool saddleBowlBreakout;
    size_t saddleBowlActive;
    double saddleConfinePositiveMinForce;
    double saddleConfinePositiveScaleRatio;
    double saddleConfinePositiveBoost;
    size_t saddleConfinePositiveMinActive;
    double saddleDynamicsTemperature;
    double saddleDynamicsStateCheckIntervalInput;
    double saddleDynamicsStateCheckInterval;
    double saddleDynamicsRecordIntervalInput;
    bool saddleDynamicsLinearInterpolation;
    double saddleDynamicsMaxInitCurvature;
    double saddleZeroModeAbortCurvature;

    // Not set.. old?
    long saddleMaxJumpAttempts;
    double saddleConvergedForce;
    double saddleDynamicsRecordInterval;
     // end unset

    // [Optimizer] //
    std::string optMethod;
    std::string optConvergenceMetric;
    std::string optConvergenceMetricLabel;
    size_t optMaxIterations;
    double optConvergedForce;
    double optMaxMove;
    double optTimeStepInput;
    double optTimeStep;
    double optMaxTimeStepInput;
    size_t optLBFGSMemory;
    double optLBFGSInverseCurvature;
    bool optLBFGSAutoScale;
    bool optLBFGSAngleReset;
    bool optLBFGSDistanceReset;
    bool optQMSteepestDecent;
    bool optCGNoOvershooting;
    bool optCGKnockOutMaxMove;
    double optCGLineConverged;
    bool optCGLineSearch;
    size_t optCGMaxIterBeforeReset;
    size_t optCGLineSearchMaxIter;
    double optSDAlpha;
    bool optSDTwoPoint;

    // Also unused?
    double optMaxTimeStep;
    double optLBFGSMaxInverseCurvature;
    // end unused

    // [Process Search] //
    bool processSearchMinimizeFirst;
    double processSearchMinimizationOffset;

    // [Dimer] //
    double dimerRotationAngle;
    bool dimerImproved;
    double dimerConvergedAngle;
    std::string dimerOptMethod;
    double dimerTorqueMin;
    double dimerTorqueMax;
    size_t dimerRotationsMin;
    size_t dimerRotationsMax;
    size_t dimerMaxIterations;
    bool dimerRemoveRotation;

    // [Lanczos] //
    double lanczosTolerance;
    size_t lanczosMaxIterations;
    bool lanczosQuitEarly;

    // [GPR Dimer] //
    double gprDimerRotationAngle;
    double gprDimerConvergedAngle;
    double gprDimerRelaxConvAngle;
    size_t gprDimerInitRotationsMax;
    size_t gprDimerRelaxRotationsMax;
    double gprDimerDivisorTdimerGP;
    size_t gprDimerMaxOuterIterations;
    size_t  gprDimerMaxInnerIterations;
    double gprDimerMidpointMaxDisp;
    std::string gprDimerRotOptMethod;
    std::string gprDimerTransOptMethod;
    double gprActiveRadius;
    double gprDimerSep;
    double gprDimerConvStep;
    double gprDimerMaxStep;
    double gprForceThreshold;
    double gprDimerRatioAtLimit;
    bool gprDimerInitRotGP;
    bool gprDimerInitTransGP;
    bool gprDimerManyIterations;

    // GPR Params
    std::string gprDimerHyperOptMethod;
    double gprDimerSigma2;
    double gprDimerJitterSigma2;
    double gprDimerNoiseSigma2;
    double gprDimerPriorMu;
    double gprDimerPriorSigma2;
    double gprDimerPriorNu;

    // GPR Optimization Parameters
    bool gprOptCheckDerivatives;
    size_t gprOptMaxIterations;
    double gprOptTolFunc;
    double gprOptTolSol;
    double gprOptLambdaLimit;
    double gprOptLambdaInit;
    // GPR Prune Parameters
    bool gprUsePrune;
    size_t gprPruneBegin;
    size_t gprPruneNVals;
    double gprPruneThreshold;

    // GPR Debugging Parameters
    size_t gprReportLevel;
    size_t  gprDebugLevel;
    std::string gprDebugOutDir;
    std::string gprDebugPosFile;
    std::string gprDebugEnergyFile;
    std::string gprDebugGradFile;
    std::string gprDebugOutExt;
    float gprDebugOffsetMidPoint;
    float gprDebugDy;
    float gprDebugDz;

    // [Hessian] //
    std::string hessianAtomList;
    double hessianZeroFreqValue;

    // [Nudged Elastic Band] //
    size_t nebImages;
    double nebSpring;
    bool nebClimbingImageMethod;
    bool nebClimbingImageConvergedOnly;
    bool nebOldTangent;
    size_t nebMaxIterations;
    bool nebDoublyNudged;
    bool nebDoublyNudgedSwitching;
    bool nebElasticBand;
    double nebConvergedForce;
    // Not set?
    std::string nebOptMethod;

    // [Molecular Dynamics] //
    double mdTimeStepInput;
    double mdTimeInput;
    // Not set?
    double mdTimeStep;
    double mdTime;
    long mdSteps;

    // [Thermostat] //
    std::string thermostat;
    double thermoAndersenAlpha;
    double thermoAndersenTcolInput;
    double thermoNoseMass;
    double thermoLangevinFrictionInput;
    // Not set?
    double thermoAndersenTcol;
    double thermoLangevinFriction;
    // std::vector<int> thermoAtoms;

    // [Parallel Replica] //
    bool parrepRefineTransition;
    bool parrepAutoStop;
    bool parrepDephaseLoopStop;
    double parrepDephaseTimeInput;
    size_t parrepDephaseLoopMax;
    double parrepStateCheckIntervalInput;
    double parrepRecordIntervalInput;
    double parrepCorrTimeInput;
    // Not set?
    double parrepStateCheckInterval;
    double parrepDephaseTime;
    double parrepRecordInterval;
    double parrepCorrTime;

    // [Temperature Accelerated Dynamics] //
    double tadLowT;
    double tadMinPrefactor;
    double tadConfidence;

    // [Replica Exchange] //
    std::string repexcTemperatureDistribution;
    size_t repexcReplicas;
    size_t repexcExchangeTrials;
    double repexcSamplingTimeInput;
    double repexcTemperatureLow;
    double repexcTemperatureHigh;
    double repexcSamplingTime;
    double repexcExchangePeriod;
    // Not set
    double repexcExchangePeriodInput;

    // [Bond Boost] // or Hyperdynamics
    std::string biasPotential;
    std::string bondBoostBALS;
    double bondBoostDVMAX;
    double bondBoostQRR;
    double bondBoostPRR;
    double bondBoostQcut; // Not set
    double bondBoostRMDTimeInput;
    double bondBoostRMDTime; // Not set

    // [Basin Hopping] //
    double basinHoppingDisplacement;
    double basinHoppingPushApartDistance;
    double basinHoppingInitialRandomStructureProbability;
    size_t basinHoppingSteps;
    size_t basinHoppingQuenchingSteps;
    bool basinHoppingSingleAtomDisplace;
    bool basinHoppingSignificantStructure;
    std::string basinHoppingDisplacementAlgorithm;
    std::string basinHoppingDisplacementDistribution;
    double basinHoppingSwapProbability;
    size_t basinHoppingJumpMax;
    size_t basinHoppingJumpSteps;
    bool basinHoppingAdjustDisplacement;
    size_t basinHoppingAdjustPeriod;
    double basinHoppingAdjustFraction;
    double basinHoppingTargetRatio;
    bool basinHoppingWriteUnique;
    double basinHoppingStopEnergy;

    // [Global Optimization] //
    std::string globalOptimizationMoveMethod;
    std::string globalOptimizationDecisionMethod;
    size_t globalOptimizationSteps;
    double globalOptimizationBeta;
    double globalOptimizationAlpha;
    double globalOptimizationMdmin;
    double globalOptimizationTargetEnergy;

    // [Monte Carlo] //
    double monteCarloStepSize;
    size_t monteCarloSteps;

    // [BGSD] //

    double alpha;
    double beta;
    double gradientfinitedifference;
    double Hforceconvergence;
    double grad2energyconvergence;
    double grad2forceconvergence;

#ifdef EONMPI
    MPI_Comm MPIClientComm;
#endif

private:
    std::string toLowerCase(string s);
};

#endif

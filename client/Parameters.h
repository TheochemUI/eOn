/*
** This file is part of eON.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eON Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eON
*/
#pragma once

#include "BaseStructures.h"
#include <string>

#ifdef EONMPI
#include "mpi.h"
#endif

/** Contains all runtime parameters and results. No functionality just
 * bookkeeping.*/
class Parameters {

public:
  Parameters();
  ~Parameters() = default;
  Parameters(const Parameters &) = default;
  int load(const std::string &filename);

  // Physical Constants
  double kB;
  double timeUnit;

  // Sections for better organization
  struct Main {
    JobType job;
    long randomSeed;
    double temperature;
    bool quiet;
    bool writeLog;
    bool checkpoint;
    std::string inpFilename;
    std::string conFilename;
    double finiteDifference;
    long maxForceCalls;
    bool removeNetForce;
  } main;

  struct Potential {
    PotType potential;
    double MPIPollPeriod;
    bool LAMMPSLogging;
    int LAMMPSThreads;
    bool EMTRasmussen;
    bool LogPotential;
    std::string extPotPath;
  } pot;

  struct AMS {              // Also for AMS_IO
    std::string engine;     // MOPAC, ADF, BAND, REAXFF, FORCEFIELD
    std::string forcefield; // OPt.ff etc. (REAXFF)
    std::string model;      // Model hamiltonian (MOPAC)
    std::string resources;  // DFTB
    std::string xc;         // Exchange (BAND, ADF)
    std::string basis;      // Basis (BAND, ADF)
  } ams;

  struct AMS_ENV {
    std::string amshome;
    std::string scm_tmpdir;
    std::string scmlicense;
    std::string scm_pythondir;
    std::string amsbin;
    std::string amsresources;
  } amsenv;

  struct XTBPot {
    std::string paramset;
    double elec_temperature;
    size_t maxiter;
    double acc;
  } xtbpot;

  struct StructureComparison {
    double
        distanceDifference; ///< The distance criterion for comparing geometries
    double
        neighborCutoff; ///< radius used in the local atomic structure analysis
    bool checkRotation;
    bool indistinguishableAtoms;
    double energyDifference;
    bool removeTranslation;
  } structcomp;

  struct ProcessSearch {
    bool minimizeFirst;
    double minimizationOffset; ///< how far from the saddle to displace the
                               ///< minimization images
  } procsearch;

  struct SaddleSearch {
    long maxJumpAttempts; /**< number of displacements to reach a convex
                          region;  if 0, a search is started after the
                          displacement */
                          // TODO(rg): Not used
    long maxIterations;   ///< max iterations for saddle point searches and
                          ///< minimization
    std::string method;
    std::string minmodeMethod; /**< algorithm to be used for lowest eigenmode
                               determination */
    std::string displaceType;  ///< displacement type to use
    double maxEnergy;          /**< energy above product state that will cause
                               termination of the saddle point search */
    double displaceMagnitude;  ///< norm of the displacement vector
    double maxSingleDisplace;  /**< maximum value of displacement in x, y and z
                               direction for atoms being displaced */
    double displaceRadius;     /**< atoms within this radius of the displacement
                               atoms are also displaced */
    double convergedForce;     /**< force convergence criterion required for a
                               saddle point search */
    double perpForceRatio; /**< proportion to keep of the perpendicular force
                           when the lowest eigenvalue is positive */
    bool nonnegativeDisplacementAbort; /** abort the saddle search if the
                                       displacement does not have a
                                       negative mode */
    long nonlocalCountAbort;      /**< abort the search if this many atoms move
                                  more than NonlocalDistanceAbort */
    double nonlocalDistanceAbort; /**< abort the search if NonlocalCountAbort
                                  atoms move more than this distance */
    bool removeRotation; /**< remove dominant rotational component when system
                         is translated */
    double
        dynamicsTemperature; ///< temperature for dynamics saddle search method
    double dynamicsStateCheckIntervalInput;
    double dynamicsStateCheckInterval; ///< how often to minimize
    double dynamicsRecordIntervalInput;
    double dynamicsRecordInterval;
    bool dynamicsLinearInterpolation;
    double dynamicsMaxInitCurvature;
    bool confinePositive;
    bool bowlBreakout;
    long bowlActive;
    double confinePositiveMinForce;
    double confinePositiveScaleRatio;
    double confinePositiveBoost;
    long confinePositiveMinActive;
    double zeroModeAbortCurvature;
  } saddle;

  struct Optimizer {
    OptType method;
    OptType refineOptMethod;       ///< used below refine threshold
    double refineThreshold;        ///< threshold to switch opt_method
    std::string convergenceMetric; ///< norm, max_atom, max_component
    std::string convergenceMetricLabel;
    size_t maxIterations; /**< maximum iterations for saddle point searches and
                          minimization */
    double
        maxMove; ///< maximum displacement vector for a step during optimization
    double convergedForce; ///< force convergence criterion
    double timeStepInput;
    double timeStep;         ///< QuickMin timestep
    double maxTimeStepInput; ///< FIRE timestep (input)
    double maxTimeStep;      ///< FIRE max timestep
    long LBFGSMemory; ///< number of previous forces to keep in the bfgs memory
    double LBFGSInverseCurvature;
    double LBFGSMaxInverseCurvature;
    bool LBFGSAutoScale;
    bool LBFGSAngleReset;
    bool LBFGSDistanceReset;
    bool QMSteepestDecent;  /**< if set the velocity will always be set to zero
                            in  quickmin */
    bool CGNoOvershooting;  /**< if set it is ensured that the approximate line
                            search in conjugate gradients never overshoot the
                            minimum along the search line */
    bool CGKnockOutMaxMove; /**< if set the old search direction is nullified
                            when steps larger than the maxMove are conducted */
    bool CGLineSearch;      ///< if set full line search is conducted
    double CGLineConverged; /**< convergence criteria for line search, ratio
                            between force component along search line and
                            the orthogonal part */
    long CGLineSearchMaxIter;  ///< maximal nr of iterations during line search
    long CGMaxIterBeforeReset; /**< max nr of cg steps before reset, if 0 no
                               resetting is done */
    double SDAlpha;
    bool SDTwoPoint;
  } optim;

  struct Dimer {
    double rotationAngle;  ///< finite difference rotation angle
    bool improved;         ///< turn on the improved dimer method
    double convergedAngle; ///< stop rotating when angle drops below this value
    long maxIterations;    ///< maximum number of rotation iterations
    std::string optMethod; ///< method to determine the next rotation direction
    bool removeRotation;   /** remove dominant rotational component when
                           estimating the eigenmode */
    // These are only for the older variant
    long rotationsMax;
    long rotationsMin;
    double torqueMax;
    double torqueMin;
  } dimer;

  struct GPRDimer {
    double rotationAngle;       ///< finite difference rotation angle
    double convergedAngle;      /**< stop rotating when angle drops below this
                                value {T_anglerot_init} */
    double relaxConvAngle;      /** stop rotating when angle drops below this
                                value during relaxation {T_anglerot_gp} */
    long initRotationsMax;      ///< {num_iter_initrot}
    long relaxRotationsMax;     ///< {num_iter_rot_gp}
    long divisorTdimerGP;       ///< {divisor_T_dimer_gp}
    long maxOuterIterations;    /**< maximum number of outer iterations or new
                                sets of observations {num_bigiter} */
    long maxInnerIterations;    /** maximum number of steps during the
                                relaxation phase {num_iter} */
    double midpointMaxDisp;     ///< {disp_max}
    std::string rotOptMethod;   /** method to determine the next rotation
                                direction {method_rot} */
    std::string transOptMethod; /**< method to determine the next rotation
                                direction {method_trans} */
    double activeRadius; /**< activation radius for inclusion in covariance
                         matrix {actdist_fro} */
    double dimerSep; /**< distance from the middle point of the dimer to the two
                     images {dimer_sep} */
    double convStep; ///< step length for convex regions {param_trans[0]}
    double maxStep;  ///< maximum step length {param_trans[1]}
    double forceThreshold; /** maximum component of the force acting on the
                           middle point of the dimer (stops when accurate
                           force components are below this value) {T_dimer} */
    double ratioAtLimit;   /** {ratio_at_limit} defines the limit for the
                           ratio of inter-atomic distances between the
                           image and its "nearest observed data point" */
    bool initRotGP;        /** initial rotations without GP (1) or with GP (0)
                           {initrot_nogp} */
    bool initTransGP;      /** initial translations without GP (1) or with GP
                           (0) {inittrans_nogp} */
    bool manyIterations;   /** indicates of the number of iterations is
                           larger than required for dimer convergence on
                           the accurate energy surface (1), otherwise (0)
                           the relaxation phase is continued from the
                           current dimer if the maximum iterations in the
                           relaxation phase is reached {islarge_num_iter} */
    std::string hyperOptMethod; ///< method to optimize hyperparameters
                                ///< {optimization_alg}
    double sigma2;              ///< GPR variance {gp_sigma2}
    double jitterSigma2;        ///< GPR jitter variance {jitter_sigma2}
    double noiseSigma2;         ///< noise Variance {sigma2}
    double priorMu;             ///< prior mean {prior_mu}
    double priorSigma2;         ///< prior variance {prior_s2}
    long priorNu;               ///< prior degrees of freedom {prior_nu}
    // GPR Optimization Params
    bool optCheckDerivatives; ///< {check_derivative}
    int optMaxIterations;     ///< {max_iter}
    double optTolFunc;        ///< {tolerance_func}
    double optTolSol;         ///< {tolerance_sol}
    long optLambdaLimit;      ///< {lambda_limit}
    long optLambdaInit;       ///< {lambda}
    bool usePrune;            ///< {use_prune}
    int pruneBegin;           ///< {start_prune_at}
    int pruneNVals;           ///< {nprune_vals}
    double pruneThreshold;    ///< {prune_threshold}
    // GPR Debugging
    int reportLevel;             ///< {report_level}
    int debugLevel;              ///< {debug_level}
    std::string debugOutDir;     ///< {debug_output_dir}
    std::string debugPosFile;    ///< {debug_output_file_R}
    std::string debugEnergyFile; ///< {debug_output_file_E}
    std::string debugGradFile;   ///< {debug_output_file_G}
    std::string debugOutExt;     ///< {debug_output_file_extension}
    double debugOffsetMidPoint;  ///< {debug_offset_from_mid_point}
    double debugDy;              ///< {debug_dy}
    double debugDz;              ///< {debug_dz}
  } gprd;

  struct Surrogate {
    bool use;
    JobType sub_job;
    double gp_uncertainity;
    bool gp_linear_path_always;
    PotType potential; // ONLY: catlearn for now
  } surrogate;

  struct CatLearn {
    std::string path;
    std::string model;
    std::string prior;
    bool use_deriv;
    bool use_fingerprint;
    bool parallel;
  } catl;

  struct ASEOrca {
    std::string orca_path;
    std::string orca_nproc;
    std::string simpleinput;
  } aseorca;

  struct Lanczos {
    double tolerance;   /** difference between the lowest eignevalues of two
                        successive iterations */
    long maxIterations; ///<< maximum number of iterations
    bool quitEarly;
  } lanczos;

  struct Prefactor {
    double defaultValue;    ///< default prefactor; calculate explicitly if zero
    double maxValue;        ///< max prefactor allowed
    double minValue;        ///< min prefactor allowed
    double withinRadius;    /** atoms within this radius of the displaced
                            atoms are put in the Hessian, unless
                            filterMode is fraction */
    double minDisplacement; /** atoms with displacement between min1 or min2
                            and the saddle point are put in the Hessian */
    ::Prefactor::RATE rate; ///< method to estimate prefactor
    ::Prefactor::TYPE configuration; /** configuration for which the frequencies
                               should be determined */
    bool allFreeAtoms; ///< use all free atom when determining the prefactor
    ::Prefactor::FILTER filterScheme; /** "cutoff" or "fraction", which use
                              prefactorMinDisplacement or
                              prefactorFilterFraction, respectively */
    double filterFraction; /** Include atoms whose summed motion comprise
                           more than prefactorFilterFraction in the
                           prefactor calculation. Prioritizes atoms
                           that move more. */
  } prefactor;

  struct Hessian {
    std::string atomList;
    double zeroFreqValue;
  } hessian;

  struct NudgedElasticBand {
    std::string optMethod;
    size_t images;
    size_t maxIterations;
    double spring;
    bool climbingImageMethod;
    bool climbingImageConvergedOnly;
    bool oldTangent;
    bool doublyNudged;
    bool doublyNudgedSwitching;
    bool elasticBand;
    double convergedForce; /** force convergence criterion required for an
                           optimization */
    // For energy weighted
    bool energyWeighted;
    double KSPMin;
    double KSPMax;
  } neb;

  struct MolDynamics {
    double timeStepInput;
    double timeStep;
    double timeInput;
    double time;
    long steps;
  } md;

  struct Thermostat {
    std::string kind; // thermostat
    double andersenAlpha;
    double andersenTcolInput;
    double andersenTcol;
    double noseMass;
    double langevinFrictionInput;
    double langevinFriction;
    // std::vector<int> thermoAtoms;
  } thermostat;

  struct ParallelReplica {
    bool refineTransition;
    bool autoStop;
    bool dephaseLoopStop;
    double dephaseTimeInput;
    double dephaseTime;
    long dephaseLoopMax;
    double stateCheckIntervalInput;
    double stateCheckInterval;
    double recordIntervalInput;
    double recordInterval;
    double corrTimeInput;
    double corrTime;
  } parrep;

  // [Temperature Accelerated Dynamics] //
  struct TAD {
    double lowT;
    double minPrefactor;
    double confidence;
  } tad;

  struct ReplicaExchange {
    std::string temperatureDistribution;
    size_t replicas;
    size_t exchangeTrials;
    double samplingTimeInput;
    double samplingTime;
    double temperatureHigh;
    double temperatureLow;
    double exchangePeriodInput;
    double exchangePeriod;
  } repexc;

  struct Hyperdynamics {
    std::string biasPotential;
    std::string BALS;
    double RMDTimeInput;
    double RMDTime;
    double DVMAX;
    double QRR;
    double PRR;
    double Qcut;
  } bondBoost;

  struct BasinHopping {
    double displacement;
    double initialRandomStructureProbability;
    double pushApartDistance;
    long steps;
    long quenchingSteps;
    bool significantStructure;
    bool singleAtomDisplace;
    std::string displacementAlgorithm;
    std::string displacementDistribution;
    double swapProbability;
    long jumpMax;
    long jumpSteps;
    bool adjustDisplacement;
    long adjustPeriod;
    double adjustFraction;
    double targetRatio;
    bool writeUnique;
    double stopEnergy;
  } bhop;

  struct GlobalOptimization {
    std::string moveMethod;
    std::string decisionMethod;
    size_t steps;
    double beta;
    double alpha;
    long mdmin;
    double targetEnergy;
  } globopt;

  struct MonteCarlo {
    double stepSize;
    size_t steps;
  } monte_carlo;

  struct BGSD {
    double alpha;
    double beta;
    double gradientfinitedifference;
    double Hforceconvergence;
    double grad2energyconvergence;
    double grad2forceconvergence;
  } bgsd;

  struct Debug {
    bool writeMovies;
    long writeMoviesInterval;
  } debug;

  // MPI stuff, not actually specified in config file
  // it is used to pass information to the GPAW MPI potential.
  int MPIPotentialRank;
#ifdef EONMPI
  MPI_Comm MPIClientComm;
#endif

private:
  std::string toLowerCase(std::string s);
};

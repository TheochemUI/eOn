#ifndef PARAMETERS_H
#define PARAMETERS_H

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
  double kB;
  double timeUnit;

  /** input parameters **/

  // [Main] //
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

  // [Potential] //
  PotType potential;
  double MPIPollPeriod;
  bool LAMMPSLogging;
  int LAMMPSThreads;
  bool EMTRasmussen;
  bool LogPotential;
  string extPotPath;

  // [AMS] and [AMS_IO] //
  string engine;     // MOPAC, ADF, BAND, REAXFF, FORCEFIELD
  string forcefield; // OPt.ff etc. (REAXFF)
  string model;      // Model hamiltonian (MOPAC)
  string resources;  // DFTB
  string xc;         // Exchange (BAND, ADF)
  string basis;      // Basis (BAND, ADF)

  // [AMS_ENV] //
  string amshome;
  string scm_tmpdir;
  string scmlicense;
  string scm_pythondir;
  string amsbin;
  string amsresources;

  // [XTBPot] //
  string xtb_paramset;
  double xtb_elec_temperature;
  size_t xtb_maxiter;
  double xtb_acc;

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
  double distanceDifference; // The distance criterion for comparing geometries
  double neighborCutoff; // radius used in the local atomic structure analysis
  bool checkRotation;
  bool indistinguishableAtoms;
  double energyDifference;
  bool removeTranslation;

  // [Process Search] //
  bool processSearchMinimizeFirst;
  double processSearchMinimizationOffset; // how far from the saddle to displace
                                          // the minimization images

  // [Saddle Search] //
  long saddleMaxJumpAttempts; // number of displacements to reach a convex
                              // region;  if 0, a search is started after the
                              // displacement
                              // TODO(rg): Not used
  long saddleMaxIterations;   // max iterations for saddle point searches and
                              // minimization
  string saddleMethod;
  string saddleMinmodeMethod;     // algorithm to be used for lowest eigenmode
                                  // determination
  string saddleDisplaceType;      // displacement type to use
  double saddleMaxEnergy;         // energy above product state that will cause
                                  // termination of the saddle point search
  double saddleDisplaceMagnitude; // norm of the displacement vector
  double saddleMaxSingleDisplace; // maximum value of displacement in x, y and z
                                  // direction for atoms being displaced
  double saddleDisplaceRadius; // atoms within this radius of the displacement
                               // atoms are also displaced
  double saddleConvergedForce; // force convergence criterion required for a
                               // saddle point search
  double saddlePerpForceRatio; // proportion to keep of the perpendicular force
                               // when the lowest eigenvalue is positive
  bool saddleNonnegativeDisplacementAbort; // abort the saddle search if the
                                           // displacement does not have a
                                           // negative mode
  long saddleNonlocalCountAbort; // abort the search if this many atoms move
                                 // more than NonlocalDistanceAbort
  double saddleNonlocalDistanceAbort; // abort the search if NonlocalCountAbort
                                      // atoms move more than this distance
  bool saddleRemoveRotation; // remove dominant rotational component when system
                             // is translated

  double saddleDynamicsTemperature; // temperature for dynamics saddle search
                                    // method
  double saddleDynamicsStateCheckIntervalInput;
  double saddleDynamicsStateCheckInterval; // how often to minimize
  double saddleDynamicsRecordIntervalInput;
  double saddleDynamicsRecordInterval;
  bool saddleDynamicsLinearInterpolation;
  double saddleDynamicsMaxInitCurvature;

  bool saddleConfinePositive;             // undocumented
  bool saddleBowlBreakout;                // undocumented
  long saddleBowlActive;                  // undocumented
  double saddleConfinePositiveMinForce;   // undocumented
  double saddleConfinePositiveScaleRatio; // undocumented
  double saddleConfinePositiveBoost;      // undocumented
  long saddleConfinePositiveMinActive;    // undocumented
  double saddleZeroModeAbortCurvature;

  // [Optimizer] //
  OptType optMethod;
  OptType refineOptMethod;     // used below refine threshold
  double refineThreshold;      // threshold to switch opt_method
  string optConvergenceMetric; // norm, max_atom, max_component
  string optConvergenceMetricLabel;
  size_t optMaxIterations; // maximum iterations for saddle point searches and
                           // minimization
  double
      optMaxMove; // maximum displacement vector for a step during optimization
  double optConvergedForce; // force convergence criterion required for an
                            // optimization
  double optTimeStepInput;
  double optTimeStep;         // time step size used in quickmin
  double optMaxTimeStepInput; // maximum time step for FIRE.
  double optMaxTimeStep;      // maximum time step for FIRE.
  long optLBFGSMemory; // number of previous forces to keep in the bfgs memory
  double optLBFGSInverseCurvature;
  double optLBFGSMaxInverseCurvature;
  bool optLBFGSAutoScale;
  bool optLBFGSAngleReset;
  bool optLBFGSDistanceReset;
  bool optQMSteepestDecent; // if set the velocity will always be set to zero in
                            // quickmin
  bool optCGNoOvershooting; // if set it is ensured that the approximate line
                            // search in conjugate gradients never overshoot the
                            // minimum along the search line
  bool
      optCGKnockOutMaxMove; // if set the old search direction is nullified when
                            // steps larger than the optMaxMove are conducted
  bool optCGLineSearch;     // if set full line search is conducted
  double optCGLineConverged;    // convergence criteria for line search, ratio
                                // between force component along search line and
                                // the orthogonal part
  long optCGLineSearchMaxIter;  // maximal nr of iterations during line search
  long optCGMaxIterBeforeReset; // max nr of cg steps before reset, if 0 no
                                // resetting is done
  double optSDAlpha;
  bool optSDTwoPoint;

  // [Dimer] //
  double dimerRotationAngle;  // finite difference rotation angle
  bool dimerImproved;         // turn on the improved dimer method
  double dimerConvergedAngle; // stop rotating when angle drops below this value
  long dimerMaxIterations;    // maximum number of rotation iterations
  string dimerOptMethod;      // method to determine the next rotation direction
  long dimerRotationsMax;     // old
  long dimerRotationsMin;     // old
  double dimerTorqueMax;      // old
  double dimerTorqueMin;      // old
  bool dimerRemoveRotation;   // remove dominant rotational component when
                              // estimating the eigenmode

  // [GPR Dimer] //
  double gprDimerRotationAngle;    // finite difference rotation angle
  double gprDimerConvergedAngle;   // stop rotating when angle drops below this
                                   // value {T_anglerot_init}
  double gprDimerRelaxConvAngle;   // stop rotating when angle drops below this
                                   // value during relaxation {T_anglerot_gp}
  long gprDimerInitRotationsMax;   // {num_iter_initrot}
  long gprDimerRelaxRotationsMax;  // {num_iter_rot_gp}
  long gprDimerDivisorTdimerGP;    // {divisor_T_dimer_gp}
  long gprDimerMaxOuterIterations; // maximum number of outer iterations or new
                                   // sets of observations {num_bigiter}
  long gprDimerMaxInnerIterations; // maximum number of steps during the
                                   // relaxation phase {num_iter}
  double gprDimerMidpointMaxDisp;  // {disp_max}
  string gprDimerRotOptMethod;     // method to determine the next rotation
                                   // direction {method_rot}
  string gprDimerTransOptMethod;   // method to determine the next rotation
                                   // direction {method_trans}
  double gprActiveRadius; // activation radius for inclusion in covariance
                          // matrix {actdist_fro}
  double gprDimerSep; // distance from the middle point of the dimer to the two
                      // images {dimer_sep}
  struct gprd_translation_t {
    double step_length;     // step length for convex regions {param_trans[0]}
    double max_step_length; // maximum step length {param_trans[1]}
    double rotrem_thresh;   // {rotation_removal_projection_threshold}
  } gprd_trans_options;
  struct early_stopping_t {
    std::string dist_metrics;
    double threshold;
  } early_stopping_options;
  double gprDimerRatioAtLimit; // {ratio_at_limit} defines the limit for the
                               // ratio of inter-atomic distances between the
                               // image and its "nearest observed data point"
  bool gprDimerInitRotGP;   // initial rotations without GP (1) or with GP (0)
                            // {initrot_nogp}
  bool gprDimerManyIterations; // indicates of the number of iterations is
                               // larger than required for dimer convergence on
                               // the accurate energy surface (1), otherwise (0)
                               // the relaxation phase is continued from the
                               // current dimer if the maximum iterations in the
                               // relaxation phase is reached {islarge_num_iter}
  // GPR Params
  double gprDimerSigma2;       // GPR variance {gp_sigma2}
  double gprDimerJitterSigma2; // GPR jitter variance {jitter_sigma2}
  double gprDimerNoiseSigma2;  // noise Variance {sigma2}
  double gprDimerPriorMu;      // prior mean {prior_mu}
  double gprDimerPriorSigma2;  // prior variance {prior_s2}
  long gprDimerPriorNu;        // prior degrees of freedom {prior_nu}
  // GPR Optimization Parameters
  bool gprUsePrune;         // {use_prune}
  int gprPruneBegin;        // {start_prune_at}
  int gprPruneNVals;        // {nprune_vals}
  double gprPruneThreshold; // {prune_threshold}
  struct gpr_fps_t {
    std::string metric;
    int history;
  } fps_options;
  struct gpr_hypopt_t {
    std::string hopt_method;
    bool check_derivative;
    double tol_func;
    double tol_sol;
    int max_iter;

    // SCG
    struct scg_t {
      double lambda_limit;
      double lambda;
    } scg;

    // ADAM
    struct adam_t {
      double lr;
      double lrd;
      double b1;
      double b2;
      double eps;
      double weight_decay;
      bool amsgrad;
    } adam;
  } gpr_hypopt_options;

  // GPR Debugging Parameters
  int gprReportLevel;            // {report_level}
  int gprDebugLevel;             // {debug_level}
  string gprDebugOutDir;         // {debug_output_dir}
  string gprDebugPosFile;        // {debug_output_file_R}
  string gprDebugEnergyFile;     // {debug_output_file_E}
  string gprDebugGradFile;       // {debug_output_file_G}
  string gprDebugOutExt;         // {debug_output_file_extension}
  double gprDebugOffsetMidPoint; // {debug_offset_from_mid_point}
  double gprDebugDy;             // {debug_dy}
  double gprDebugDz;             // {debug_dz}

  // GP Surrogate Parameters
  bool use_surrogate;
  JobType sub_job;
  double gp_uncertainity;
  bool gp_linear_path_always;
  PotType surrogatePotential; // ONLY: catlearn for now

  // [CatLearn]
  std::string catl_path;
  std::string catl_model;
  std::string catl_prior;
  bool catl_use_deriv;
  bool catl_use_fingerprint;
  bool catl_parallel;

  // [ASE_ORCA] //
  std::string orca_path;
  std::string orca_nproc;
  std::string orca_sline; // Other catchall values

  // [ASE_NWCHEM] //
  std::string nwchem_path;
  std::string nwchem_nproc;
  std::string nwchem_multiplicity; // 1 for singlet, 2 for doublet
  double nwchem_scf_thresh;
  long nwchem_scf_maxiter;

  // [Metatomic] //
  struct metatomic_options_t {
    std::string model_path;  // Path to the TorchScript model file.
    std::string device;      // "cpu", "cuda", "mps", or empty to auto-detect.
    std::string length_unit; // The unit of length used in the simulation (e.g.,
                             // "angstrom").
    std::string extensions_directory; // Path for TorchScript extensions.
    bool check_consistency;           // To enable model's internal checks.
  } metatomic_options;

  // [Lanczos] //
  double lanczosTolerance;   // difference between the lowest eignevalues of two
                             // successive iterations
  long lanczosMaxIterations; // maximum number of iterations
  bool lanczosQuitEarly;

  // [Prefactor] //
  double
      prefactorDefaultValue; // default prefactor; calculate explicitly if zero
  double prefactorMaxValue;  // max prefactor allowed
  double prefactorMinValue;  // min prefactor allowed
  double prefactorWithinRadius; // atoms within this radius of the displaced
                                // atoms are put in the Hessian, unless
                                // filterMode is fraction
  double
      prefactorMinDisplacement;  // atoms with displacement between min1 or min2
                                 // and the saddle point are put in the Hessian
  string prefactorRate;          // method to estimate prefactor
  string prefactorConfiguration; // configuration for which the frequencies
                                 // should be determined
  bool
      prefactorAllFreeAtoms; // use all free atom when determining the prefactor
  string prefactorFilterScheme;   // "cutoff" or "fraction", which use
                                  // prefactorMinDisplacement or
                                  // prefactorFilterFraction, respectively.
  double prefactorFilterFraction; // Include atoms whose summed motion comprise
                                  // more than prefactorFilterFraction in the
                                  // prefactor calculation. Prioritizes atoms
                                  // that move more.

  // [Hessian] //
  string hessianAtomList;
  double hessianZeroFreqValue;

  // [Nudged Elastic Band] //
  long nebImages;
  long nebMaxIterations;
  double nebSpring;
  bool nebClimbingImageMethod;
  bool nebClimbingImageConvergedOnly;
  bool nebOldTangent;
  bool nebDoublyNudged;
  bool nebDoublyNudgedSwitching;
  string nebOptMethod;
  bool nebElasticBand;
  double nebConvergedForce; // force convergence criterion required for an
                            // optimization
  double nebciAfter;        // force convergence before ci-neb is used
  // For energy weighted
  double nebKSPMin;
  double nebKSPMax;
  bool nebEnergyWeighted;

  // For hybrid dimer-NEB
  bool nebciWithMMF;
  double nebciMMFAfter;
  long nebciMMFnSteps;

  // Initial path
  string nebIpath; // file containing list of .con files for the initial path
  // Minimize endpoints
  bool nebMinimEP;

  // [Molecular Dynamics] //
  double mdTimeStepInput;
  double mdTimeStep;
  double mdTimeInput;
  double mdTime;
  long mdSteps;

  // [Parallel Replica] //
  bool parrepRefineTransition;
  bool parrepAutoStop;
  bool parrepDephaseLoopStop;
  double parrepDephaseTimeInput;
  double parrepDephaseTime;
  long parrepDephaseLoopMax;
  double parrepStateCheckIntervalInput;
  double parrepStateCheckInterval;
  double parrepRecordIntervalInput;
  double parrepRecordInterval;
  double parrepCorrTimeInput;
  double parrepCorrTime;

  // [Temperature Accelerated Dynamics] //
  double tadLowT;
  double tadMinPrefactor;
  double tadConfidence;

  // [Thermostat] //
  string thermostat;
  double thermoAndersenAlpha;
  double thermoAndersenTcolInput;
  double thermoAndersenTcol;
  double thermoNoseMass;
  double thermoLangevinFrictionInput;
  double thermoLangevinFriction;
  // std::vector<int> thermoAtoms;

  // [Replica Exchange] //
  string repexcTemperatureDistribution;
  long repexcReplicas;
  long repexcExchangeTrials;
  double repexcSamplingTimeInput;
  double repexcSamplingTime;
  double repexcTemperatureHigh;
  double repexcTemperatureLow;
  double repexcExchangePeriodInput;
  double repexcExchangePeriod;

  // [Bond Boost] //
  string biasPotential;
  string bondBoostBALS;
  double bondBoostRMDTimeInput;
  double bondBoostRMDTime;
  double bondBoostDVMAX;
  double bondBoostQRR;
  double bondBoostPRR;
  double bondBoostQcut;

  // [Basin Hopping] //
  double basinHoppingDisplacement;
  double basinHoppingInitialRandomStructureProbability;
  double basinHoppingPushApartDistance;
  long basinHoppingSteps;
  long basinHoppingQuenchingSteps;
  bool basinHoppingSignificantStructure;
  bool basinHoppingSingleAtomDisplace;
  string basinHoppingDisplacementAlgorithm;
  string basinHoppingDisplacementDistribution;
  double basinHoppingSwapProbability;
  long basinHoppingJumpMax;
  long basinHoppingJumpSteps;
  bool basinHoppingAdjustDisplacement;
  long basinHoppingAdjustPeriod;
  double basinHoppingAdjustFraction;
  double basinHoppingTargetRatio;
  bool basinHoppingWriteUnique;
  double basinHoppingStopEnergy;

  // [Global Optimization] //
  string globalOptimizationMoveMethod;
  string globalOptimizationDecisionMethod;
  long globalOptimizationSteps;
  double globalOptimizationBeta;
  double globalOptimizationAlpha;
  long globalOptimizationMdmin;
  double globalOptimizationTargetEnergy;

  // [Monte Carlo] //
  double monteCarloStepSize;
  int monteCarloSteps;

  // [BGSD] //

  double alpha;
  double beta;
  double gradientfinitedifference;
  double Hforceconvergence;
  double grad2energyconvergence;
  double grad2forceconvergence;

  // MPI stuff, not actually specified in config file
  // it is used to pass information to the GPAW MPI potential.
  int MPIPotentialRank;
#ifdef EONMPI
  MPI_Comm MPIClientComm;
#endif

  // [Debug] //
  bool writeMovies;
  long writeMoviesInterval;
  bool estNEBeig;
  string nebMMF;

private:
  string toLowerCase(string s);
};

#endif

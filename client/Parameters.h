#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <ctype.h>
#include <vector>

using namespace std;

#ifdef EONMPI
    #include "mpi.h"
#endif

/** Contains all runtime parameters and results. No functionality just bookkeeping.*/
class Parameters {

public:
    Parameters();
    ~Parameters();
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
    string job;
    long   randomSeed;
    double temperature;
    bool   quiet;
    bool   writeLog;
    bool   checkpoint;
    string iniFilename;
    string conFilename;
    double finiteDifference;
    long   maxForceCalls;
    bool   removeNetForce;

    // [Potential] //
    string potential;
    double MPIPollPeriod;
    bool   LAMMPSLogging;
    int    LAMMPSThreads;
    bool   EMTRasmussen;
    bool   LogPotential;
    string extPotPath;

    // [Structure Comparison] //
    double distanceDifference; // The distance criterion for comparing geometries
    double neighborCutoff; // radius used in the local atomic structure analysis
    bool   checkRotation;
    bool   indistinguishableAtoms;
	double energyDifference;
    bool   removeTranslation;

    // [Process Search] //
    bool   processSearchMinimizeFirst;
    double processSearchMinimizationOffset; // how far from the saddle to displace the minimization images

    // [Saddle Search] //
    long   saddleMaxJumpAttempts; // number of displacements to reach a convex region;  if 0, a search is started after the displacement
    long   saddleMaxIterations; // max iterations for saddle point searches and minimization
    string saddleMethod;
    string saddleMinmodeMethod; // algorithm to be used for lowest eigenmode determination
    string saddleDisplaceType; // displacement type to use
    double saddleMaxEnergy; // energy above product state that will cause termination of the saddle point search
    double saddleDisplaceMagnitude; // norm of the displacement vector
    double saddleMaxSingleDisplace; // maximum value of displacement in x, y and z direction for atoms being displaced
    double saddleDisplaceRadius; // atoms within this radius of the displacement atoms are also displaced
    double saddleConvergedForce; // force convergence criterion required for a saddle point search
    double saddlePerpForceRatio; // proportion to keep of the perpendicular force when the lowest eigenvalue is positive 
    bool   saddleNonnegativeDisplacementAbort; // abort the saddle search if the displacement does not have a negative mode
    long   saddleNonlocalCountAbort; // abort the search if this many atoms move more than NonlocalDistanceAbort
    double saddleNonlocalDistanceAbort; // abort the search if NonlocalCountAbort atoms move more than this distance
    bool   saddleRemoveRotation; // remove dominant rotational component when system is translated

    double saddleDynamicsTemperature; //temperature for dynamics saddle search method
    double saddleDynamicsStateCheckIntervalInput; 
    double saddleDynamicsStateCheckInterval; //how often to minimize 
    double saddleDynamicsRecordIntervalInput;
    double saddleDynamicsRecordInterval;
    bool   saddleDynamicsLinearInterpolation;
    double saddleDynamicsMaxInitCurvature;

    bool   saddleConfinePositive; // undocumented
    bool   saddleBowlBreakout; // undocumented
    long   saddleBowlActive; // undocumented
    double saddleConfinePositiveMinForce; // undocumented
    double saddleConfinePositiveScaleRatio; // undocumented
    double saddleConfinePositiveBoost; // undocumented
    long   saddleConfinePositiveMinActive; // undocumented
    double saddleZeroModeAbortCurvature;

    // [Optimizer] //
    string optMethod;
    string optConvergenceMetric; // norm, max_atom, max_component
    string optConvergenceMetricLabel;
    long   optMaxIterations; // maximum iterations for saddle point searches and minimization
    double optMaxMove; // maximum displacement vector for a step during optimization
    double optConvergedForce; // force convergence criterion required for an optimization
    double optTimeStepInput;
    double optTimeStep; // time step size used in quickmin
    double optMaxTimeStepInput; // maximum time step for FIRE.
    double optMaxTimeStep; // maximum time step for FIRE.
    long   optLBFGSMemory; // number of previous forces to keep in the bfgs memory
    double optLBFGSInverseCurvature;
    double optLBFGSMaxInverseCurvature;
    bool   optLBFGSAutoScale;
    bool   optLBFGSAngleReset;
    bool   optLBFGSDistanceReset;
    bool   optQMSteepestDecent; // if set the velocity will always be set to zero in quickmin
    bool   optCGNoOvershooting; // if set it is ensured that the approximate line search in conjugate gradients never overshoot the minimum along the search line
    bool   optCGKnockOutMaxMove; // if set the old search direction is nullified when steps larger than the optMaxMove are conducted 
    bool   optCGLineSearch; // if set full line search is conducted
    double optCGLineConverged; // convergence criteria for line search, ratio between force component along search line and the orthogonal part
    long   optCGLineSearchMaxIter; // maximal nr of iterations during line search
    long   optCGMaxIterBeforeReset; // max nr of cg steps before reset, if 0 no resetting is done
    double optSDAlpha;
    bool   optSDTwoPoint;
    
    // [Dimer] //
    double dimerRotationAngle; // finite difference rotation angle
    bool   dimerImproved; // turn on the improved dimer method
    double dimerConvergedAngle; // stop rotating when angle drops below this value
    long   dimerMaxIterations; // maximum number of rotation iterations
    string dimerOptMethod; // method to determine the next rotation direction
    long   dimerRotationsMax; // old
    long   dimerRotationsMin; // old
    double dimerTorqueMax; // old
    double dimerTorqueMin; // old
    bool   dimerRemoveRotation; // remove dominant rotational component when estimating the eigenmode

    // [Lanczos] //
    double lanczosTolerance; // difference between the lowest eignevalues of two successive iterations
    long   lanczosMaxIterations; // maximum number of iterations
    bool   lanczosQuitEarly;

    // [Prefactor] //
    double prefactorDefaultValue; // default prefactor; calculate explicitly if zero
    double prefactorMaxValue; // max prefactor allowed
    double prefactorMinValue; // min prefactor allowed
    double prefactorWithinRadius; // atoms within this radius of the displaced atoms are put in the Hessian, unless filterMode is fraction
    double prefactorMinDisplacement;// atoms with displacement between min1 or min2 and the saddle point are put in the Hessian
    string prefactorRate;// method to estimate prefactor
    string prefactorConfiguration;// configuration for which the frequencies should be determined
    bool   prefactorAllFreeAtoms;// use all free atom when determining the prefactor
    string prefactorFilterScheme; // "cutoff" or "fraction", which use prefactorMinDisplacement or prefactorFilterFraction, respectively.
    double prefactorFilterFraction; // Include atoms whose summed motion comprise more than prefactorFilterFraction in the prefactor calculation. Prioritizes atoms that move more.
      
    // [Hessian] //
    string hessianAtomList;
    double hessianZeroFreqValue;

    // [Nudged Elastic Band] //
    long   nebImages;
    long   nebMaxIterations;
    double nebSpring;
    bool   nebClimbingImageMethod;
    bool   nebClimbingImageConvergedOnly;
    bool   nebOldTangent;
    bool   nebDoublyNudged;
    bool   nebDoublyNudgedSwitching;
    string nebOptMethod;
    bool   nebElasticBand;
    double nebConvergedForce; // force convergence criterion required for an optimization

    // [Molecular Dynamics] //
    double mdTimeStepInput;
    double mdTimeStep;
    double mdTimeInput;  
    double mdTime;  
    long   mdSteps;

    // [Parallel Replica] //
    bool   parrepRefineTransition;
    bool   parrepAutoStop;
    bool   parrepDephaseLoopStop;
    double parrepDephaseTimeInput;
    double parrepDephaseTime;
    long   parrepDephaseLoopMax;
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
    //std::vector<int> thermoAtoms;

    // [Replica Exchange] //
    string repexcTemperatureDistribution;
    long   repexcReplicas;
    long   repexcExchangeTrials;
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
    long   basinHoppingSteps;
    long   basinHoppingQuenchingSteps;
    bool   basinHoppingSignificantStructure;
    bool   basinHoppingSingleAtomDisplace;
    string basinHoppingDisplacementAlgorithm;
    string basinHoppingDisplacementDistribution;
    double basinHoppingSwapProbability;
    long   basinHoppingJumpMax;
    long   basinHoppingJumpSteps;
    bool   basinHoppingAdjustDisplacement;
    long   basinHoppingAdjustPeriod;
    double basinHoppingAdjustFraction;
    double basinHoppingTargetRatio;
    bool   basinHoppingWriteUnique;
    double basinHoppingStopEnergy;

    // [Global Optimization] //
    string globalOptimizationMoveMethod;
    string globalOptimizationDecisionMethod;
    long   globalOptimizationSteps;
    double globalOptimizationBeta;
    double globalOptimizationAlpha;
    long   globalOptimizationMdmin;
    double globalOptimizationTargetEnergy;

    // [Monte Carlo] //
    double monteCarloStepSize;
    int    monteCarloSteps;

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
    bool   writeMovies;
    long   writeMoviesInterval;

private:

    string toLowerCase(string s);

};

#endif

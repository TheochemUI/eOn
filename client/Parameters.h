//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <ctype.h>

#ifdef EONMPI
    #include "mpi.h"
#endif

using namespace std;

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

/** input parameters **/

    // [Main] //
    string job;
    long   randomSeed;
    double temperature;
    bool   quiet;
    bool   checkpoint;
    string iniFilename;
    string conFilename;
    double finiteDifference;
    long maxForceCalls;

    // [Potential] //
    string potential;
    double MPIPollPeriod;
    bool   LAMMPSLogging;
    bool   EMTRasmussen;
    bool   LogPotential;

    // [Structure Comparison] //
    double distanceDifference; // The distance criterion for comparing geometries
    double neighborCutoff; // radius used in the local atomic structure analysis
    bool   checkRotation;
    bool   indistinguishableAtoms;

	double energyDifference;

    // [Process Search] //
    bool   processSearchMinimizeFirst;
    double processSearchMinimizationOffset; // how far from the saddle to displace the minimization images

    // [Saddle Search]
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

    bool   saddleConfinePositive; // undocumented
    double saddleConfinePositiveMinForce; // undocumented
    double saddleConfinePositiveScaleRatio; // undocumented
    double saddleConfinePositiveBoost; // undocumented
    long   saddleConfinePositiveMinActive; // undocumented

    // [Optimizer]
    string optMethod;
    long   optMaxIterations; // maximum iterations for saddle point searches and minimization
    double optMaxMove; // maximum displacement vector for a step during optimization
    double optConvergedForce; // force convergence criterion required for an optimization
    double optTimeStep; // time step size used in quickmin
    double optMaxTimeStep; // maximum time step for FIRE.
    long   optLBFGSMemory; // number of previous forces to keep in the bfgs memory
    double optLBFGSInverseCurvature;
    bool   optQMSteepestDecent; // if set the velocity will always be set to zero in quickmin

    // [Dimer]
    double dimerRotationAngle; // finite difference rotation angle
    bool   dimerImproved; // turn on the improved dimer method
    double dimerConvergedAngle; // stop rotating when angle drops below this value
    long   dimerMaxIterations; // maximum number of rotation iterations
    string dimerOptMethod; // method to determine the next rotation direction
    long   dimerRotationsMax; // old
    long   dimerRotationsMin; // old
    double dimerTorqueMax; // old
    double dimerTorqueMin; // old

    // [Lanczos]
    double lanczosTolerance; // difference between the lowest eignevalues of two successive iterations
    long   lanczosMaxIterations; // maximum number of iterations

    // [Prefactor] //
    double prefactorDefaultValue; // default prefactor; calculate explicitly if zero
    double prefactorMaxValue; // max prefactor allowed
    double prefactorMinValue; // min prefactor allowed
    double prefactorWithinRadius; // atoms within this radius of the displaced atoms are put in the Hessian
    double prefactorMinDisplacement;// atoms with displacement between min1 or min2 and the saddle point are put in the Hessian
    string prefactorRate;// method to estimate prefactor
    string prefactorConfiguration;// configuration for which the frequencies should be determined
    bool   prefactorAllFreeAtoms;// use all free atom when determining the prefactor
      
    // [Hessian]
    string hessianAtomList;
    double hessianZeroFreqValue;

    // [Nudged Elastic Band]
    long   nebImages;
    long   nebMaxIterations;
    double nebSpring;
    bool   nebClimbingImageMethod;
    bool   nebOldTangent;
    string nebOptMethod;

    // [Molecular Dynamics]
    double mdTimeStepInput;
    double mdTimeStep;
    double mdTime;  
    long   mdSteps;

    // [Parallel Replica]
    bool   parrepRefineTransition;
    bool   parrepAutoStop;
    bool   parrepDephaseLoopStop;
    double parrepDephaseTime;
    long   parrepDephaseLoopMax;
    double parrepStateCheckInterval;
    double parrepRecordInterval;
    double parrepCorrTime;

    // [Thermostat]
    string thermostat;
    double thermoAndersenAlpha;
    double thermoAndersenTcol;
    double thermoNoseMass;
    double thermoLangvinFriction;

    // [Replica Exchange]
    string repexcTemperatureDistribution;
    long   repexcReplicas;
    long   repexcExchangeTrials;
    double repexcSamplingTime;
    double repexcTemperatureHigh;
    double repexcTemperatureLow;
    double repexcExchangePeriod;

    // [Bond Boost]
    string biasPotential;
    string bondBoostBALS;
    double bondBoostRMDTime;
    double bondBoostDVMAX;
    double bondBoostQRR;
    double bondBoostPRR;
    double bondBoostQcut;

    // [Basin Hopping]
    double basinHoppingMaxDisplacement;
    long   basinHoppingSteps;
    long   basinHoppingQuenchingSteps;
    bool   basinHoppingSignificantStructure;
    bool   basinHoppingSingleAtomDisplace;
    string basinHoppingMaxDisplacementAlgorithm;
    string basinHoppingDisplacementDistribution;
    double basinHoppingSwapProbability;
    long   basinHoppingJumpMax;
    long   basinHoppingJumpSteps;
    bool   basinHoppingInitialMD;
    double basinHoppingInitialMDTemperature;
    long   basinHoppingAdjustPeriod;
    double basinHoppingAdjustFraction;
    double basinHoppingTargetRatio;

	// [Global Optimization]
	string globalOptimizationMoveMethod;
	string globalOptimizationDecisionMethod;
	long globalOptimizationSteps;
	double globalOptimizationBeta;
	double globalOptimizationAlpha;
	long globalOptimizationMdmin;

    // MPI stuff, not actually specified in config file
    // this is used to pass information to the GPAW MPI
    // potential.
    int MPIPotentialRank;
    #ifdef EONMPI
        MPI_Comm MPIClientComm;
    #endif

    // [Debug]
    long   boincProgressMax;
    bool   writeMovies;
    long   writeMoviesInterval;
private:
    string toLowerCase(string s);
    

};

#endif

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

    // [Potential] //
    string potential;
    bool MPIPotentialAggressive;

    // [Structure Comparison] //
    double distanceDifference; // The distance criterion for comparing geometries
    double neighborCutoff; // radius used in the local atomic structure analysis
    bool checkRotation;

    // [Process Search] //
    bool   processSearchMinimizeFirst;
    double processSearchDefaultPrefactor; // default prefactor; calculate explicitly if zero
    double processSearchPrefactorMax; // max prefactor allowed
    double processSearchPrefactorMin; // min prefactor allowed
    double processSearchMinimizationOffset; // how far from the saddle to displace the minimization images

    // [Saddle Search]
//    bool   saddleDisplace; // do saddle search displacements client-side
    long   saddleMaxJumpAttempts; // number of initial displacements of the system to reach a convex region;  if 0, a search is started after the displacement
    long   saddleMaxIterations; // max iterations for saddle point searches and minimization
    string saddleMinmodeMethod; // algorithm to be used for lowest eigenmode determination
    string saddleDisplaceType; // displacement type to use
    double saddleMaxStepSize; // maximum length of the norm of the displacement when positive eigenvalue
    double saddleMaxEnergy; // energy above product state that will cause termination of the saddle point search
    double saddleDisplaceMagnitude; // norm of the displacement vector
    double saddleMaxSingleDisplace; // maximum value of displacement in x, y and z direction for atoms being displaced
    double saddleDisplaceRadius; // atoms within this radius of the displacement atoms are also displaced
    double saddleConvergedForce; // force convergence criterion required for a saddle point search
    double saddlePerpForceRatio; // proportion to keep of the perpendicular force when the lowest eigenvalue is positive 

    bool saddleConfinePositive; // undocumented
    double saddleConfinePositiveMinMove; // undocumented
    double saddleConfinePositiveScaleRatio; // undocumented
    long saddleConfinePositiveMaxActiveAtoms; // undocumented

    // [Optimizers]
    string optMethod;
    long   optMaxIterations; // maximum iterations for saddle point searches and minimization
    double optMaxMove; // maximum displacement vector for a step during optimization
    double optConvergedForce; // force convergence criterion required for an optimization
    double optTimeStep; // time step size used in quickmin
    bool optVariableTimeStep; // if quickmin time step should be adjusted
    long optLBFGSMemory; // number of previous forces to keep in the bfgs memory

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

    // [Hessian]
    string hessianType;
    double hessianMinDisplacement; // atoms with displacement between min1 or min2 and the saddle point are put in the Hessian
    double hessianWithinRadius; // atoms within this radius of the displaced atoms are put in the Hessian

    // [Nudged Elastic Band]
    long   nebImages;
    double nebSpring;
    bool   nebClimbingImageMethod;
    bool   nebOldTangent;
    string nebOptMethod;

    // [Molecular Dynamics]
    double mdTimeStep;
    double mdMaxMovedDist;
    bool   mdRefine;
    bool   mdAutoStop;
    long   mdRecordAccuracy;
    long   mdRefineAccuracy;
    long   mdSteps;
    long   mdDephaseSteps;
    long   mdCheckFreq;
    long   mdRelaxSteps;
    bool   mdDephaseLoopStop;
    long   mdDephaseLoopMax;

    // [Thermostat]
    string thermostat;
    double thermoAndersenAlpha;
    double thermoAndersenTcol;
    double thermoNoseMass;
    double thermoLangvinFriction;

    // [Thermostat]
    long drBalanceSteps;
    long drSamplingSteps;
    double drTargetTemperature;

    // [Bond Boost]
    string biasPotential;
    long   bondBoostRMDS;
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
    bool   basinHoppingMDFirst;
    double basinHoppingMDTemp;

    // MPI stuff, not actually specified in config file
    // this is used to pass information to the GPAW MPI
    // potential.
    int MPIPotentialRank;

    // [Debug]
    bool   writeMovies;
    long   writeMoviesSteps;
private:
    string toLowerCase(string s);

};

#endif

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

#include "Constants.h"

/** Contains all runtime parameters and results. No functionality just bookkeeping.*/
class Parameters {

public:
    Parameters();
    ~Parameters();
    int load(string filename);
    int load(FILE *file);

/** input parameters.*/

    // [Main] //
    long   job;
    long   randomSeed;
    string potential;
    double temperature;

    // [Structure Comparison] //
    double distanceDifference; // The distance criterion for comparing geometries
    double neighborCutoff; // radius used in the local atomic structure analysis
    double structureComparisonEnergyDifference;
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
    long   saddleMinmodeMethod; // algorithm to be used for lowest eigenmode determination
    long   saddleDisplaceType; // displacement type to use
    double saddleMaxStepSize; // maximum length of the norm of the displacement when positive eigenvalue
    double saddleMaxEnergy; // energy above product state that will cause termination of the saddle point search
    double saddleDisplaceMagnitude; // norm of the displacement vector
    double saddleMaxSingleDisplace; // maximum value of displacement in x, y and z direction for atoms being displaced
    double saddleDisplaceRadius; // atoms within this radius of the displacement atoms are also displaced
    double saddlePerpForceRatio; // proportion to keep of the perpendicular force when the lowest eigenvalue is positive 
    int    saddleMaxLocalizedAtoms; // maximum number of atoms to consider part of the mode

    // [Optimizers]
    string optMethod;
    long   optMaxIterations; // maximum iterations for saddle point searches and minimization
    double optMaxMove; // maximum displacement vector for a step during optimization
    double optConvergedForce; // force convergence criterion required for an optimization
    double optFiniteDist; // finite difference step size used in conjugate gradients
    double optTimeStep; // time step size used in quickmin

    // [Dimer]
    double dimerSeparation; // distance between the two dimer images
    double dimerRotationAngle; // finite difference rotation angle
    bool   dimerImproved; // turn on the improved dimer method
    double dimerConvergedRotation; // stop rotating when angle drops below this value
    long   dimerMaxIterations; // maximum number of rotation iterations
    long   dimerOptimizer; // method to determine the next rotation direction
    // old parameters
    long   dimerRotationsMax;
    long   dimerRotationsMin;
    double dimerTorqueMax;
    double dimerTorqueMin;

    // [Lanczos]
    double lanczosFiniteDist ; // finite difference distance
    double lanczosTolerance; // difference between the lowest eignevalues of two successive iterations
    long   lanczosMaxIterations; // maximum number of iterations

    // [Hessian]
    int    hessianType;
    double hessianFiniteDist; // finite difference distance
    double hessianMinDisplacement; // atoms with displacement between min1 or min2 and the saddle point are put in the Hessian
    double hessianWithinRadius; // atoms within this radius of the displaced atoms are put in the Hessian

    // [Displacement Sampling] 
    long   displaceNSamples;
    long   displaceIterMax; // maximum number of rotations to perform on the dimer
    double displaceTorqueConvergence; // convergence criteria of the dimer rotation
    double displaceMaxCurvature; // maximum curvature which is considered good; used to avoid soft negative curvatures
    double displaceMaxDE; // maximum dE for which a sample is considered good
    string displaceCutoffs;
    string displaceMagnitudes;

    // [Nudged Elastic Band]
    long   nebImages;
    double nebSpring;
    bool   nebClimb;
    string nebOptMethod;
    long   nebOptMaxIterations; // maximum iterations for saddle point searches and minimization
    double nebOptMaxMove; // maximum displacement vector for a step during optimization
    double nebOptConvergedForce; // force convergence criterion required for an optimization
    double nebOptFiniteDist; // finite difference step size used in conjugate gradients
    double nebOptTimeStep; // time step size used in quickmin

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
    long   thermostat;
    double thermoAndersenAlpha;
    double thermoAndersenTcol;
    double thermoNoseMass;
    double thermoLangvinFriction;

    // [Thermostat]
    long drBalanceSteps;
    long drSamplingSteps;
    double drTargetTemperature;

    // [Bond Boost]
    long   biasPotential; 
    long   bondBoostRMDS;
    double bondBoostDVMAX;
    double bondBoostQRR; 
    double bondBoostPRR; 
    double bondBoostQcut;

    // [Basin Hopping]
    double basinHoppingStepSize;
    long   basinHoppingSteps;
    bool   basinHoppingStayMinimized;
    bool   basinHoppingSingleAtomDisplace;

    // [Debug]
    bool   writeMovies;

private:
    string toLowerCase(string s);

};
#endif

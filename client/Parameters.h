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
    long   job; // indicate which type of job will be run
    long   randomSeed; // seed for random generator
    long   potential; // tag to describe which potential to use
    double temperature;

    // [Structure Comparison] //
    double distanceDifference; // The distance criterion for comparing geometries
    double neighborCutoff; // radius used in the local atomic structure analysis
    double structureComparisonEnergyDifference;

    // [Process Search] //
    bool   processSearchMinimizeFirst;
    double processSearchDefaultPrefactor; // default prefactor; calculate explicitly if zero
    double processSearchPrefactorMax; // max prefactor allowed
    double processSearchPrefactorMin; // min prefactor allowed


    // [Saddle Search]
    bool   saddleDisplace; // do saddle search displacements client-side
    long   saddleMaxJumpAttempts; // number of initial displacements to bring the system directly to a convex region;  if 0, a search is always started after the displacement
    long   saddleMaxIterations; // max iterations for saddle point searches and minimization
    long   saddleMinmodeMethod; // the algorithm to be used for lowest eigenmode determination
    long   saddleDisplaceType; // displacement type to use
    double saddleMaxStepSize; // Max length of the norm of the displacement when positive eigenvalue
    double saddleMaxEnergy; // Energy above product state that will cause termination of the saddle point search
    double saddleDisplaceMagnitude; // The norm of the displacement vector
    double saddleMaxSingleDisplace; // max value of displacement in x, y and z direction for atoms being displaced
    double saddleDisplaceRadius; // Atoms within this radius of the displacement atoms are also displaced
    double saddlePerpForceRatio; // proportion to keep of the perpendicular force when the lowest eigenvalue is positive 

    // [Optimizers]
    long   optMaxIterations; // max iterations for saddle point searches and minimization
    double optMaxMove; // maximum displacement vector for a step during optimization
    double optConvergedForce; // force convergence criterion required for an optimization
    double optFiniteDiffStep; // finite difference step size used in conjugate gradients
    double optTimeStep; // time step size used in quickmin

    // [Dimer]
    double dimerSeparation; // distance between the two dimer images
    double dimerRotationAngle; // finite difference rotation angle
    double dimerMaxIterations;
    long   dimerRotationsMax;
    long   dimerRotationsMin;
    double dimerWindowMax;
    double dimerWindowMin;
   
    // [Lanczos]
    double lanczosConvergence; // difference between the lowest eignevalues of two successive iterations
    int    lanczosIteration; // maximum number of iterations

    // [Hessian]
    int    hessianType;
    double hessianMinDisplacement; // atomic displacement between min1 and the saddle point or min2 and the saddle point causing the atom to be accounted for in the Hessian
    double hessianWithinRadius; // atoms within this radius of one the atom considered displace are also accounted for in the Hessian
    double hessianPrefactorMax; // max prefactor allowed
    double hessianPrefactorMin; // min prefactor allowed

    // [Displacement Sampling] 
    long   displaceNSamples;
    long   displaceIterMax; // maximum number of rotations to perform on the dimer
    double displaceTorqueConvergence; // convergence criteria of the dimer rotation
    double displaceMaxCurvature; // maximum curvature for which a sample is considered good; used to avoid shallow but negative curvatures
    double displaceMaxDE; // maximum dE for which a sample is considered good
    string displaceCutoffs;
    string displaceMagnitudes;

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

    // [Bond Boost]
    long   biasPotential; 
    long   bondBoostRMDS;
    double bondBoostDVMAX;
    double bondBoostQRR; 
    double bondBoostPRR; 
    double bondBoostQcut;

    // Basin Hopping
    double basinHoppingStepSize;
    long   basinHoppingSteps;
 
    // [Debug] //
    bool   saveStdout;

private:
    string toLowerCase(string s);

};
#endif

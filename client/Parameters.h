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

    long   job; // indicate which type of job will be run
    long   randomSeed; // seed for random generator
    long   potential; // tag to describe which potential to use
    bool   minimizeOnly; // only perform minimization, not saddle search
    bool   minimizeBox; // also minimize the box dimensions if minimize_only_ is true
    bool   calculatePrefactor; // flag to determine if prefactors should be determined. 
    double convergedForceRelax; // force criterion for minimization convergence
    long   maximumIterations; // max iterations for saddle point searches and minimization

    bool   saveStdout;

    double cgCurvatureStep; // finite difference step size used in conjugate gradients
    double cgMaxMoveFullRelax; // maximum displacement vector for a step during minimization
    double qmTimeStep; // time step size used in Quickmin
    double qmMaxMove; // time step size used in Quickmin

    double maxDifferencePos; // The distance criterion for comparing geometries
    double neighborCutoff; // radius used in the local atomic structure analysis

    bool   processSearchMinimizeFirst;

    bool   saddleDisplace; // do saddle search displacements client-side
    long   saddleMaxJumpAttempts; // number of initial displacements to bring the system directly to a convex region;  if 0, a search is always started after the displacement
    long   saddleMaxIterations; // max iterations for saddle point searches and minimization
    long   saddleMinmodeMethod; // the algorithm to be used for lowest eigenmode determination
    long   saddleDisplaceType; // displacement type to use
    double saddleConvergedForce; // force convergence criterion required for a saddle point search
    double saddleMaxStepSize; // Max length of the norm of the displacement when positive eigenvalue
    double saddleMaxEnergy; // Energy above product state that will cause termination of the saddle point search
    double saddleNormDisplace; // The norm of the displacement vector
    double saddleMaxSingleDisplace; // max value of displacement in x, y and z direction for atoms being displaced
    double saddleWithinRadiusDisplace; // Atoms within this radius this of the one defining the center of the displacement are also being displaced with the value sizePerturbation_SP_ 
    double saddlePerpForceRatio; // proportion to keep of the perpendicular force when the lowest eigenvalue is positive 

    int    hessianType;
    double hessianMinDisplacement; // atomic displacement between min1 and the saddle point or min2 and the saddle point causing the atom to be accounted for in the Hessian
    double hessianWithinRadius; // atoms within this radius of one the atom considered displace are also accounted for in the Hessian
    double hessianPrefactorMax; // max prefactor allowed
    double hessianPrefactorMin; // min prefactor allowed

    long   dimerRotationsMax;
    long   dimerRotationsMin;
    double dimerWindowMax;
    double dimerWindowMin;
    //long dimerRotationsNewSearch; // number of iteration before starting a new saddle point search used in dimer
    double dimerSeparation; // distance between the two dimer images
    double dimerRotationAngle; // finite difference rotation angle
    double dimerMaxIterations;
    
    long   displaceNSamples;
    long   displaceIterMax; // maximum number of rotations to perform on the dimer
    double displaceTorqueConvergence; // convergence criteria of the dimer rotation
    double displaceMaxCurvature; // maximum curvature for which a sample is considered good; used to avoid shallow but negative curvatures
    double displaceMaxDE; // maximum dE for which a sample is considered good
    string displaceCutoffs;
    string displaceMagnitudes;
    
    double lanczosConvergence; // difference between the lowest eignevalues of two successive iterations
    int    lanczosIteration; // maximum number of iterations

    double mdTimeStep;
    double mdTemperature;
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

    long   thermostat;
    double thermoAndersenAlpha;
    double thermoAndersenTcol;
    double thermoNoseMass;

    long   biasPotential; 
    long   bondBoostRMDS;
    double bondBoostDVMAX;
    double bondBoostQRR; 
    double bondBoostPRR; 
    double bondBoostQcut;
 
private:
    string toLowerCase(string s);

};
#endif

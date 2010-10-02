/*
 *===============================================
 *  Parameters.h
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 10/30/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */

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

/** All input parameters. If one would use values different from the default values these should be specified in the file with the name set in (Constants::PARMS_FILE_NAME).*/
struct Input {
    long randomSeed; // seed for random generator
    long reactantStateTag; // tag to describe to which reactant state the saddle point connects
    long potentialTag; // tag to describe which potential to use. Compare with values in Constants.cpp
    long potentialNoTranslation; // translation will be removed, handled in the potential class   
    long minimizeOnly; // only perform minimization, not saddle search
    long minimizeBox; // also minimize the box dimensions if minimize_only_ is true
    long getPrefactorsTag; // tag to describe if the prefactors should be determined. 
    double convergedRelax; // converge criterion during relaxation [eV/A]
    long maximumIterations; // max iterations for saddle point searches and minimization

    double cgCurvatureStep; // finite difference step size used in conjugate gradients
    double cgMaxMoveFullRelax; // maximum displacement vector for a step during minimization
    double qmTimeStep; // time step size used in Quickmin

    double maxDifferencePos; // The distance criterion for comparing geometries
    double neighborCutoff; // radius used in the local atomic structure analysis

    bool saddleRefine; // refine saddle point
    long saddleMaxJumpAttempts; // how many times the initial displacement should try to bring the system directly to convex region. If 0 a search is started after the displacement no matter what
    long saddleNrOfTriesToDetermine; // max number of attempts to obtain a converged saddle point
    long saddleMaxIterations; // max iterations for saddle point searches and minimization [GH: fix comment]
    long saddleMaxIterationsConcave; // max iterations for saddle point searches and minimization [GH: fix comment]
    long saddleLowestEigenmodeDetermination; // the algorithm to be used for lowest eigenmode determination; compare with values in Constants.cpp (now in SaddlePoint.h)
    long saddleTypePerturbation; // displacement type to use; compare with values in Constants.cpp (now in SaddlePoint.h)
    double saddleConverged; // converge criterion during saddle point search [eV/A]
    double saddleMaxStepSize; // Max length of the norm of the displacement when positive eigenvalue [A]
    double saddleMaxEnergy; // Energy above product state that will cause termination of the saddle point search [eV]
    double saddleNormPerturbation; // The norm or the perturbation vector [A]
    double saddleMaxSinglePerturbation; // max value of displacement in x, y and z direction for atoms being perturbated [A]
    double saddleWithinRadiusPerturbated; // Atoms within this radius this of the one defining the center of the displacement are also being dispalced with the value sizePerturbation_SP_ [A]
	double saddlePerpendicularForceRatio; // proportion to keep of the perpendicular force when the lowest eigenvalue is positive 

    long hessianMaxSize; // if specified, the size of the hessian determined will be equal or smaller than this value
	double hessianMinDisplacement; // atomic displacement between min1 and the saddle point or min2 and the saddle point causing the atom to be accounted for in the Hessian [A]
    double hessianWithinRadiusDisplaced; // atoms within this radius of one the atom considered displace are also accounted for in the Hessian [A]
    double hessianPrefactorMax; // max prefactor allowed
    double hessianPrefactorMin; // min prefactor allowed

    long dimerRotations; // number of rotation iterations during the eigenmode estimation used in dimer
    //long dimerRotationsNewSearch; // number of iteration before starting a new saddle point search used in dimer
    double dimerSeparation; // distance between the two dimer images
    double dimerRotationAngle; // finite difference rotation angle
    double dimerConvergenceLimit;
    double dimerMaxIterations;

    double lanczosMaxIterations;
    double lanczosConvergenceLimit;
};

/** All results being obtained. The results are store in the file with the name set in (Constants::RESULTS_FILE_NAME).*/
struct Output {
    long terminationReason; // tag to indicate if the calculation converged. Compare with values in Constants.h
    long forceCalls_; // total number of force calls
    long forceCallsSaddlePointConcave_; // number of force calls used during the saddle point determination in the concave
    long forceCallsSaddlePointConvex_; // number of force calls used during the saddle point determination in the convex region
    long forceCallsPrefactors_; // number of force calls used during the prefactor determination
    double potentialEnergySP; // energy of the saddle point [eV]
    double potentialEnergyMin1; // energy of min1 [eV]
    double potentialEnergyMin2; // energy of min2 [eV]
    double barrierReac_Prod; // barrier (dE[1->2]) for process from min1 -> min2 [eV]
    double barrierProd_Reac; // barrier (dE[2->1]) for process from min2 -> min1 [eV]
    double prefactorReac_Prod; // prefactor (v[1->2]) for process from min1 -> min2 [1/s]
    double prefactorProd_Reac; // prefactor (v[2->1]) for process from min2 -> min1 [1/s]
    double displacementSaddleDistance; // distance between the displacement and the discovered saddle
};

/** Contains all runtime parameters and results. No functionality just bookkeeping.*/
class Parameters : public Input, public Output{
public:
    Parameters();
    ~Parameters();
    void load(string filename);
    void load(FILE *file);
    void saveOutput(FILE *file);
    void saveInput(FILE *file);
    void printInput();
    void printOutput();

    enum JobType {PROCESS_SEARCH, SADDLE_SEARCH, MINIMIZATION, UNKNOWN_JOBTYPE};
    JobType jobType;
    
    // Accounting for the force calls
    long getForceCalls();
    long getForceCallsSaddlePoint();
    long getForceCallsSaddlePointConcave();
    long getForceCallsSaddlePointConvex();
    long getForceCallsPrefactors();
    void addForceCalls(long forceCalls);
    void addForceCallsSaddlePoint(long forceCallsSaddlePoint, double eigenvalue);
    void addForceCallsPrefactors(long forceCallsPrefactors);
    void resetForceCalls();
    void resetForceCallsSaddlePoint();
    void resetForceCallsPrefactors();
};
#endif

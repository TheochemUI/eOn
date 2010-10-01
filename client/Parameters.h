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
    long randomSeed_; // seed for random generator
    long reactantStateTag_; // tag to describe to which reactant state the saddle point connects
    long potentialTag_; // tag to describe which potential to use. Compare with values in Constants.cpp
    long potentialNoTranslation_; // translation will be removed, handled in the potential class   
    long minimize_only_; // only perform minimization, not saddle search
    long minimize_box_; // also minimize the box dimensions if minimize_only_ is true
    long getPrefactorsTag_; // tag to describe if the prefactors should be determined. 
    long typePerturbation_SP_; // displacement type to use; compare with values in Constants.cpp (now in SaddlePoint.h)
    bool refine_SP_; // refine saddle point
    long lowestEigenmodeDetermination_SP_; // the algorithm to be used for lowest eigenmode determination; compare with values in Constants.cpp (now in SaddlePoint.h)
    double converged_Relax_; // converge criterion during relaxation [eV/A]

    double cgCurvatureStep_; // finite difference step size used in conjugate gradients
    double cgMaxMoveFullRelax_; // maximum displacement vector for a step during minimization
    double qmTimeStep_; // time step size used in Quickmin
    double maxDifferencePos_; // The distance criterion for comparing geometries
    double neighborCutoff_; // radius used in the local atomic structure analysis

    double converged_SP_; // converge criterion during saddle point search [eV/A]
    long maxJumpAttempts_SP_; // how many times the initial displacement should try to bring the system directly to convex region. If 0 a search is started after the displacement no matter what
    long nrOfTriesToDetermineSaddlePoint_SP_; // max number of attempts to obtain a converged saddle point
    double maxStepSize_SP_;///< Max length of the norm of the displacement when positive eigenvalue [A]
//    double maxStepSizeConcave_SP_; // max length of the norm of the displacement when positive eigenvalue [A]
//    double maxStepSizeConvex_SP_; // max length of the norm of the displacement when negative eigenvalue [A]
    double maxEnergy_SP_;///< Energy above product state that will cause termination of the saddle point search [eV]
    double normPerturbation_SP_;///< The norm or the perturbation vector [A]
    double maxSinglePerturbation_SP_; // max value of displacement in x, y and z direction for atoms being perturbated [A]
    double withinRadiusPerturbated_SP_; // Atoms within this radius this of the one defining the center of the displacement are also being dispalced with the value sizePerturbation_SP_ [A]
    long maximumIterations_; // max iterations for saddle point searches and minimization

    long maxIterations_; // maximum number of iterations in a saddle search
    long maxIterationsConcave_; // maximum number of iterations in the concave region

	double perpendicularForceRatio_; // proportion to keep of the perpendicular force when the lowest eigenvalue is positive 
    long maxSize_Hessian_; // if specified, the size of the hessian determined will be equal or smaller than this value
	double minDisplacement_Hessian_; // atomic displacement between min1 and the saddle point or min2 and the saddle point causing the atom to be accounted for in the Hessian [A]
    double withinRadiusDisplaced_Hessian_; // atoms within this radius of one the atom considered displace are also accounted for in the Hessian [A]

    double prefactorMax_; // max prefactor allowed
    double prefactorMin_; // min prefactor allowed

    long rotations_Dimer_; // number of rotation iterations during the eigenmode estimation used in Dimer
    //long rotationsNewSearch_Dimer_; // number of iteration before starting a new saddle point search used in Dimer
    double separation_Dimer_; // distance between the two dimer images
    double rotationAngle_Dimer_; // finite difference rotation angle
    double convergenceLimit_Lanczos_;
    double iterationLimit_Lanczos_;
};

/** All results being obtained. The results are store in the file with the name set in (Constants::RESULTS_FILE_NAME).*/
struct Output {
    long terminationReason_; // tag to indicate if the calculation converged. Compare with values in Constants.h
    long forceCalls_; // total number of force calls
    long forceCallsSaddlePointConcave_; // number of force calls used during the saddle point determination in the concave
    long forceCallsSaddlePointConvex_; // number of force calls used during the saddle point determination in the convex region
    long forceCallsPrefactors_; // number of force calls used during the prefactor determination
    double potentialEnergySP_; // energy of the saddle point [eV]
    double potentialEnergyMin1_; // energy of min1 [eV]
    double potentialEnergyMin2_; // energy of min2 [eV]
    double barrierReac_Prod_; // barrier (dE[1->2]) for process from min1 -> min2 [eV]
    double barrierProd_Reac_; // barrier (dE[2->1]) for process from min2 -> min1 [eV]
    double prefactorReac_Prod_; // prefactor (v[1->2]) for process from min1 -> min2 [1/s]
    double prefactorProd_Reac_; // prefactor (v[2->1]) for process from min2 -> min1 [1/s]
    double displacement_saddle_distance_; // distance between the displacement and the discovered saddle
};

/** Contains all runtime parameters and results. No functionality just bookkeeping.*/
class Parameters : private Input, private Output{
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
    JobType job_Type_;
    
    // Passing the input parameters
    long getRandomSeed();
    void setRandomSeed(long randomSeed);
    long getReactantStateTag();
    long getPotentialTag();
    long getPotentialNoTranslation();
    bool getMinimizeOnly();
    bool getMinimizeBox();
    long getPrefactorsTag();
    long getTypePerturbation_SP();
    bool getRefineSP();
    long getLowestEigenmodeDetermination_SP();

    double getConverged_Relax();

    double getCgCurvatureStep();
    double getCgMaxMoveFullRelax();
    double getQmTimeStep();
    double getMaxDifferencePos();
    double getNeighborCutoff();

    double getConverged_SP();
    long getMaxJumpAttempts_SP();
    long getNrOfTriesToDetermineSaddlePoints_SP();
    double getEnergyIntoMode_SP();
    double getMaxStepSize_SP();
//    double getMaxStepSizeConcave_SP();
//    double getMaxStepSizeConvex_SP();
    double getMaxEnergy_SP();
    double getNormPerturbation_SP();
    double getWithinRadiusPerturbated_SP();
    double getMaxSinglePerturbation_SP();
    long getMaximumIterations();
    long getMaxIterations();
    long getMaxIterationsConcave();

    double getPerpendicularForceRatio();

	long getMaxSize_Hessian();
    double getMinDisplacement_Hessian();
    double getWithinRadiusDisplaced_Hessian();
    double getPrefactorMax();
    double getPrefactorMin();


    double getSeparation_Dimer();
    double getRotationAngle_Dimer();

    long getRotations_Dimer();
    //long getRotationsNewSearch_Dimer();
    
    double getConvergenceLimit_Lanczos();
    double getIterationLimit_Lanczos();

    // Setting results values
    void setTerminationReason(long terminationReason);
    void setPotentialEnergySP(double potentialEnergySP);
    void setPotentialEnergyMin1(double potentialEnergyMin1);
    void setPotentialEnergyMin2(double potentialEnergyMin2);
    void setPrefactorReac_Prod(double prefactorReac_Prod);
    void setPrefactorProd_Reac(double prefactorProd_Reac);
    void setBarrierReac_Prod(double barrierReac_Prod);
    void setBarrierProd_Reac(double barrierProd_Reac);
    void setDisplacementSaddleDistance(double displacement_saddle_distance);

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
    
    long linesInFile(FILE *file); // Determines the number of lines in the specified file. The number of lines should be used to initialize a characther 2D array being passed as an arguement in loadParameters.
};
#endif

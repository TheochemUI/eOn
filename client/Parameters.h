/*
 *===============================================
 *  Parameters.h
 *  eon2
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

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <ctype.h>

#include "Constants.h"

/** All input parameters. If one would use values different from the default values these should be specified in the file with the name set in (Constants::PARMS_FILE_NAME).*/
struct Input {
    long randomSeed_;///< Seed for random generator.
    long reactantStateTag_;///< Tag to describe to which reactant state the saddle point connects.
    long potentialTag_;///< Tag to describe which potential to use. Compare with values in Constants.cpp
    long potentialNoTranslation_;///< Translation will be removed, handled in the potential class   
    long minimize_only_; /* RT: only perform minimization, not saddle search. */
    long minimize_box_; /* RT: also minimize the box dimensions if minimize_only_ is true. */
    long getPrefactorsTag_;///< Tag to describe if the prefactors should be determined.  
    long typePertubation_SP_;///< Displacement type to use. Compare with values in Constants.cpp
    bool refine_SP_;///< Refine saddle point.
    long lowestEigenmodeDetermination_SP_;///< The algorithm to be used for lowest eigenmode determination. Compare with values in Constants.cpp
    double converged_Relax_;///< Converge criterion during relaxation [eV/A].
    double converged_SP_;///< Converge criterion during saddle point search [eV/A].
    long maxJumpAttempts_SP_;///< How many times the initial displacement should try to bring the system directly to convex region. If 0 a search is started after the displacement no matter what.
    long nrOfTriesToDetermineSaddlePoint_SP_;///< Maximal number of attempts to obtain a converged saddle point;
    double maxStepSizeConcave_SP_;///< Max length of the norm of the displacement when positive eigenvalue [A].
    double maxStepSizeConvex_SP_;///< Max length of the norm of the displacement when negative eigenvalue [A].
    double maxEnergy_SP_;///< Energy above product state that will cause termination of the saddle point search [eV].
    double normPertubation_SP_;///< The norm or the pertubation vector [A].
    double maxSinglePertubation_SP_;///< Max value of displacement in x, y and z direction for atoms being pertubated [A].
    double withinRadiusPertubated_SP_;///< Atoms within this radius this of the one defining the center of the displacement are also being dispalced with the value sizePertubation_SP_ [A].
    long maximumIterations_;///< Maximum of iterations for saddle point searches and minimisation.
    double minDisplacement_Hessian_;///< Atomic displacement between min1 and the saddle point or min2 and the saddle point causing the atom to be accounted for in the Hessian [A].
    double withinRadiusDisplaced_Hessian_;///< Atoms within this radius of one the atom considered displace are also accounted for in the Hessian [A].
    long rotations_Dimer_;///< The number of rotation iterations during the eigenmode estimation used in Dimer.
    long rotationsNewSearch_Dimer_;///< The number of iteration before starting a new saddle point search used in Dimer.
    double dimer_dR_; /* RT: The distance between the two dimer images. */
    double lanczosConvergenceLimit_;
    double iterationLimit_;
};

/** All results being obtained. The results are store in the file with the name set in (Constants::RESULTS_FILE_NAME).*/
struct Output {
    long terminationReason_;///< Tag to indicate if the calculation converged. Compare with values in Constants.h
    long forceCalls_;///< The total number of force calls.
    long forceCallsSaddlePointConcave_;///< The total number of force calls used during the saddle point determination in the concave.
    long forceCallsSaddlePointConvex_;///< The total number of force calls used during the saddle point determination in the convex region.
    long forceCallsPrefactors_;///< The total number of force calls used during the prefactor determination.
    double potentialEnergySP_;///< Energy of the saddle point [eV].
    double potentialEnergyMin1_;///< Energy of min1 [eV].
    double potentialEnergyMin2_;///< Energy of min2 [eV].
    double barrierReac_Prod_;///< Barrier (dE[1->2]) for process from min1 -> min2 [eV].
    double barrierProd_Reac_;///< Barrier (dE[2->1]) for process from min2 -> min1 [eV].
    double prefactorReac_Prod_;///< Prefactor (v[1->2]) for process from min1 -> min2 [1/s].
    double prefactorProd_Reac_;///< Prefactor (v[2->1]) for process from min2 -> min1 [1/s].
    double displacement_saddle_distance_; // RT: The distance between the displacement and the discovered saddle.
};

/** Contains all runtime parameters and results. No functionality just bookkeeping.*/
class Parameters : private Input, private Output{
public:
    Parameters();
    ~Parameters();
    void load(FILE *file);
    void saveOutput(FILE *file);
    void saveInput(FILE *file);
    void printInput();
    void printOutput();
    
    // Passing the input parameters.
    long getRandomSeed();
    void setRandomSeed(long randomSeed);
    long getReactantStateTag();
    long getPotentialTag();
    long getPotentialNoTranslation();
    bool getMinimizeOnly();
    bool getMinimizeBox();
    long getPrefactorsTag();
    long getTypePertubation_SP();
    bool getRefineSP();
    long getLowestEigenmodeDetermination_SP();

    double getConverged_Relax();
    double getConverged_SP();
    long getMaxJumpAttempts_SP();
    long getNrOfTriesToDetermineSaddlePoints_SP();
    double getEnergyIntoMode_SP();
    double getMaxStepSizeConcave_SP();
    double getMaxStepSizeConvex_SP();
    double getMaxEnergy_SP();
    double getNormPertubation_SP();
    double getWithinRadiusPertubated_SP();
    double getMaxSinglePertubation_SP();
    long getMaximumIterations();
    double getMinDisplacement_Hessian();
    double getWithinRadiusDisplaced_Hessian();
    double getDimerDR();

    long getRotations_Dimer();
    long getRotationsNewSearch_Dimer();
    
    double getLanczosConvergenceLimit();
    double getLanczosIterationLimit();

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
    
private:        
    long linesInFile(FILE *file);///< Determines the number of lines in the specified file. The number of lines should be used to initialize a characther 2D array being passed as an arguement in loadParameters.
    void loadParameters(FILE *file, char **parms, double *values, long nLines);///< Load the first word in each line into parms. The next 'word' is assumed to be the values correspond to the parameter and is stored in values. Line 0: first word->parms[0] second 'word'->values[0]. Line 1 : first word->parms[1] second 'word'->values[1]. And so on. Lines stating with # are considered user comments and are ignored ('# comment')   
};
#endif

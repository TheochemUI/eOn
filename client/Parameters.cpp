/*
 *===============================================
 *  Created by Andreas Pedersen on 10/30/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#include "Parameters.h"
#include "INIFile.h"

Parameters::Parameters(){

    // Default values
    randomSeed = -1;
    reactantStateTag = 0;
    potentialTag = 1;
    potentialNoTranslation = 0;
    getPrefactorsTag = 0;
    minimizeOnly = 0;
    minimizeBox = 0;

    jobType = PROCESS_SEARCH;
    
    // Parameters, default values if not read in from parameters_passed.dat

    // Value used in the Relaxation   
    convergedRelax = 0.005;

    // Constants used in the client
    maxDifferencePos = 0.1; // The distance criterion for comparing geometries
    neighborCutoff = 0.33;

    // Values used in the Saddle Point determination   
    saddleTypePerturbation = 1;
    saddleRefine = false;
    saddleLowestEigenmodeDetermination = 1;
    saddleConverged = 0.025;
    saddleMaxJumpAttempts = 0;
    saddleNrOfTriesToDetermine = 1;
    saddleMaxStepSize = 0.2;
    saddleMaxEnergy = 20.0;
    saddleNormPerturbation = 0.1;
    saddleWithinRadiusPerturbated = 4.0;
    saddleMaxSinglePerturbation = 0.1;
    saddleMaxIterations = 512;
    saddleMaxIterations = 256;
    saddleMaxIterationsConcave = 256;
    saddlePerpendicularForceRatio = 0.0;

    // Values used in the Hessian determination   
	hessianMaxSize = 0;
	hessianMinDisplacement = 0.25;
    hessianWithinRadiusDisplaced = 5.0;
    hessianPrefactorMax = 10e20;
    hessianPrefactorMin = 10e8;

    // Values used in the Dimer method
    dimerRotations = 1;
    dimerSeparation = 0.0001;
    dimerRotationAngle = 0.005;
    
    // Defaults for Lanczos
    lanczosMaxIterations = 50;
    lanczosConvergenceLimit = 1e-4;

    // Initializing the cummulative output
    forceCalls_ = 0;
    forceCallsSaddlePointConcave_ = 0;
    forceCallsSaddlePointConvex_ = 0;
    forceCallsPrefactors_ = 0;
    displacementSaddleDistance = 0;

    // Parameters used by the optimizers
    cgCurvatureStep = 0.001; // finite difference step size used in conjugate gradients
    cgMaxMoveFullRelax = 0.2; // maximum displacement vector for a step during minimization
    qmTimeStep = 0.1; // time step size used in Quickmin.

    return;
}

Parameters::~Parameters(){
    return;
}

void Parameters::load(string filename)
{
    FILE *parametersFile;

    parametersFile = fopen(filename.c_str(), constants::READ.c_str());
    load(parametersFile);
    fclose(parametersFile);
}

void Parameters::load(FILE *file){
    
    //If we fail to parse the file as an INI, we will need to rewind the file. So, we store the current stream position.
    fpos_t pos;
    fgetpos(file, &pos); 

    CIniFile ini;
    ini.CaseInsensitive();
    if(ini.ReadFile(file))
    {
        //If we succesfully read the file, then parse it as an INI
        randomSeed = ini.GetValueL("Default", "RANDOM_SEED", randomSeed);
        reactantStateTag = ini.GetValueL("Default", "REACTANT_STATE_TAG", reactantStateTag);
        potentialTag = ini.GetValueL("Default", "POTENTIAL_TAG", potentialTag);
        potentialNoTranslation = ini.GetValueL("Default", "POTENTIAL_NO_TRANSLATION", potentialNoTranslation);
        minimizeOnly = ini.GetValueL("Default", "MINIMIZE_ONLY", minimizeOnly);
        minimizeBox = ini.GetValueL("Default", "MINIMIZE_BOX", minimizeBox);
        getPrefactorsTag = ini.GetValueL("Default", "GET_PREFACTORS_TAG", getPrefactorsTag);
        convergedRelax = ini.GetValueF("Default", "CONVERGED_RELAX", convergedRelax);
        maximumIterations = ini.GetValueL("Default", "MAXIMUM_ITERATIONS", maximumIterations);

        string jobTypeString;
        jobTypeString = ini.GetValue("Default", "JOB_TYPE", "ProcessSearch");

        if (jobTypeString == "ProcessSearch") {
            jobType = PROCESS_SEARCH;
        }else if (jobTypeString == "SaddleSearch") {
            jobType = SADDLE_SEARCH;
        }else if (jobTypeString == "Minimization") {
            jobType = MINIMIZATION;
        }else{
            jobType = UNKNOWN_JOBTYPE;
        }
        
        saddleTypePerturbation = ini.GetValueL("Saddle_Point", "TYPE_PERTURBATION", saddleTypePerturbation);
        saddleLowestEigenmodeDetermination = ini.GetValueL("Saddle_Point", "LOWEST_EIGENMODE_DETERMINATION", saddleLowestEigenmodeDetermination);
        saddleRefine = ini.GetValueB("Saddle_Point", "REFINE", saddleRefine);
        saddleConverged = ini.GetValueF("Saddle_Point", "CONVERGED", saddleConverged);
        saddleMaxJumpAttempts = ini.GetValueL("Saddle_Point", "MAX_JUMP_ATTEMPTS", saddleMaxJumpAttempts);
        saddleNrOfTriesToDetermine = ini.GetValueL("Saddle_Point", "NR_OF_TRIES_TO_DETERMINE", saddleNrOfTriesToDetermine);
        saddleMaxStepSize = ini.GetValueF("Saddle_Point", "MAX_STEP_SIZE", saddleMaxStepSize);
        saddleMaxEnergy = ini.GetValueF("Saddle_Point", "MAX_ENERGY", saddleMaxEnergy);
        saddleNormPerturbation = ini.GetValueF("Saddle_Point", "NORM_PERTURBATION", saddleNormPerturbation);
        saddleMaxSinglePerturbation = ini.GetValueF("Saddle_Point", "MAX_SINGLE_PERTURBATION", saddleMaxSinglePerturbation);
        saddleWithinRadiusPerturbated = ini.GetValueF("Saddle_Point", "WITHIN_RADIUS_PERTURBATED", saddleWithinRadiusPerturbated);
        saddleMaxIterations = ini.GetValueL("Saddle_Point", "MAX_ITERATIONS", saddleMaxIterations);
        saddlePerpendicularForceRatio = ini.GetValueL("Default", "PERPENDICULAR_FORCE_RATIO", saddlePerpendicularForceRatio);

        hessianMaxSize = ini.GetValueL("Hessian", "MAX_SIZE", hessianMaxSize);
        hessianWithinRadiusDisplaced = ini.GetValueF("Hessian", "WITHIN_RADIUS_DISPLACED", hessianWithinRadiusDisplaced);
        hessianMinDisplacement = ini.GetValueF("Hessian", "MIN_DISPLACEMENT", hessianMinDisplacement);
        
        dimerRotations = ini.GetValueL("Dimer", "ROTATIONS", dimerRotations);
        dimerSeparation = ini.GetValueF("Dimer", "SEPARATION", dimerSeparation);

        lanczosConvergenceLimit = ini.GetValueF("Lanczos", "CONVERGENCE", lanczosConvergenceLimit);
        lanczosMaxIterations = ini.GetValueF("Lanczos", "ITERATIONS", lanczosMaxIterations);
    }
    else
    {
        // GH: old parameter style is no more
        // should print some error or warning here
    }
}

// Accounting for the force calls
long Parameters::getForceCalls(){
    return forceCalls_;
}
long Parameters::getForceCallsSaddlePoint(){
    return forceCallsSaddlePointConcave_+forceCallsSaddlePointConvex_;
}
long Parameters::getForceCallsSaddlePointConcave(){
    return forceCallsSaddlePointConcave_;
}
long Parameters::getForceCallsSaddlePointConvex(){
    return forceCallsSaddlePointConvex_;
}
long Parameters::getForceCallsPrefactors(){
    return forceCallsPrefactors_;
}
void Parameters::addForceCalls(long forceCalls){
    forceCalls_ += forceCalls;
    return;
}
void Parameters::addForceCallsSaddlePoint(long forceCallsSaddlePoint, double eigenvalue){
    if(0 < eigenvalue)
        forceCallsSaddlePointConcave_ += forceCallsSaddlePoint;
    else
        forceCallsSaddlePointConvex_ += forceCallsSaddlePoint;
    return;
}
void Parameters::addForceCallsPrefactors(long forceCallsPrefactors){
    forceCallsPrefactors_ += forceCallsPrefactors;
    return;
}
void Parameters::resetForceCalls(){
    forceCalls_ = 0;
    return;
}
void Parameters::resetForceCallsSaddlePoint(){
    forceCallsSaddlePointConcave_ = 0;
    forceCallsSaddlePointConvex_ = 0;
    return;
}
void Parameters::resetForceCallsPrefactors(){
    forceCallsPrefactors_ = 0;
}


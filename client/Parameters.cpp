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

//using namespace constants;

Parameters::Parameters(){

    // Default values
    randomSeed_ = -1;
    reactantStateTag_ = 0;
    potentialTag_ = 1;
    potentialNoTranslation_ = 0;
    getPrefactorsTag_ = 0;
    typePerturbation_SP_ = 1;
    refine_SP_=false;
    lowestEigenmodeDetermination_SP_ = 1;
    minimize_only_ = 0;
    minimize_box_ = 0;
    
    // Tweakable parameters, default values if not read in from parameters_passed.dat

    // Value used in the Relaxation   
    converged_Relax_ = 0.005;

    // Values used in the Saddle Point determination   
    converged_SP_ = 0.025;
    maxJumpAttempts_SP_ = 0;
    nrOfTriesToDetermineSaddlePoint_SP_ = 1;
    maxStepSize_SP_ = 0.2;
//    maxStepSizeConcave_SP_ = 0.2;
//    maxStepSizeConvex_SP_ = 0.1;
    maxEnergy_SP_ = 20.0;
    normPerturbation_SP_ = 0.1;
    withinRadiusPerturbated_SP_ = 4.0;
    maxSinglePerturbation_SP_ = 0.1;
    maximumIterations_ = 512;

    // Values used in the Hessian determination   
	maxSize_Hessian_ = 0;
	minDisplacement_Hessian_ = 0.25;
    withinRadiusDisplaced_Hessian_ = 5.0;

    // Values used in the Dimer method
    rotations_Dimer_ = 1;
    //rotationsNewSearch_Dimer_ = 10; // GH: not using this any more
    separation_Dimer_ = 0.0001;
    rotationAngle_Dimer_ = 0.005;
    
    // Defaults for Lanczos
    iterationLimit_Lanczos_=50;
    convergenceLimit_Lanczos_=1e-4;

    // Initializing the cummulative output
    forceCalls_ = 0;
    forceCallsSaddlePointConcave_ = 0;
    forceCallsSaddlePointConvex_ = 0;
    forceCallsPrefactors_ = 0;
    displacement_saddle_distance_ = 0;

    return;
}


Parameters::~Parameters(){
    return;
}


void Parameters::load(FILE *file){
    
    //If we fail to parse the file as an INI, we will need to rewind the
    //file. So, we store the current stream position.
    fpos_t pos;
    fgetpos(file, &pos); 

    CIniFile ini;
    ini.CaseInsensitive();
    if(ini.ReadFile(file))
    {
        //If we succesfully read the file, then parse it as an INI
        randomSeed_ = ini.GetValueL("Default", "RANDOM_SEED", randomSeed_);
        reactantStateTag_ = ini.GetValueL("Default", "REACTANT_STATE_TAG", reactantStateTag_);
        potentialTag_ = ini.GetValueL("Default", "POTENTIAL_TAG", potentialTag_);
        potentialNoTranslation_ = ini.GetValueL("Default", "POTENTIAL_NO_TRANSLATION", potentialNoTranslation_);
        minimize_only_ = ini.GetValueL("Default", "MINIMIZE_ONLY", minimize_only_);
        minimize_box_ = ini.GetValueL("Default", "MINIMIZE_BOX", minimize_box_);
        getPrefactorsTag_ = ini.GetValueL("Default", "GET_PREFACTORS_TAG", getPrefactorsTag_);
        converged_Relax_ = ini.GetValueF("Default", "CONVERGED_RELAX", converged_Relax_);
        maximumIterations_ = ini.GetValueL("Default", "MAXIMUM_ITERATIONS", maximumIterations_);
        
        
        typePerturbation_SP_ = ini.GetValueL("Saddle_Point", "TYPE_PERTURBATION", typePerturbation_SP_);
        lowestEigenmodeDetermination_SP_ = ini.GetValueL("Saddle_Point", "LOWEST_EIGENMODE_DETERMINATION", lowestEigenmodeDetermination_SP_);
        refine_SP_ = ini.GetValueB("Saddle_Point", "REFINE", refine_SP_);
        converged_SP_ = ini.GetValueF("Saddle_Point", "CONVERGED", converged_SP_);
        maxJumpAttempts_SP_ = ini.GetValueL("Saddle_Point", "MAX_JUMP_ATTEMPTS", maxJumpAttempts_SP_);
        nrOfTriesToDetermineSaddlePoint_SP_ = ini.GetValueL("Saddle_Point", "NR_OF_TRIES_TO_DETERMINE_SADDLE_POINT", nrOfTriesToDetermineSaddlePoint_SP_);
        maxStepSize_SP_ = ini.GetValueF("Saddle_Point", "MAX_STEP_SIZE", maxStepSize_SP_);
        maxEnergy_SP_ = ini.GetValueF("Saddle_Point", "MAX_ENERGY", maxEnergy_SP_);
        normPerturbation_SP_ = ini.GetValueF("Saddle_Point", "NORM_PERTURBATION", normPerturbation_SP_);
        maxSinglePerturbation_SP_ = ini.GetValueF("Saddle_Point", "MAX_SINGLE_PERTURBATION", maxSinglePerturbation_SP_);
        withinRadiusPerturbated_SP_ = ini.GetValueF("Saddle_Point", "WITHIN_RADIUS_PERTURBATED", withinRadiusPerturbated_SP_);


        maxSize_Hessian_ = ini.GetValueL("Hessian", "MAX_SIZE", maxSize_Hessian_);
        withinRadiusDisplaced_Hessian_ = ini.GetValueF("Hessian", "WITHIN_RADIUS_DISPLACED", withinRadiusDisplaced_Hessian_);
        minDisplacement_Hessian_ = ini.GetValueF("Hessian", "MIN_DISPLACEMENT", minDisplacement_Hessian_);

        
        rotations_Dimer_ = ini.GetValueL("Dimer", "ROTATIONS", rotations_Dimer_);
        separation_Dimer_ = ini.GetValueF("Dimer", "SEPARATION", separation_Dimer_);

        convergenceLimit_Lanczos_ = ini.GetValueF("Lanczos", "CONVERGENCE", convergenceLimit_Lanczos_);
        iterationLimit_Lanczos_ = ini.GetValueF("Lanczos", "ITERATIONS", iterationLimit_Lanczos_);        
    }
    else
    {
        //Otherwise, parse it as an old-style parameters file
        printf("Reading old-style parameters file\n"); 
        //Rewind the file to before CIniFile::ReadFile was called
        fsetpos(file, &pos);

        char **parms;
        double *values;
        long nLines, i;
        
        nLines = linesInFile(file);
        
        values = new double[nLines];
        parms = new char*[nLines];
        for(int i=0; i<nLines; i++)
            parms[i] = new char[STRING_SIZE];
        i = 0;

        // Important that the variable value is float
        // The used function sscanf only support this type 
        float value;
        char parm[STRING_SIZE];
        char lineAll[STRING_SIZE];
        for (i=0; i<nLines; i++){
            fgets(lineAll, STRING_SIZE, file);
            std::sscanf(lineAll, "%s %e", parm, &value);            
            for(int j = 0; parm[j]; j++)
                parm[j] = toupper(parm[j]);
                
            values[i] = value;
            strcpy(parms[i], parm);
        }
        
        for(i=0; i<nLines; i++){
            // Note strcmp() return 0 if strings are equal
            if(!strcmp(parms[i], "RANDOM_SEED"))
                randomSeed_ = long(values[i]);
            else if(!strcmp(parms[i], "REACTANT_STATE_TAG"))
                reactantStateTag_ = long(values[i]);
            else if(!strcmp(parms[i], "POTENTIAL_TAG"))
                potentialTag_ = long(values[i]);
            else if(!strcmp(parms[i], "POTENTIAL_NO_TRANSLATION"))
                potentialNoTranslation_ = long(values[i]);
            else if(!strcmp(parms[i], "MINIMIZE_ONLY"))
                minimize_only_ = long(values[i]);
            else if(!strcmp(parms[i], "MINIMIZE_BOX"))
                minimize_box_ = long(values[i]);
            else if(!strcmp(parms[i], "GET_PREFACTORS_TAG"))
                getPrefactorsTag_ = long(values[i]);
            else if(!strcmp(parms[i], "TYPE_PERTURBATION_SP"))
                typePerturbation_SP_ = (long) values[i];
            else if(!strcmp(parms[i], "LOWEST_EIGENMODE_DETERMINATION_SP"))
                lowestEigenmodeDetermination_SP_ = long(values[i]);
            else if(!strcmp(parms[i], "REFINE_SP")) {
                refine_SP_ = (bool) values[i];
            }

            // Tweakable parameters

            // Relaxation related
            else if(!strcmp(parms[i], "CONVERGED_RELAX"))
                converged_Relax_ = values[i];

            // Saddle Point related
            else if(!strcmp(parms[i], "CONVERGED_SP"))
                converged_SP_ = values[i];
            else if(!strcmp(parms[i], "MAX_JUMP_ATTEMPTS_SP"))
                maxJumpAttempts_SP_ = (long) values[i];
            else if(!strcmp(parms[i], "NR_OF_TRIES_TO_DETERMINE_SADDLE_POINT_SP"))
                nrOfTriesToDetermineSaddlePoint_SP_ = long(values[i]);
            else if(!strcmp(parms[i], "MAX_STEP_SIZE_SP"))
                maxStepSize_SP_ = values[i];
/*            else if(!strcmp(parms[i], "MAX_STEP_SIZE_CONCAVE_SP"))
                maxStepSizeConcave_SP_ = values[i];
            else if(!strcmp(parms[i], "MAX_STEP_SIZE_CONVEX_SP"))
                maxStepSizeConvex_SP_ = values[i]; */
            else if(!strcmp(parms[i], "MAX_ENERGY_SP"))
                maxEnergy_SP_ = values[i];
            else if(!strcmp(parms[i], "NORM_PERTURBATION_SP"))
                normPerturbation_SP_ = values[i];
            else if(!strcmp(parms[i], "WITHIN_RADIUS_PERTURBATED_SP"))
                withinRadiusPerturbated_SP_ = values[i];
            else if(!strcmp(parms[i], "MAX_SINGLE_PERTURBATION_SP"))
                maxSinglePerturbation_SP_ = values[i];
            else if(!strcmp(parms[i], "MAXIMUM_ITERATIONS"))
                maximumIterations_ = (long)values[i];

            // Hessian related
            else if(!strcmp(parms[i], "MAX_SIZE_HESSIAN"))
                maxSize_Hessian_ = values[i];
            else if(!strcmp(parms[i], "MIN_DISPLACEMENT_HESSIAN"))
                minDisplacement_Hessian_ = values[i];
            else if(!strcmp(parms[i], "WITHIN_RADIUS_DISPLACED_HESSIAN"))
                withinRadiusDisplaced_Hessian_ = values[i];

            // Dimer related
            else if(!strcmp(parms[i], "ROTATIONS_DIMER"))
                rotations_Dimer_ = (long) values[i];
/*            else if(!strcmp(parms[i], "ROTATIONS_NEW_SEARCH_DIMER"))
                rotationsNewSearch_Dimer_ = (long) values[i];*/
            else if(!strcmp(parms[i], "SEPARATION_DIMER"))
                separation_Dimer_ = (double) values[i];

            else if(!strcmp(parms[i], "convergenceLimit_Lanczos_"))
                convergenceLimit_Lanczos_ = (double) values[i];
            else if(!strcmp(parms[i], "iterationLimit_Lanczos_"))
                iterationLimit_Lanczos_ = (double) values[i];
            // Lines with user comment are started with #
            else if(parms[i][0]=='#'){}
            else
                std::cout<<"Unknown property: "<<parms[i]<<"\n";
        }
        delete [] values;
        for(i=0; i<nLines; i++)
            delete [] parms[i];
        delete [] parms;
    }
    return;
}
void Parameters::saveOutput(FILE *file){
    fprintf(file, "%ld termination_reason\n", terminationReason_);
    fprintf(file, "%ld random_seed\n", randomSeed_);
    fprintf(file, "%ld reactant_state_tag\n", reactantStateTag_);
    fprintf(file, "%ld potential_tag\n", potentialTag_);
    fprintf(file, "%ld force_calls\n", forceCalls_);
    fprintf(file, "%ld force_calls_saddle_point_concave\n", forceCallsSaddlePointConcave_);
    fprintf(file, "%ld force_calls_saddle_point_convex\n", forceCallsSaddlePointConvex_);
    fprintf(file, "%f potential_energy_saddle\n", potentialEnergySP_);
    fprintf(file, "%f potential_energy_reactant\n", potentialEnergyMin1_);
    fprintf(file, "%f potential_energy_product\n", potentialEnergyMin2_);
    fprintf(file, "%f barrier_reactant_to_product\n", barrierReac_Prod_);
    fprintf(file, "%f barrier_product_to_reactant\n", barrierProd_Reac_);
    fprintf(file, "%f displacement_saddle_distance\n", displacement_saddle_distance_);

    if(getPrefactorsTag_)
    {
        fprintf(file, "%ld force_calls_prefactors\n", forceCallsPrefactors_);
        fprintf(file, "%f prefactor_reactant_to_product\n", prefactorReac_Prod_);
        fprintf(file, "%f prefactor_product_to_reactant\n", prefactorProd_Reac_);
    }
    
    return;
}


void Parameters::saveInput(FILE *file){
    fprintf(file, "REACTANT_STATE_TAG %ld\n\n",reactantStateTag_);
    fprintf(file, "# describtion RANDOM_SEED\n# -1 = time is used to seed the random generator\nRANDOM_SEED %ld\n\n", randomSeed_);
    fprintf(file, "# describtion POTENTIAL_TAG\n# 0 = a new potential provided by the user, code will crash at runtime if used and no new potential has been defined!\n# 1 = lennard-jones potential\n# 2 = morse\n# 3 = emt\n# 4 = edip\n# 5 = VASP\n# 6 = Tersoff\n# 7 = Stilling Weber\n# 8 = Lenosky\nPOTENTIAL_TAG %ld\n\n",potentialTag_);
    fprintf(file, "# describtion POTENTIAL_NO_TRANSLATION\n# 0 = the system is free to translate\n# 1 = translation of the system is prevented\nPOTENTIAL_NO_TRANSLATION %ld\n\n", potentialNoTranslation_);
    fprintf(file, "# describtion GET_PREFACTORS_TAG\n# 0 = processes prefactors are not determined\n# 1 = processes prefactors are determined\nGET_PREFACTORS_TAG %ld\n\n", getPrefactorsTag_);
    fprintf(file, "# describtion LOWEST_EIGENMODE_DETERMINATION_SP\n# 1 = dimer method\n# 2 = lanczos method\nLOWEST_EIGENMODE_DETERMINATION_SP %ld\n\n", lowestEigenmodeDetermination_SP_);
    fprintf(file, "# description TYPE_PERTURBATION_SP\n# 0 = no displacement\n# 1 = displace atoms not being fcc or hcp coordinated\n# 2 = displace atom with lowest coordination\n# 3 = displace last atom\nTYPE_PERTURBATION_SP %ld\n\n", typePerturbation_SP_);
    fprintf(file, "#  Refine saddle point\nREFINE_SP %d\n\n", refine_SP_);
    fprintf(file, "MINIMIZE_ONLY %ld\n\n", minimize_only_);
    fprintf(file, "MINIMIZE_BOX %ld\n\n", minimize_box_);
    
    fprintf(file, "#### Tweakable Parameters ####\n");
    fprintf(file, "# relaxation related\n");
    fprintf(file, "converged_Relax %f\n\n", converged_Relax_);
    
    fprintf(file, "# saddle point related\n");
    fprintf(file, "converged_SP %f\n", converged_SP_);
    fprintf(file, "maxJumpAttempts_SP %ld\n", maxJumpAttempts_SP_);
    fprintf(file, "nr_Of_Tries_To_Determine_Saddle_Point_SP %ld\n", nrOfTriesToDetermineSaddlePoint_SP_);
    fprintf(file, "max_Step_Size_SP %f\n", maxStepSize_SP_);
//    fprintf(file, "max_Step_Size_Concave_SP %f\n", maxStepSizeConcave_SP_);
//    fprintf(file, "max_Step_Size_Convex_SP %f\n", maxStepSizeConvex_SP_);
    fprintf(file, "max_Energy_SP %f\n", maxEnergy_SP_);
    fprintf(file, "norm_Perturbation_SP %f\n", normPerturbation_SP_);
    fprintf(file, "within_Radius_Perturbated_SP %f\n", withinRadiusPerturbated_SP_);
    fprintf(file, "max_Single_Perturbation_SP %f\n\n", maxSinglePerturbation_SP_);
    
    fprintf(file, "# hessian related\n");
    fprintf(file, "max_Size_Hessian %ld\n", maxSize_Hessian_);
    fprintf(file, "min_Displacement_Hessian %f\n", minDisplacement_Hessian_);
    fprintf(file, "within_Radius_Displaced_Hessian %f\n\n", withinRadiusDisplaced_Hessian_);

    fprintf(file, "# dimer related\n");
    fprintf(file, "rotations_Dimer %ld\n", rotations_Dimer_);
//    fprintf(file, "rotations_New_Search_Dimer %ld\n", rotationsNewSearch_Dimer_);
    fprintf(file, "rotationAngle_Dimer %f\n\n", rotationAngle_Dimer_);
    return;
}


void Parameters::printInput(){
    std::cout
    <<"---------Client side input---------\n"
    <<randomSeed_<<" : RANDOM_SEED\n"
    <<reactantStateTag_<<" : REACTANT_STATE_TAG\n"
    <<potentialTag_<<" : POTENTIAL_TAG\n"
    <<potentialNoTranslation_<<" : POTENTIAL_NO_TRANSLATION\n"
    <<getPrefactorsTag_<<" : GET_PREFACTORS_TAG\n"
    <<lowestEigenmodeDetermination_SP_<<" : LOWEST_EIGENMODE_DETERMINATION_SP\n" 
    <<typePerturbation_SP_<<" : TYPE_PERTURBATION_SP\n\n"
    << refine_SP_ <<" : REFINE_SP\n\n"
    << minimize_only_ <<" : MINIMIZE_ONLY\n\n"
    << minimize_box_ <<" : MINIMIZE_BOX\n\n"

    <<"#### Tweakable Parameters ####\n"
    <<"# relaxation related\n"
    <<converged_Relax_<<" : converged_Relax\n\n"

    <<"# saddle point related\n"
    <<maxJumpAttempts_SP_<<" : max_Jump_Attempts_SP\n"
    <<nrOfTriesToDetermineSaddlePoint_SP_<<" : nr_Of_Tries_To_Determine_Saddle_Point_SP\n"
    <<converged_SP_<<" : converged_SP\n"
    <<maxStepSize_SP_<<" : max_Step_Size_SP\n"
//    <<maxStepSizeConcave_SP_<<" : max_Step_Size_Concave_SP\n"
//    <<maxStepSizeConvex_SP_<<" : max_StepSize_Convex_SP\n"
    <<maxEnergy_SP_<<" : max_Energy_SP\n"
    <<normPerturbation_SP_<<" : norm_Perturbation_SP\n"
    <<withinRadiusPerturbated_SP_<<" : within_Radius_Perturbated_SP\n"
    <<maxSinglePerturbation_SP_<<" : max_Single_Perturbation_SP\n\n"
    
    <<"# hessian related\n"
    <<maxSize_Hessian_<<" : max_Size_Hessian\n"
    <<minDisplacement_Hessian_<<" : min_Displacement_Hessian\n"
    <<withinRadiusDisplaced_Hessian_<<" : within_Radius_Displaced_Hessian\n\n"
    
    <<"# dimer related\n"
    <<rotations_Dimer_<<" : rotations_Dimer\n"
//    <<rotationsNewSearch_Dimer_<<" : rotationsNewSearch_Dimer\n"
    <<rotationAngle_Dimer_<<" : rotationAngle_Dimer\n"
    <<"----------------------------------\n";
    return;        
}


void Parameters::printOutput(){
    std::cout
    <<"---------Client side data---------\n"
    <<randomSeed_<<" randomSeed\n"
    <<terminationReason_<<" terminationReason\n"
    <<reactantStateTag_<<" reactantStateTag\n"
    <<potentialTag_<<" potentialTag\n"
    <<potentialNoTranslation_<<" potentialNoTranslation\n"
    <<forceCalls_<<" forceCalls\n"
    <<forceCallsSaddlePointConcave_<<" forceCallsSaddlePointConcave\n"
    <<forceCallsSaddlePointConvex_<<" forceCallsSaddlePointConvex\n"
    <<potentialEnergySP_<<" potentialEnergySP\n"
    <<potentialEnergyMin1_<<" potentialEnergyMin1\n"
    <<potentialEnergyMin2_<<" potentialEnergyMin2\n"
    <<barrierReac_Prod_<<" barrierReac_Prod\n"
    <<barrierProd_Reac_<<" barrierProd_Reac\n"
    <<displacement_saddle_distance_<<" displacement_saddle_distance\n";
    
    if(getPrefactorsTag_){
        std::cout
        <<forceCallsPrefactors_<<" forceCallsPrefactors\n"
        <<prefactorReac_Prod_<<" prefactorReac_Prod\n"
        <<prefactorProd_Reac_<<" prefactorProd_Reac\n";
    }
    std::cout
    <<"----------------------------------\n";
    return;
}


// Passing the input parameters
long Parameters::getRandomSeed(){
    return randomSeed_;
}
void Parameters::setRandomSeed(long randomSeed){
    randomSeed_=randomSeed;
}
long Parameters::getReactantStateTag(){
    return reactantStateTag_;
}
long Parameters::getPotentialTag(){
    return potentialTag_;
}
long Parameters::getPotentialNoTranslation(){
    return potentialNoTranslation_;
}
bool Parameters::getMinimizeOnly(){
    return minimize_only_;
}
bool Parameters::getMinimizeBox(){
    return minimize_box_;
}
long Parameters::getPrefactorsTag(){
    return getPrefactorsTag_;
}
// Tweakable parameters

// Relaxation related
double Parameters::getConverged_Relax(){
    return converged_Relax_;
}
// Saddle Point determination related
double Parameters::getConverged_SP(){
    return converged_SP_;
}
long Parameters::getMaxJumpAttempts_SP(){
    return maxJumpAttempts_SP_;
}
long Parameters::getTypePerturbation_SP(){
    return typePerturbation_SP_;
}
bool Parameters::getRefineSP(){
    return refine_SP_;
}
long Parameters::getLowestEigenmodeDetermination_SP(){
    return lowestEigenmodeDetermination_SP_;
}
long Parameters::getNrOfTriesToDetermineSaddlePoints_SP(){
    return nrOfTriesToDetermineSaddlePoint_SP_;
}
double Parameters::getMaxStepSize_SP(){
    return maxStepSize_SP_;
}
/*double Parameters::getMaxStepSizeConcave_SP(){
    return maxStepSizeConcave_SP_;
}
double Parameters::getMaxStepSizeConvex_SP(){
    return maxStepSizeConvex_SP_;
}*/
double Parameters::getMaxEnergy_SP(){
    return maxEnergy_SP_;
}
double Parameters::getNormPerturbation_SP(){
    return normPerturbation_SP_;
}
double Parameters::getWithinRadiusPerturbated_SP(){
    return withinRadiusPerturbated_SP_;
}
double Parameters::getMaxSinglePerturbation_SP(){
    return maxSinglePerturbation_SP_;
}
/// Limit on the number of iterations that may be performed by the saddle point searches and minimization
long Parameters::getMaximumIterations() {
      return maximumIterations_;
}
// Hessian related
long Parameters::getMaxSize_Hessian(){
    return maxSize_Hessian_;
}
double Parameters::getMinDisplacement_Hessian(){
    return minDisplacement_Hessian_;
}
double Parameters::getWithinRadiusDisplaced_Hessian(){
    return withinRadiusDisplaced_Hessian_;
}
// Dimer related
long Parameters::getRotations_Dimer(){
    return rotations_Dimer_;
}
/*long Parameters::getRotationsNewSearch_Dimer(){
    return rotationsNewSearch_Dimer_;
}*/
double Parameters::getSeparation_Dimer(){
    return separation_Dimer_;
}
double Parameters::getRotationAngle_Dimer(){
    return rotationAngle_Dimer_;
}
// Lanczos related
double Parameters::getConvergenceLimit_Lanczos(){
      return convergenceLimit_Lanczos_;
}
double Parameters::getIterationLimit_Lanczos(){
      return iterationLimit_Lanczos_;
}

// Setting results values
void Parameters::setTerminationReason(long terminationReason){
    terminationReason_ = terminationReason;
    return;
}
void Parameters::setPotentialEnergySP(double potentialEnergySP){
    potentialEnergySP_ = potentialEnergySP;
    return;
}
void Parameters::setPotentialEnergyMin1(double potentialEnergyMin1){
    potentialEnergyMin1_ = potentialEnergyMin1;
    return;
}
void Parameters::setPotentialEnergyMin2(double potentialEnergyMin2){
    potentialEnergyMin2_ = potentialEnergyMin2;
    return;
}
void Parameters::setPrefactorReac_Prod(double prefactorReac_Prod){
    prefactorReac_Prod_ = prefactorReac_Prod;
    return;
}
void Parameters::setPrefactorProd_Reac(double prefactorProd_Reac){
    prefactorProd_Reac_ = prefactorProd_Reac;
    return;
}
void Parameters::setBarrierReac_Prod(double barrierReac_Prod){
    barrierReac_Prod_ = barrierReac_Prod;
    return;
}
void Parameters::setBarrierProd_Reac(double barrierProd_Reac){
    barrierProd_Reac_ = barrierProd_Reac;
    return;
}
void Parameters::setDisplacementSaddleDistance(double displacement_saddle_distance){
    displacement_saddle_distance_ = displacement_saddle_distance;
    return;
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

long Parameters::linesInFile(FILE *file){
    int nLines = 0;
    char ch, prev = '\n'; 
    // read all chars in the file
    while((ch = fgetc(file)) != EOF){
        if(ch == '\n'){
            nLines = nLines+1;
        }
        prev = ch; // keep a copy to later test whether
    }
    if(prev!='\n') // ...the last line did not end in a newline.
        nLines = nLines+1;

    rewind(file);
    return (nLines);
}


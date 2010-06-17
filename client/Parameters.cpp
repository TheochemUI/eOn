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
#include <cstring>
#include "Parameters.h"

//using namespace constants;

Parameters::Parameters(){    
    // Default values
    randomSeed_ = -1;
    reactantStateTag_ = 0;
    potentialTag_ = 1;
    potentialNoTranslation_ = 0;
    getPrefactorsTag_ = 0;
    typePertubation_SP_ = 1;
    refine_SP_=false;
    lowestEigenmodeDetermination_SP_ = 1;
    minimize_only_ = 0;
    minimize_box_ = 0;
    
    // Tweak able parameters, default values if not read in from parameters_passed.dat
    // Value used in the Relaxation   
    converged_Relax_ = 0.005;

    // Values used in the Saddle Point determination   
    converged_SP_ = 0.025;
    maxJumpAttempts_SP_ = 0;
    nrOfTriesToDetermineSaddlePoint_SP_ = 1;
    maxStepSizeConcave_SP_ = 0.2;
    maxStepSizeConvex_SP_ = 0.1;
    maxEnergy_SP_ = 20.0;
    normPertubation_SP_ = 0.1;
    withinRadiusPertubated_SP_ = 4.0;
    maxSinglePertubation_SP_ = 0.1;
    maximumIterations_ = 512;
    // Values used in the Hessian determination   
    minDisplacement_Hessian_ = 0.25;
    withinRadiusDisplaced_Hessian_ = 5.0;

    // Values used in the Dimer method
    rotations_Dimer_ = 1;
    rotationsNewSearch_Dimer_ = 10;
    dimer_dR_ = 0.0001;

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
    char **parms;
    double *values;
    long nLines, i;
    
    nLines = linesInFile(file);
    
    values = new double[nLines];
    parms = new char*[nLines];
    for(int i=0; i<nLines; i++)
        parms[i] = new char[STRING_SIZE];
    
    loadParameters(file, parms, values, nLines);
    
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
        else if(!strcmp(parms[i], "TYPE_PERTUBATION_SP"))
            typePertubation_SP_ = (long) values[i];
        else if(!strcmp(parms[i], "LOWEST_EIGENMODE_DETERMINATION_SP"))
            lowestEigenmodeDetermination_SP_ = long(values[i]);
        else if(!strcmp(parms[i], "REFINE_SP")) {
            refine_SP_ = (bool) values[i];
        }
        // Tweak able parameters
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
        else if(!strcmp(parms[i], "MAX_STEP_SIZE_CONCAVE_SP"))
            maxStepSizeConcave_SP_ = values[i];
        else if(!strcmp(parms[i], "MAX_STEP_SIZE_CONVEX_SP"))
            maxStepSizeConvex_SP_ = values[i];            
        else if(!strcmp(parms[i], "MAX_ENERGY_SP"))
            maxEnergy_SP_ = values[i];
        else if(!strcmp(parms[i], "NORM_PERTUBATION_SP"))
            normPertubation_SP_ = values[i];
        else if(!strcmp(parms[i], "WITHIN_RADIUS_PERTUBATED_SP"))
            withinRadiusPertubated_SP_ = values[i];
        else if(!strcmp(parms[i], "MAX_SINGLE_PERTUBATION_SP"))
            maxSinglePertubation_SP_ = values[i];
        else if(!strcmp(parms[i], "MAXIMUM_ITERATIONS"))
              maximumIterations_ = (long)values[i];
             // Hessian related
        else if(!strcmp(parms[i], "MIN_DISPLACEMENT_HESSIAN"))
            minDisplacement_Hessian_ = values[i];
        else if(!strcmp(parms[i], "WITHIN_RADIUS_DISPLACED_HESSIAN"))
            withinRadiusDisplaced_Hessian_ = values[i];
        // Dimer related
        else if(!strcmp(parms[i], "ROTATIONS_DIMER"))
            rotations_Dimer_ = (long) values[i];
        else if(!strcmp(parms[i], "ROTATIONS_NEW_SEARCH_DIMER"))
            rotationsNewSearch_Dimer_ = (long) values[i];
        else if(!strcmp(parms[i], "DIMER_DR"))
            dimer_dR_ = (double) values[i];

        // Lines with user comment are started with #
        else if(parms[i][0]=='#'){}
        else
            std::cout<<"Unknown property: "<<parms[i]<<"\n";
    }
    delete [] values;
    for(i=0; i<nLines; i++)
        delete [] parms[i];
    delete [] parms;
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
    fprintf(file, "# description TYPE_PERTUBATION_SP\n# 0 = no displacement\n# 1 = displace atoms not being fcc or hcp coordinated\n# 2 = displace atom with lowest coordination\n# 3 = displace last atom\nTYPE_PERTUBATION_SP %ld\n\n", typePertubation_SP_);
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
    fprintf(file, "max_Step_Size_Concave_SP %f\n", maxStepSizeConcave_SP_);
    fprintf(file, "max_Step_Size_Convex_SP %f\n", maxStepSizeConvex_SP_);
    fprintf(file, "max_Energy_SP %f\n", maxEnergy_SP_);
    fprintf(file, "norm_Pertubation_SP %f\n", normPertubation_SP_);
    fprintf(file, "within_Radius_Pertubated_SP %f\n", withinRadiusPertubated_SP_);
    fprintf(file, "max_Single_Pertubation_SP %f\n\n", maxSinglePertubation_SP_);
    
    fprintf(file, "# hessian related\n");
    fprintf(file, "min_Displacement_Hessian %f\n", minDisplacement_Hessian_);
    fprintf(file, "within_Radius_Displaced_Hessian %f\n\n", withinRadiusDisplaced_Hessian_);

    fprintf(file, "# dimer related\n");
    fprintf(file, "rotations_Dimer %ld\n", rotations_Dimer_);
    fprintf(file, "rotations_New_Search_Dimer %ld\n", rotationsNewSearch_Dimer_);
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
    <<typePertubation_SP_<<" : TYPE_PERTUBATION_SP\n\n"
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
    <<maxStepSizeConcave_SP_<<" : max_Step_Size_Concave_SP\n"
    <<maxStepSizeConvex_SP_<<" : max_StepSize_Convex_SP\n"
    <<maxEnergy_SP_<<" : max_Energy_SP\n"
    <<normPertubation_SP_<<" : norm_Pertubation_SP\n"
    <<withinRadiusPertubated_SP_<<" : within_Radius_Pertubated_SP\n"
    <<maxSinglePertubation_SP_<<" : max_Single_Pertubation_SP\n\n"
    
    <<"# hessian related\n"
    <<minDisplacement_Hessian_<<" : min_Displacement_Hessian\n"
    <<withinRadiusDisplaced_Hessian_<<" : within_Radius_Displaced_Hessian\n\n"
    
    <<"# dimer related\n"
    <<rotations_Dimer_<<" :rotations_Dimer\n"
    <<rotationsNewSearch_Dimer_<<" :rotationsNewSearch_Dimer\n"
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


// Passing the input parameters.
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
// Tweak able paramters
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
long Parameters::getTypePertubation_SP(){
    return typePertubation_SP_;
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
double Parameters::getMaxStepSizeConcave_SP(){
    return maxStepSizeConcave_SP_;
}
double Parameters::getMaxStepSizeConvex_SP(){
    return maxStepSizeConvex_SP_;
}
double Parameters::getMaxEnergy_SP(){
    return maxEnergy_SP_;
}
double Parameters::getNormPertubation_SP(){
    return normPertubation_SP_;
}
double Parameters::getWithinRadiusPertubated_SP(){
    return withinRadiusPertubated_SP_;
}
double Parameters::getMaxSinglePertubation_SP(){
    return maxSinglePertubation_SP_;
}
/// Limit on the number of iterations that may be performed by the saddle point searches and minimisation
long Parameters::getMaximumIterations() {
      return maximumIterations_;
}
// Hessian related
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
long Parameters::getRotationsNewSearch_Dimer(){
    return rotationsNewSearch_Dimer_;
}

double Parameters::getDimerDR()
{
    return dimer_dR_;
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

void Parameters::loadParameters(FILE *file, char **parms, double *values, long nLines){
    int i = 0;

    // Important that the variable value is float
    // The used function sscanf only support this type 
    float value;
    char parm[STRING_SIZE];
    char lineAll[STRING_SIZE];
    while (nLines){
        fgets(lineAll, STRING_SIZE, file);
        std::sscanf(lineAll, "%s %e", parm, &value);            
        for(int j = 0; parm[j]; j++)
            parm[j] = toupper(parm[j]);
            
        values[i] = value;
        strcpy(parms[i], parm);
        i = i+1;
        nLines = nLines-1;
    }
    return;
}


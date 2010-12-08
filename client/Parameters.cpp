//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <errno.h>
#include <time.h>
#include "Parameters.h"
#include "INIFile.h"
#include "HelperFunctions.h"
#include "Hessian.h"
#include "Job.h"
#include "Dynamics.h"
#include "Potentials.h"
#include "SaddlePoint.h"

Parameters::Parameters(){

    // Default values

    jobType = Job::PROCESS_SEARCH;
    randomSeed = -1;
    potentialType = Potential::POT_LJ;
    getPrefactor = false;
    minimizeOnly = false;
    minimizeBox = false;
    maxDifferencePos = 0.1;

    // for finding Epicenters
    neighborCutoff = 3.3;

    //Debug Options
    saveStdout = false;

    // default parameters for process search
    processSearchMinimizeFirst = false;

    // default parameters for saddle point determination
    saddleDisplacementType = SaddlePoint::DISP_MIN_COORDINATED;
    saddleRefine = false;
    saddleMinModeMethod = SaddlePoint::MINMODE_DIMER;
    saddleConverged = 0.025;
    saddleMaxJumpAttempts = 0;
    saddleMaxStepSize = 0.2;
    saddleMaxEnergy = 20.0;
    saddleNormPerturbation = 0.1;
    saddleWithinRadiusPerturbated = 4.0;
    saddleMaxSinglePerturbation = 0.1;
    saddleMaxIterations = 1000;
    saddlePerpForceRatio = 0.0;

    // default parameters for Hessian determination
    hessianType = Hessian::REACTANT;
    hessianMinDisplacement = 0.25;
    hessianWithinRadius = 5.0;
    hessianPrefactorMax = 10e20;
    hessianPrefactorMin = 10e8;

    // default parameters the dimer method
    dimerSeparation = 0.0001;
    dimerRotationAngle = 0.005;
    dimerWindowHigh = 1.0;
    dimerWindowLow = 0.1;
    dimerRotationsHigh = 8;
    dimerRotationsLow = 1;

    // default parameters used by the optimizers
    maximumIterations=512;
    cgCurvatureStep = 0.001;
    cgMaxMoveFullRelax = 0.2;
    qmTimeStep = 0.1;
    convergedRelax = 0.005;

    // default parameters used by parallel repica dynamics
    mdTimeStep = 0.1;
    mdTemperature = 300.0;
    mdSteps = 1000;
    mdMaxMovedDist = 2.0;
    mdRefine = false;
    mdAutoStop = false;
    mdRecordAccuracy = 1;
    mdRefineAccuracy = 10;
    mdCheckFreq = 1000;
    mdRelaxSteps = 500;
    mdDephaseSteps = 200;
    mdDephaseLoopStop = false;
    mdDephaseLoopMax = 5;

    // default parameters used by hyperdynamics
    bondBoost = false ;
    bondBoostDVMAX = 0.0;
    bondBoostQRR = 0.0001; // can not be set to 0
    bondBoostPRR = 0.95;
    bondBoostQcut = 3.0;
    bondBoostRMDS = 0;

    // default parameters used by thermostat
    thermoType = Dynamics::ANDERSEN;
    thermoAndersenAlpha = 0.2; // collision strength
    thermoAndersenTcol = 10; // collision frequency in unit of dt
    thermoNoseMass = 0;
    return;

    // Default parameters for the displacement sampling job
    displaceNSamples = 32; // The number of samples to take.
    displaceIterMax = 32; // The maximum number of rotations to perform on the dimer.
    displaceTorqueConvergence = 0.01; // The convergence criteria of the dimer rotation.
    displaceMaxCurvature = -0.1; // The maximum curvature for which a sample is considered good. Used to avoid shallow but negative curvatures.
    displaceMaxDE = 10.0; // The maximum dE for which a sample is considered good. XXX: Should use saddleMaxEnergy?
    displaceCutoffs = "0.0 3.3";
    displaceMagnitudes = "0.0625 0.125 0.25";
}

Parameters::~Parameters(){
    return;
}

string Parameters::toLowerCase(string s)
{
    for (string::size_type i = 0; i < s.length(); ++i) {
      s[i] = tolower(s[i]);
    }
    return s;
}

int Parameters::load(string filename)
{
    FILE *fh;

    fh = fopen(filename.c_str(), "rb");
    if (fh == NULL) {
        return 1;
    }
    int error = load(fh);
    fclose(fh);
    return error;
}

int Parameters::load(FILE *file){
    
    CIniFile ini;
    ini.CaseInsensitive();
    int error=0;
    printf("List of non-default parameters: ");
    if(ini.ReadFile(file))
    {
        printf("\n\n");
        // if we succesfully read the file, then parse it as an INI

        randomSeed = ini.GetValueL("Default", "random_seed", randomSeed);
        // Initialize random generator
        if(randomSeed < 0){
            unsigned i = time(NULL);
            randomSeed = i;
            helper_functions::random(i);
        }else{
            helper_functions::random(randomSeed);
        }
        printf("Random seed is: %ld\n", randomSeed);

        string potentialTypeString = ini.GetValue("Default", "potential_type");
        potentialTypeString = toLowerCase(potentialTypeString);
        if(potentialTypeString == "user"){
            potentialType = Potential::POT_USER;
        }else if(potentialTypeString == "lj"){
            potentialType = Potential::POT_LJ;
        }else if(potentialTypeString == "morse"){
            potentialType = Potential::POT_MORSE;
        }else if(potentialTypeString == "emt"){
            potentialType = Potential::POT_EMT;
        }else if(potentialTypeString == "edip"){
            potentialType = Potential::POT_EDIP;
        }else if(potentialTypeString == "vasp"){
            potentialType = Potential::POT_VASP;
        }else if(potentialTypeString == "tersoff"){
            potentialType = Potential::POT_TERSOFF;
        }else if(potentialTypeString == "sw"){
            potentialType = Potential::POT_SW;
        }else if(potentialTypeString == "lenosky"){
            potentialType = Potential::POT_LENOSKY;
        }else if(potentialTypeString == "ljbinary"){
            potentialType = Potential::POT_LJBINARY;
        }else if(potentialTypeString == "aluminum"){
            potentialType = Potential::POT_ALUMINUM;
        }else if(potentialTypeString == "eam"){
            potentialType = Potential::POT_EAM;
        }else if(potentialTypeString == "qsc"){
            potentialType = Potential::POT_QSC;
        }else if(potentialTypeString == "zpice"){
            potentialType = Potential::POT_ZPICE;
        }else if(potentialTypeString == "tip4p"){
            potentialType = Potential::POT_TIP4P;
        }else if(potentialTypeString == "bopfox"){
            potentialType = Potential::POT_BOPFOX;
        }else if(potentialTypeString == "bop"){
            potentialType = Potential::POT_BOP;
        }else{
            fprintf(stderr, "Unknown potential_type: %s\n", potentialTypeString.c_str());
            error = 1;
        }

        minimizeOnly = ini.GetValueB("Default", "minimize_only", minimizeOnly);
        minimizeBox = ini.GetValueB("Default", "minimze_box", minimizeBox);
        convergedRelax = ini.GetValueF("Default", "converged_relax", convergedRelax);
        maximumIterations = ini.GetValueL("Default", "maximum_iterations", maximumIterations);
        cgCurvatureStep = ini.GetValueF("CG","curvature_step", cgCurvatureStep);
        cgMaxMoveFullRelax = ini.GetValueF("CG","max_move", cgMaxMoveFullRelax);
        qmTimeStep = ini.GetValueF("QM","time_step", qmTimeStep);
        
        maxDifferencePos = ini.GetValueF("Default", "max_difference_pos", maxDifferencePos);
        neighborCutoff = ini.GetValueF("Default", "neighbor_cutoff", neighborCutoff);

        string jobTypeString;
        jobTypeString = ini.GetValue("Default", "job_type", "process_search");
        jobTypeString = toLowerCase(jobTypeString);
        if (jobTypeString == "process_search") {
            jobType = Job::PROCESS_SEARCH;
        }else if (jobTypeString == "saddle_search") {
            jobType = Job::SADDLE_SEARCH;
        }else if (jobTypeString == "minimization") {
            jobType = Job::MINIMIZATION;
        }else if (jobTypeString == "hessian") {
            jobType = Job::HESSIAN;
        }else if (jobTypeString == "parallel_replica"){
            jobType = Job::PARALLEL_REPLICA;
        }else if (jobTypeString == "replica_exchange"){
            jobType = Job::REPLICA_EXCHANGE;         
        }else if (jobTypeString == "dimer_dr"){
            jobType = Job::DIMER_DR;
        }else if (jobTypeString == "dimer_rotation"){
            jobType = Job::DIMER_ROTATION;
        }else if (jobTypeString == "displacement_sampling"){
            jobType = Job::DISPLACEMENT_SAMPLING;
        }else if (jobTypeString == "test"){
            jobType = Job::TEST;
        }
        else{
            fprintf(stderr, "Unknown job_type: %s\n", jobTypeString.c_str());
            error = 1;
        }

        processSearchMinimizeFirst = ini.GetValueB("Process_Search", "minimize_first", processSearchMinimizeFirst);
        saveStdout= ini.GetValueB("Debug", "save_stdout", saveStdout);

        string minModeMethodString = ini.GetValue("Saddle_Point", "minmode_method");
        minModeMethodString = toLowerCase(minModeMethodString);
        if(minModeMethodString == "dimer"){
            saddleMinModeMethod = SaddlePoint::MINMODE_DIMER;
        }else if(minModeMethodString == "lanczos"){
            saddleMinModeMethod = SaddlePoint::MINMODE_LANCZOS;
        }

        string displacementTypeString = ini.GetValue("Saddle_Point", "displacement_type");
        displacementTypeString = toLowerCase(displacementTypeString);
        if(displacementTypeString == "none"){ 
            saddleDisplacementType = SaddlePoint::DISP_NONE;
        }else if(displacementTypeString == "not_fcc_or_hcp"){
            saddleDisplacementType = SaddlePoint::DISP_NOT_FCC_OR_HCP;
        }else if(displacementTypeString == "min_coordinated"){
            saddleDisplacementType = SaddlePoint::DISP_MIN_COORDINATED;
        }else if(displacementTypeString == "not_last_atom"){
            saddleDisplacementType = SaddlePoint::DISP_LAST_ATOM;
        }

        saddleRefine = ini.GetValueB("Saddle_Point", "refine", saddleRefine); 
        saddleConverged = ini.GetValueF("Saddle_Point", "converged", saddleConverged);
        saddleMaxJumpAttempts = ini.GetValueL("Saddle_Point", "max_jump_attempts", saddleMaxJumpAttempts);
        saddleMaxStepSize = ini.GetValueF("Saddle_Point", "max_step_size", saddleMaxStepSize);
        saddleMaxEnergy = ini.GetValueF("Saddle_Point", "max_energy", saddleMaxEnergy);
        saddleNormPerturbation = ini.GetValueF("Saddle_Point", "norm_perturbation", saddleNormPerturbation);
        saddleMaxSinglePerturbation = ini.GetValueF("Saddle_Point", "max_single_perturbation", saddleMaxSinglePerturbation);
        saddleWithinRadiusPerturbated = ini.GetValueF("Saddle_Point", "within_radius_perturbated", saddleWithinRadiusPerturbated);
        saddleMaxIterations = ini.GetValueL("Saddle_Point", "max_iterations", saddleMaxIterations);
        saddlePerpForceRatio = ini.GetValueF("Saddle_Point", "perp_force_ratio", saddlePerpForceRatio);

        string hessianTypeString = ini.GetValue("Hessian", "type", "reactant");
        hessianTypeString = toLowerCase(hessianTypeString);
        if(hessianTypeString == "reactant"){
            hessianType = Hessian::REACTANT;
        }else if(hessianTypeString == "saddle"){
            hessianTypeString = Hessian::SADDLE;
        }else if(hessianTypeString == "product"){
            hessianTypeString = Hessian::PRODUCT;
        }

        hessianWithinRadius = ini.GetValueF("Hessian", "within_radius", hessianWithinRadius);
        hessianMinDisplacement = ini.GetValueF("Hessian", "min_displacement", hessianMinDisplacement);
 
        dimerSeparation = ini.GetValueF("Dimer", "separation", dimerSeparation);
        dimerRotationAngle = ini.GetValueF("Dimer", "finite_diff_angle", dimerRotationAngle);
        dimerRotationsHigh = ini.GetValueL("Dimer", "rotations_high", dimerRotationsHigh);
        dimerRotationsLow = ini.GetValueL("Dimer", "rotations_low", dimerRotationsLow);
        dimerWindowHigh = ini.GetValueF("Dimer", "window_high", dimerWindowHigh);
        dimerWindowLow = ini.GetValueF("Dimer", "window_low", dimerWindowLow);

        lanczosConvergence = ini.GetValueF("Lanczos", "convergence", lanczosConvergence);
        lanczosIteration = ini.GetValueL("Lanczos", "iteration", lanczosIteration);

        displaceNSamples = ini.GetValueL("DisplacementSampling", "samples", displaceNSamples);
        displaceIterMax = ini.GetValueL("DisplacementSampling", "iter_max", displaceIterMax);
        displaceTorqueConvergence = ini.GetValueF("DisplacementSampling", "torque_convergence", displaceTorqueConvergence);
        displaceMaxCurvature = ini.GetValueF("DisplacementSampling", "max_curvature", displaceMaxCurvature);
        displaceMaxDE = ini.GetValueF("DisplacementSampling", "max_de", displaceMaxDE);
        displaceCutoffs = ini.GetValue("DisplacementSampling", "cutoffs", displaceCutoffs);
        displaceMagnitudes = ini.GetValue("DisplacementSampling", "magnitudes", displaceMagnitudes);

        mdTimeStep = ini.GetValueF("Dynamics","time_step",mdTimeStep);
        mdTemperature = ini.GetValueF("Dynamics","temperature",mdTemperature);
        mdSteps = ini.GetValueL("Dynamics","steps",mdSteps);
        mdDephaseSteps = ini.GetValueL("Dynamics","dephase_steps",mdDephaseSteps);
        mdMaxMovedDist = ini.GetValueF("Dynamics","parrep_max_moved_dist",mdMaxMovedDist);  
        mdRefine = ini.GetValueB("Dynamics","refine",mdRefine);
        mdAutoStop = ini.GetValueB("Dynamics","auto_stop",mdAutoStop);
        mdRecordAccuracy = ini.GetValueL("Dynamics","record_accuracy",mdRecordAccuracy);
        mdRefineAccuracy = ini.GetValueL("Dynamics","refine_accuracy",mdRefineAccuracy);
        mdCheckFreq = ini.GetValueL("Dynamics","check_freq", mdCheckFreq);
        mdRelaxSteps = ini.GetValueL("Dynamics","new_relax_step",mdRelaxSteps);
        mdDephaseLoopStop = ini.GetValueB("Dynamics","dephase_loop_stop",mdDephaseLoopStop);
        mdDephaseLoopMax = ini.GetValueL("Dynamics","dephase_loop_max",mdDephaseLoopMax);
         
        string thermoTypeString;
        thermoTypeString = ini.GetValue("Dynamics", "thermo_type", "andersen");
        thermoTypeString = toLowerCase(thermoTypeString);
        if (thermoTypeString == "andersen") {
            thermoType = Dynamics::ANDERSEN;
        }else if (thermoTypeString == "nosehover") {
            thermoType = Dynamics::NOSE_HOVER;
        }

        thermoAndersenAlpha = ini.GetValueF("Dynamics","andersen_alpha",thermoAndersenAlpha);
        thermoAndersenTcol = ini.GetValueF("Dynamics","andersen_tcol",thermoAndersenTcol);
        thermoNoseMass = ini.GetValueF("Dynamics","nose_mass",thermoNoseMass);

        bondBoost = ini.GetValueB("Hyper","bond_boost",bondBoost);
        bondBoostRMDS = ini.GetValueL("Hyper","RMDS",bondBoostRMDS);
        bondBoostDVMAX = ini.GetValueF("Hyper","DVMAX",bondBoostDVMAX);
        bondBoostQRR = ini.GetValueF("Hyper","QRR",bondBoostQRR );
        bondBoostPRR = ini.GetValueF("Hyper","PRR",bondBoostPRR );
        bondBoostQcut= ini.GetValueF("Hyper","Qcut",bondBoostQcut);

    }
    else
    {
        fprintf(stderr, "Couldn't parse the ini file. Perhaps you are "
                        "using the old style config?\n");
        error = 1;
    }
    return error;
}

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

    // [Main]
    job = Job::PROCESS_SEARCH;
    randomSeed = -1;
    potential = Potential::POT_LJ_RH;
    calculatePrefactor = true;
    minimizeOnly = false;
    minimizeBox = false;
    maxDifferencePos = 0.1;

    // for finding Epicenters
    neighborCutoff = 3.3;

    // [Debug] options
    saveStdout = false;

    // [Process Search] default parameters
    processSearchMinimizeFirst = false;

    // [Saddle Search] default parameters
    saddleDisplaceType = SaddlePoint::DISP_MIN_COORDINATED;
    saddleMinmodeMethod = SaddlePoint::MINMODE_DIMER;
    saddleConvergedForce = 0.025;
    saddleMaxStepSize = 0.2;
    saddleMaxEnergy = 20.0;
    saddleMaxIterations = 1000;
    // undocumented
    saddleDisplace = false;
    saddleNormDisplace = 0.1;
    saddleWithinRadiusDisplace = 4.0;
    saddleMaxSingleDisplace = 0.1;
    saddleMaxJumpAttempts = 0;
    saddlePerpForceRatio = 0.0;

    // [Hessian] determination
    hessianType = Hessian::REACTANT;
    hessianMinDisplacement = 0.25;
    hessianWithinRadius = 5.0;
    hessianPrefactorMax = 10e20;
    hessianPrefactorMin = 10e8;

    // [Dimer] method
    dimerSeparation = 0.0001;
    dimerRotationAngle = 0.005;
    dimerWindowMax = 1.0;
    dimerWindowMin = 0.1;
    dimerRotationsMax = 8;
    dimerRotationsMin = 1;

    // [Minimizers]
    maximumIterations=512;
    convergedForceRelax = 0.005;

    // [Conjugate Gradients]
    cgCurvatureStep = 0.001;
    cgMaxMoveFullRelax = 0.2;

    //[Quick Min]
    qmTimeStep = 0.1;
    qmMaxMove = 0.2;

    // [Parallel Replica]
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

    // [Hyperdynamics]
    bondBoost = false ;
    bondBoostDVMAX = 0.0;
    bondBoostQRR = 0.0001; // can not be set to 0
    bondBoostPRR = 0.95;
    bondBoostQcut = 3.0;
    bondBoostRMDS = 0;

    // [Thermostat]
    thermostat = Dynamics::ANDERSEN;
    thermoAndersenAlpha = 0.2; // collision strength
    thermoAndersenTcol = 10; // collision frequency in unit of dt
    thermoNoseMass = 0;

    // [Displacement Sampling]
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
    //printf("List of non-default parameters: ");
    if(ini.ReadFile(file))
    {
        printf("\n\n");
        // if we succesfully read the file, then parse it as an INI

        randomSeed = ini.GetValueL("Main", "random_seed", randomSeed);
        // Initialize random generator
        if(randomSeed < 0){
            unsigned i = time(NULL);
            randomSeed = i;
            helper_functions::random(i);
        }else{
            helper_functions::random(randomSeed);
        }
        printf("Random seed is: %ld\n", randomSeed);

        string potentialString = ini.GetValue("Main", "potential");
        potentialString = toLowerCase(potentialString);
        if(potentialString == "lj_rh"){
            potential = Potential::POT_LJ_RH;
        }else if(potentialString == "morse_pt"){
            potential = Potential::POT_MORSE_PT;
        }else if(potentialString == "emt"){
            potential = Potential::POT_EMT;
        }else if(potentialString == "edip"){
            potential = Potential::POT_EDIP;
        }else if(potentialString == "vasp"){
            potential = Potential::POT_VASP;
        }else if(potentialString == "tersoff_si"){
            potential = Potential::POT_TERSOFF_SI;
        }else if(potentialString == "sw_si"){
            potential = Potential::POT_SW_SI;
        }else if(potentialString == "lenosky_si"){
            potential = Potential::POT_LENOSKY_SI;
        }else if(potentialString == "ljbinary"){
            potential = Potential::POT_LJBINARY;
        }else if(potentialString == "eam_al"){
            potential = Potential::POT_EAM_AL;
        }else if(potentialString == "qsc"){
            potential = Potential::POT_QSC;
        }else if(potentialString == "zpice"){
            potential = Potential::POT_ZPICE;
        }else if(potentialString == "tip4p"){
            potential = Potential::POT_TIP4P;
        }else if(potentialString == "bopfox"){
            potential = Potential::POT_BOPFOX;
        }else if(potentialString == "bop"){
            potential = Potential::POT_BOP;
        }else{
            fprintf(stderr, "Unknown potential: %s\n", potentialString.c_str());
            error = 1;
        }

        minimizeOnly = ini.GetValueB("Main", "minimize_only", minimizeOnly);
        minimizeBox = ini.GetValueB("Main", "minimze_box", minimizeBox);
        
        maxDifferencePos = ini.GetValueF("Default", "max_difference_pos", maxDifferencePos);
        neighborCutoff = ini.GetValueF("Default", "neighbor_cutoff", neighborCutoff);

        string jobString;
        jobString = ini.GetValue("Main", "job", "process_search");
        jobString = toLowerCase(jobString);
        if (jobString == "process_search") {
            job = Job::PROCESS_SEARCH;
        }else if (jobString == "saddle_search") {
            job = Job::SADDLE_SEARCH;
        }else if (jobString == "minimization") {
            job = Job::MINIMIZATION;
        }else if (jobString == "hessian") {
            job = Job::HESSIAN;
        }else if (jobString == "parallel_replica"){
            job = Job::PARALLEL_REPLICA;
        }else if (jobString == "replica_exchange"){
            job = Job::REPLICA_EXCHANGE;         
        }else if (jobString == "dimer_dr"){
            job = Job::DIMER_DR;
        }else if (jobString == "dimer_rotation"){
            job = Job::DIMER_ROTATION;
        }else if (jobString == "displacement_sampling"){
            job = Job::DISPLACEMENT_SAMPLING;
        }else if (jobString == "test"){
            job = Job::TEST;
        }
        else{
            fprintf(stderr, "Unknown job_type: %s\n", jobString.c_str());
            error = 1;
        }

        saveStdout= ini.GetValueB("Debug", "save_stdout", saveStdout);

        processSearchMinimizeFirst = ini.GetValueB("Process Search", "minimize_first", processSearchMinimizeFirst);

        string minmodeMethodString = ini.GetValue("Saddle Search", "minmode_method", "dimer");
        minmodeMethodString = toLowerCase(minmodeMethodString);
        if(minmodeMethodString == "dimer"){
            saddleMinmodeMethod = SaddlePoint::MINMODE_DIMER;
        }else if(minmodeMethodString == "lanczos"){
            saddleMinmodeMethod = SaddlePoint::MINMODE_LANCZOS;
        }
        string displaceString = ini.GetValue("Saddle Search", "displace", "none"); // undocumented
        displaceString = toLowerCase(displaceString);
        if(displaceString == "none"){ 
            saddleDisplaceType = SaddlePoint::DISP_NONE;
        }else if(displaceString == "not_fcc_or_hcp"){
            saddleDisplaceType = SaddlePoint::DISP_NOT_FCC_OR_HCP;
        }else if(displaceString == "min_coordinated"){
            saddleDisplaceType = SaddlePoint::DISP_MIN_COORDINATED;
        }else if(displaceString == "not_last_atom"){
            saddleDisplaceType = SaddlePoint::DISP_LAST_ATOM;
        }
        saddleConvergedForce = ini.GetValueF("Saddle Search", "converged_force", saddleConvergedForce);
        saddleMaxStepSize = ini.GetValueF("Saddle Search", "max_step_size", saddleMaxStepSize);
        saddleMaxEnergy = ini.GetValueF("Saddle Search", "max_energy", saddleMaxEnergy);
        saddleMaxIterations = ini.GetValueL("Saddle Search", "max_iterations", saddleMaxIterations);

        saddleDisplace = ini.GetValueB("Saddle Search", "displace", saddleDisplace); //undocumented
        saddleMaxJumpAttempts = ini.GetValueL("Saddle Search", "max_jump_attempts", saddleMaxJumpAttempts); //undocumented
        saddleMaxSingleDisplace = ini.GetValueF("Saddle Search", "max_single_displace", saddleMaxSingleDisplace); //undocumented
        saddleWithinRadiusDisplace= ini.GetValueF("Saddle Search", "within_radius_displace", saddleWithinRadiusDisplace); //undocumented
        saddleNormDisplace = ini.GetValueF("Saddle Search", "norm_displace", saddleNormDisplace); //undocumented
        saddlePerpForceRatio = ini.GetValueF("Saddle Search", "perp_force_ratio", saddlePerpForceRatio); //undocumented

        string hessianString = ini.GetValue("Hessian", "type", "reactant");
        hessianString = toLowerCase(hessianString);
        if(hessianString == "reactant"){
            hessianType = Hessian::REACTANT;
        }else if(hessianString == "saddle"){
            hessianType = Hessian::SADDLE;
        }else if(hessianString == "product"){
            hessianType = Hessian::PRODUCT;
        }

        convergedForceRelax = ini.GetValueF("Minimizer", "converged_force", convergedForceRelax);
        maximumIterations = ini.GetValueL("Minimizer", "max_iterations", maximumIterations);
        cgCurvatureStep = ini.GetValueF("Conjugate Gradients","curvature_step", cgCurvatureStep);
        cgMaxMoveFullRelax = ini.GetValueF("Conjugate Gradients","max_move", cgMaxMoveFullRelax);
        qmTimeStep = ini.GetValueF("Quick Min","time_step", qmTimeStep);
        qmMaxMove = ini.GetValueF("Quick Min","max_move", qmMaxMove);

        hessianWithinRadius = ini.GetValueF("Hessian", "within_radius", hessianWithinRadius);
        hessianMinDisplacement = ini.GetValueF("Hessian", "min_displacement", hessianMinDisplacement);
        hessianPrefactorMax = ini.GetValueF("Hessian", "prefactor_max", hessianPrefactorMax);
        hessianPrefactorMin = ini.GetValueF("Hessian", "prefactor_min", hessianPrefactorMin);
 
        dimerSeparation = ini.GetValueF("Dimer", "separation", dimerSeparation);
        dimerRotationAngle = ini.GetValueF("Dimer", "finite_diff_angle", dimerRotationAngle);
        dimerRotationsMax = ini.GetValueL("Dimer", "rotations_max", dimerRotationsMax);
        dimerRotationsMin = ini.GetValueL("Dimer", "rotations_min", dimerRotationsMin);
        dimerWindowMax = ini.GetValueF("Dimer", "window_max", dimerWindowMax);
        dimerWindowMin = ini.GetValueF("Dimer", "window_min", dimerWindowMin);

        lanczosConvergence = ini.GetValueF("Lanczos", "convergence", lanczosConvergence);
        lanczosIteration = ini.GetValueL("Lanczos", "iteration", lanczosIteration);

        displaceNSamples = ini.GetValueL("DisplacementSampling", "samples", displaceNSamples);
        displaceIterMax = ini.GetValueL("DisplacementSampling", "iter_max", displaceIterMax);
        displaceTorqueConvergence = ini.GetValueF("DisplacementSampling", "torque_convergence", displaceTorqueConvergence);
        displaceMaxCurvature = ini.GetValueF("DisplacementSampling", "max_curvature", displaceMaxCurvature);
        displaceMaxDE = ini.GetValueF("DisplacementSampling", "max_de", displaceMaxDE);
        displaceCutoffs = ini.GetValue("DisplacementSampling", "cutoffs", displaceCutoffs);
        displaceMagnitudes = ini.GetValue("DisplacementSampling", "magnitudes", displaceMagnitudes);

        mdTimeStep = ini.GetValueF("Dynamics", "time_step", mdTimeStep);
        mdTemperature = ini.GetValueF("Dynamics", "temperature", mdTemperature);
        mdSteps = ini.GetValueL("Dynamics", "steps", mdSteps);
        mdDephaseSteps = ini.GetValueL("Dynamics", "dephase_steps", mdDephaseSteps);
        mdMaxMovedDist = ini.GetValueF("Dynamics", "parrep_max_moved_dist", mdMaxMovedDist);  
        mdRefine = ini.GetValueB("Dynamics", "refine", mdRefine);
        mdAutoStop = ini.GetValueB("Dynamics", "auto_stop", mdAutoStop);
        mdRecordAccuracy = ini.GetValueL("Dynamics", "record_accuracy", mdRecordAccuracy);
        mdRefineAccuracy = ini.GetValueL("Dynamics", "refine_accuracy", mdRefineAccuracy);
        mdCheckFreq = ini.GetValueL("Dynamics", "check_freq", mdCheckFreq);
        mdRelaxSteps = ini.GetValueL("Dynamics", "new_relax_step", mdRelaxSteps);
        mdDephaseLoopStop = ini.GetValueB("Dynamics", "dephase_loop_stop", mdDephaseLoopStop);
        mdDephaseLoopMax = ini.GetValueL("Dynamics", "dephase_loop_max", mdDephaseLoopMax);
         
        string thermostatString;
        thermostatString = ini.GetValue("Dynamics", "thermostat", "andersen");
        thermostatString = toLowerCase(thermostatString);
        if (thermostatString == "andersen") {
            thermostat = Dynamics::ANDERSEN;
        }else if (thermostatString == "nosehover") {
            thermostat = Dynamics::NOSE_HOVER;
        }

        thermoAndersenAlpha = ini.GetValueF("Dynamics","andersen_alpha",thermoAndersenAlpha);
        thermoAndersenTcol = ini.GetValueF("Dynamics","andersen_tcol",thermoAndersenTcol);
        thermoNoseMass = ini.GetValueF("Dynamics","nose_mass",thermoNoseMass);

        bondBoost = ini.GetValueB("Hyperdynamics","bond_boost",bondBoost);
        bondBoostRMDS = ini.GetValueL("Hyperdynamics","RMDS",bondBoostRMDS);
        bondBoostDVMAX = ini.GetValueF("Hyperdynamics","DVMAX",bondBoostDVMAX);
        bondBoostQRR = ini.GetValueF("Hyperdynamics","QRR",bondBoostQRR );
        bondBoostPRR = ini.GetValueF("Hyperdynamics","PRR",bondBoostPRR );
        bondBoostQcut= ini.GetValueF("Hyperdynamics","Qcut",bondBoostQcut);

    }
    else
    {
        fprintf(stderr, "Couldn't parse the ini file. Perhaps you are "
                        "using the old style config?\n");
        error = 1;
    }
    return error;
}

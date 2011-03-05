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
#include "BondBoost.h"
#include "SaddleSearch.h"
#include "ImprovedDimer.h"

Parameters::Parameters(){

    // [Main] //
    job = Job::PROCESS_SEARCH;
    randomSeed = -1;
    potential = "lj";
    temperature = 300.0;

    // [Structure Comparison] //
    distanceDifference = 0.1;
    neighborCutoff = 3.3;
    structureComparisonEnergyDifference = 0.01;
    checkRotation = false;

    // [Debug] //
    writeMovies = false;

    // [Process Search] //
    processSearchMinimizeFirst = false;
    processSearchDefaultPrefactor = 0;
    processSearchMinimizationOffset = 0.2;

    // [Saddle Search] //
    saddleDisplaceType = SaddlePoint::DISP_LOAD;
    saddleMinmodeMethod = SaddlePoint::MINMODE_DIMER;
    saddleMaxStepSize = 0.2;
    saddleMaxEnergy = 20.0;
    saddleMaxIterations = 1000;
    saddleDisplaceRadius = 4.0;
    saddleDisplaceMagnitude = 0.1;
    saddleMaxSingleDisplace = 10.;
//    saddleMaxJumpAttempts = 0; // undocumented
    saddlePerpForceRatio = 0.0; // undocumented
    saddleLocalizationAtoms = 0; // undocumented

    // [Optimizers] //
    optMethod = "cg";
    optMaxIterations = 1000;
    optConvergedForce = 0.005;
    optMaxMove = 0.2;
    optFiniteDiffStep = 0.001;
    optTimeStep = 0.1;

    // [Dimer] //
    dimerSeparation = 0.001;
    dimerRotationAngle = 0.005;
    dimerImproved = false;
    dimerConvergedRotation = 5.0; // degrees
    dimerOptimizer = ImprovedDimer::OPT_SD;
    // old, for comparison with improved dimer method, will be removed later
    dimerTorqueMin = 0.1;
    dimerTorqueMax = 1.0;
    dimerRotationsMin = 1;
    dimerRotationsMax = 8;

    // [Lanczos] //
    lanczosFiniteDiff = 0.001;
    lanczosTolerance = 0.001;
    lanczosMaxIterations = 20;

    // [Hessian] //
    hessianType = Hessian::REACTANT;
    hessianMinDisplacement = 0.25;
    hessianWithinRadius = 5.0;
    hessianPrefactorMax = 10e20;
    hessianPrefactorMin = 10e8;

    // [Parallel Replica] //
    mdTimeStep = 1;
    mdSteps = 1000;
    mdMaxMovedDist = 2.0;
    mdRefine = false;
    mdAutoStop = false;
    mdRecordAccuracy = 1;
    mdRefineAccuracy = 1;
    mdCheckFreq = 500;
    mdRelaxSteps = 500;
    mdDephaseSteps = 200;
    mdDephaseLoopStop = false;
    mdDephaseLoopMax = 5;

    // [Distributed Replica] //
    drBalanceSteps = 500;
    drSamplingSteps = 500;
    drTargetTemperature = 300.0;

    // [Hyperdynamics] //
    biasPotential = Hyperdynamics::NONE;
    bondBoostDVMAX = 0.0;
    bondBoostQRR = 0.0001; // can not be set to 0
    bondBoostPRR = 0.95;
    bondBoostQcut = 3.0;
    bondBoostRMDS = 0;

    // [Thermostat]
    thermostat = Dynamics::ANDERSEN;
    thermoAndersenAlpha = 0.2; // collision strength
    thermoAndersenTcol = 10; // collision frequency in unit of dt
    thermoNoseMass = 1.0;
    thermoLangvinFriction = 0.005;

    // [Displacement Sampling] //
    displaceNSamples = 32; // number of samples to take
    displaceIterMax = 32; // maximum number of rotations to perform on the dimer
    displaceTorqueConvergence = 0.01; // convergence criteria of the dimer rotation
    displaceMaxCurvature = -0.1; // maximum curvature which considered good; avoid shallow but negative curvatures
    displaceMaxDE = 10.0; // maximum dE which is considered good; should use saddleMaxEnergy?
    displaceCutoffs = "0.0 3.3";
    displaceMagnitudes = "0.0625 0.125 0.25";

    // [Basin Hopping] //
    basinHoppingStepSize = 0.1;
    basinHoppingSteps = 10000;
    basinHoppingSingleAtomDisplace = false;
    basinHoppingStayMinimized = false;
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
        // if we succesfully read the file, then parse it as an INI

        // [Main] //

        string jobString;
        jobString = ini.GetValue("Main", "job", "process_search");
        jobString = toLowerCase(jobString);
        if (jobString == "process_search") {
            job = Job::PROCESS_SEARCH;
        }else if (jobString == "saddle_search") {
            job = Job::SADDLE_SEARCH;
        }else if (jobString == "minimization") {
            job = Job::MINIMIZATION;
        }else if (jobString == "point") {
            job = Job::POINT;
        }else if (jobString == "parallel_replica"){
            job = Job::PARALLEL_REPLICA;
        }else if (jobString == "distributed_replica"){
            job = Job::DISTRIBUTED_REPLICA;
        }else if (jobString == "basin_hopping") {
            job = Job::BASIN_HOPPING;
        }else if (jobString == "hessian") {
            job = Job::HESSIAN;
        }else if (jobString == "finite_difference"){
            job = Job::FINITE_DIFFERENCE;
        }else if (jobString == "dimer_rotation"){
            job = Job::DIMER_ROTATION;
        }else if (jobString == "displacement_sampling"){
            job = Job::DISPLACEMENT_SAMPLING;
        }else if (jobString == "test"){
            job = Job::TEST;
        }else{
            fprintf(stderr, "Unknown job_type: %s\n", jobString.c_str());
            error = 1;
        }

        temperature = ini.GetValueF("Main", "temperature", temperature);
        randomSeed = ini.GetValueL("Main", "random_seed", randomSeed);
        // Initialize random generator
        if(randomSeed < 0){
            unsigned i = time(NULL);
            randomSeed = i;
            helper_functions::random(i);
        }else{
            helper_functions::random(randomSeed);
        }

        potential = toLowerCase(ini.GetValue("Main", "potential"));
        
        // [Structure Comparison] //

        distanceDifference = ini.GetValueF("Structure Comparison", "distance_difference", distanceDifference);
        neighborCutoff = ini.GetValueF("Structure Comparison", "neighbor_cutoff", neighborCutoff);
        structureComparisonEnergyDifference = ini.GetValueF("Structure Comparison", "energy_difference", structureComparisonEnergyDifference);
        checkRotation = ini.GetValueB("Structure Comparison", "check_rotation", checkRotation);

        // [Process Search] //

        processSearchMinimizeFirst = ini.GetValueB("Process Search", "minimize_first", processSearchMinimizeFirst);
        processSearchDefaultPrefactor = ini.GetValueF("Process Search", "default_prefactor", processSearchDefaultPrefactor);
        processSearchPrefactorMax = ini.GetValueF("Process Search", "prefactor_max", processSearchPrefactorMax);
        processSearchPrefactorMin = ini.GetValueF("Process Search", "prefactor_min", processSearchPrefactorMin);
        processSearchMinimizationOffset = ini.GetValueF("Process Search", "minimization_offset", processSearchMinimizationOffset);
        // [Saddle Search] //

        string minmodeMethodString = ini.GetValue("Saddle Search", "min_mode_method", "dimer");
        minmodeMethodString = toLowerCase(minmodeMethodString);
        if(minmodeMethodString == "dimer"){
            saddleMinmodeMethod = SaddlePoint::MINMODE_DIMER;
        }else if(minmodeMethodString == "lanczos"){
            saddleMinmodeMethod = SaddlePoint::MINMODE_LANCZOS;
        }else if(minmodeMethodString == "exact") {
            saddleMinmodeMethod = SaddlePoint::MINMODE_EXACT;
        }
        string displaceString = ini.GetValue("Saddle Search", "displace_type", "none");
        displaceString = toLowerCase(displaceString);
//        if(displaceString == "none")
//        { 
//            saddleDisplaceType = SaddlePoint::DISP_NONE;
//        }
//        else if(displaceString == "not_FCC_HCP_coordinated_client")
        if(displaceString == "client_not_fcc_hcp_coordinated")
        {
            saddleDisplaceType = SaddlePoint::DISP_NOT_FCC_OR_HCP;
        }
        else if(displaceString == "client_least_coordinated")
        {
            saddleDisplaceType = SaddlePoint::DISP_MIN_COORDINATED;
        }
        else if(displaceString == "client_last_atom")
        {
            saddleDisplaceType = SaddlePoint::DISP_LAST_ATOM;
        }
        else if(displaceString == "client_random")
        {
            saddleDisplaceType = SaddlePoint::DISP_RANDOM;
        }        
        // default is a displacement made on the server
        else
        { 
            saddleDisplaceType = SaddlePoint::DISP_LOAD;
        }
        
        saddleDisplaceMagnitude = ini.GetValueF("Saddle Search", "displace_magnitude", saddleDisplaceMagnitude);
        saddleDisplaceRadius = ini.GetValueF("Saddle Search", "displace_radius", saddleDisplaceRadius);
        saddleMaxEnergy = ini.GetValueF("Saddle Search", "max_energy", saddleMaxEnergy);
        saddleMaxStepSize = ini.GetValueF("Saddle Search", "max_step_size", saddleMaxStepSize);
        saddleMaxIterations = ini.GetValueL("Saddle Search", "max_iterations", saddleMaxIterations);
//        saddleMaxJumpAttempts = ini.GetValueL("Saddle Search", "max_jump_attempts", saddleMaxJumpAttempts); //undocumented
        saddleMaxSingleDisplace = ini.GetValueF("Saddle Search", "max_single_displace", saddleMaxSingleDisplace);
        saddlePerpForceRatio = ini.GetValueF("Saddle Search", "perp_force_ratio", saddlePerpForceRatio); //undocumented
        saddleLocalizationAtoms = ini.GetValueF("Saddle Search", "localization_atoms", saddleLocalizationAtoms); //undocumented

        // [Optimizers] //

        optMethod = toLowerCase(ini.GetValue("Optimizers", "method", optMethod));
        optConvergedForce = ini.GetValueF("Optimizers", "converged_force", optConvergedForce);
        optMaxIterations = ini.GetValueL("Optimizers", "max_iterations", optMaxIterations);
        optMaxMove = ini.GetValueF("Optimizers","max_move", optMaxMove);
        optFiniteDiffStep = ini.GetValueF("Optimizers","finite_diff_step", optFiniteDiffStep);
        optTimeStep = ini.GetValueF("Optimizers","time_step", optTimeStep);

        // [Dimer] //

        dimerSeparation = ini.GetValueF("Dimer", "separation", dimerSeparation);
        dimerRotationAngle = ini.GetValueF("Dimer", "finite_diff_angle", dimerRotationAngle);
        dimerImproved = ini.GetValueB("Dimer", "improved", dimerImproved);
        dimerConvergedRotation = ini.GetValueF("Dimer", "converged_rotation", dimerConvergedRotation);
        dimerMaxIterations = ini.GetValueL("Dimer", "max_iterations", dimerMaxIterations);
        string dimerOptString;
        dimerOptString = ini.GetValue("Dimer", "optimizer", "steepest_descent");
        dimerOptString = toLowerCase(dimerOptString);
        if (dimerOptString == "steepest_descent") {
            dimerOptimizer = ImprovedDimer::OPT_SD;
        }else if (dimerOptString == "conjugate_gradient") {
            dimerOptimizer = ImprovedDimer::OPT_CG;
        }else if (dimerOptString == "lbfgs") {
            dimerOptimizer = ImprovedDimer::OPT_LBFGS;
        }
        // old, for comparison with Improved dimer method, will be removed later
        dimerRotationsMin = ini.GetValueL("Dimer", "rotations_min", dimerRotationsMin);
        dimerRotationsMax = ini.GetValueL("Dimer", "rotations_max", dimerRotationsMax);
        dimerTorqueMin = ini.GetValueF("Dimer", "torque_min", dimerTorqueMin);
        dimerTorqueMax = ini.GetValueF("Dimer", "torque_max", dimerTorqueMax);

        // [Lanczos] //

        lanczosFiniteDiff = ini.GetValueF("Lanczos", "finite_diff_step", lanczosFiniteDiff);
        lanczosTolerance = ini.GetValueF("Lanczos", "tolerance", lanczosTolerance);
        lanczosMaxIterations = ini.GetValueL("Lanczos", "max_iterations", lanczosMaxIterations);

        // [Hessian] //

        string hessianString = ini.GetValue("Hessian", "type", "reactant");
        hessianString = toLowerCase(hessianString);
        if(hessianString == "reactant"){
            hessianType = Hessian::REACTANT;
        }else if(hessianString == "saddle"){
            hessianType = Hessian::SADDLE;
        }else if(hessianString == "product"){
            hessianType = Hessian::PRODUCT;
        }
        hessianWithinRadius = ini.GetValueF("Hessian", "within_radius", hessianWithinRadius);
        hessianMinDisplacement = ini.GetValueF("Hessian", "min_displacement", hessianMinDisplacement);
 
        // [Displacement Sampling] //

        displaceNSamples = ini.GetValueL("Displacement Sampling", "samples", displaceNSamples);
        displaceIterMax = ini.GetValueL("Displacement Sampling", "iter_max", displaceIterMax);
        displaceTorqueConvergence = ini.GetValueF("Displacement Sampling", "torque_convergence", displaceTorqueConvergence);
        displaceMaxCurvature = ini.GetValueF("Displacement Sampling", "max_curvature", displaceMaxCurvature);
        displaceMaxDE = ini.GetValueF("Displacement Sampling", "max_de", displaceMaxDE);
        displaceCutoffs = ini.GetValue("Displacement Sampling", "cutoffs", displaceCutoffs);
        displaceMagnitudes = ini.GetValue("Displacement Sampling", "magnitudes", displaceMagnitudes);

        mdTimeStep = ini.GetValueF("Dynamics", "time_step", mdTimeStep);
        mdTimeStep = mdTimeStep * 0.09823; //transfer the time unit from fs to 10.18 fs 
        mdSteps = ini.GetValueL("Dynamics", "steps", mdSteps);
        mdDephaseSteps = ini.GetValueL("Dynamics", "dephase_steps", mdDephaseSteps);
        mdRefine = ini.GetValueB("Dynamics", "refine_transition_time", mdRefine);
        mdAutoStop = ini.GetValueB("Dynamics", "auto_stop", mdAutoStop);
        mdRecordAccuracy = ini.GetValueL("Dynamics", "record_resolution", mdRecordAccuracy);
        mdRefineAccuracy = ini.GetValueL("Dynamics", "bisection_accuracy", mdRefineAccuracy);
        mdCheckFreq = ini.GetValueL("Dynamics", "check_period", mdCheckFreq);
        mdRelaxSteps = ini.GetValueL("Dynamics", "post_transition_steps", mdRelaxSteps);
        mdDephaseLoopStop = ini.GetValueB("Dynamics", "dephase_loop_stop", mdDephaseLoopStop);
        mdDephaseLoopMax = ini.GetValueL("Dynamics", "dephase_loop_max", mdDephaseLoopMax);

        // [Distributed Replica] //
        drBalanceSteps = ini.GetValueL("Distributed Replica", "balance_steps", drBalanceSteps);
        drSamplingSteps = ini.GetValueL("Distributed Replica", "sampling_steps", drSamplingSteps);
        drTargetTemperature = ini.GetValueF("Distributed Replica", "target_temperature", drTargetTemperature);

        // [Thermostat] //

        string thermostatString;
        thermostatString = ini.GetValue("Dynamics", "thermostat", "andersen");
        thermostatString = toLowerCase(thermostatString);
        if (thermostatString == "andersen") {
            thermostat = Dynamics::ANDERSEN;
        }else if (thermostatString == "nose_hoover") {
            thermostat = Dynamics::NOSE_HOOVER;
        }else if (thermostatString == "langevin") {
            thermostat = Dynamics::LANGEVIN;
        }
        
        thermoAndersenAlpha = ini.GetValueF("Dynamics","andersen_alpha",thermoAndersenAlpha);
        thermoAndersenTcol = ini.GetValueF("Dynamics","andersen_collision_steps",thermoAndersenTcol);
        thermoNoseMass = ini.GetValueF("Dynamics","nose_mass",thermoNoseMass);
        thermoLangvinFriction = ini.GetValueF("Dynamics","langevin_friction",thermoLangvinFriction);
         
        string hyperString;
        hyperString = ini.GetValue("Hyperdynamics","bias_potential","none");
        hyperString = toLowerCase(hyperString);
        if (hyperString == "none"){
            biasPotential = Hyperdynamics::NONE;
        }else if (hyperString == "bond_boost"){
            biasPotential = Hyperdynamics::BOND_BOOST;
        }

        // [Hyperdynamics] //
        bondBoostRMDS = ini.GetValueL("Hyperdynamics","bb_rmd_steps",bondBoostRMDS);
        bondBoostDVMAX = ini.GetValueF("Hyperdynamics","bb_dvmax",bondBoostDVMAX);
        bondBoostQRR = ini.GetValueF("Hyperdynamics","bb_stretch_threshold",bondBoostQRR );
        bondBoostPRR = ini.GetValueF("Hyperdynamics","bb_ds_curvature",bondBoostPRR );
        bondBoostQcut= ini.GetValueF("Hyperdynamics","bb_rcut",bondBoostQcut);

        // [Basin Hopping] //
        basinHoppingStepSize = ini.GetValueF("Basin Hopping", "step_size", basinHoppingStepSize);
        basinHoppingSteps = ini.GetValueF("Basin Hopping", "steps", basinHoppingSteps);
        basinHoppingSingleAtomDisplace = ini.GetValueB("Basin Hopping", "single_atom_displace", basinHoppingSingleAtomDisplace);
        basinHoppingStayMinimized = ini.GetValueB("Basin Hopping", "stay_minimized", basinHoppingStayMinimized);

        // [Debug] //
        writeMovies= ini.GetValueB("Debug", "write_movies", writeMovies);

    }
    else
    {
        fprintf(stderr, "Couldn't parse the ini file. Perhaps you are "
                        "using the old style config?\n");
        error = 1;
    }
    return error;
}

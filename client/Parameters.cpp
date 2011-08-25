//----------------------------------------------------------------------------------e
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
#include "NudgedElasticBand.h"
#include "Potential.h"

Parameters::Parameters(){

    // [Main] //
    job = Job::PROCESS_SEARCH;
    randomSeed = -1;
    potential = Potential::POT_LJ;
    temperature = 300.0;
    checkpointForceCalls = -1;
    quiet = false;

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
    processSearchPrefactorMax = 10e20;
    processSearchPrefactorMin = 10e8;

    // [Saddle Search] //
    saddleDisplaceType = SaddleSearch::DISP_LOAD;
    saddleMinmodeMethod = SaddleSearch::MINMODE_DIMER;
    saddleMaxStepSize = 0.2;
    saddleMaxEnergy = 20.0;
    saddleMaxIterations = 1000;
    saddleDisplaceRadius = 4.0;
    saddleDisplaceMagnitude = 0.1;
    saddleMaxSingleDisplace = 10.;
    saddleConvergedForce = 0.005;
    saddlePerpForceRatio = 0.0; // undocumented


    // [Optimizers] //
    optMethod = "cg";
    optMaxIterations = 1000;
    optConvergedForce = 0.005;
    optMaxMove = 0.2;
    optFiniteDist = 0.001;
    optTimeStep = 0.1;

    // [Dimer] //
    dimerSeparation = 0.001;
    dimerRotationAngle = 0.005;
    dimerImproved = false;
    dimerConvergedAngle = 5.0; // degrees
    dimerOptMethod = ImprovedDimer::OPT_CG;
    dimerTorqueMin = 0.1; // old dimer
    dimerTorqueMax = 1.0; // old dimer
    dimerRotationsMin = 1; // old dimer
    dimerRotationsMax = 8; // old dimer

    // [Lanczos] //
    lanczosFiniteDist = 0.001;
    lanczosTolerance = 0.001;
    lanczosMaxIterations = 20;

    // [Hessian] //
    hessianType = Hessian::REACTANT;
    hessianFiniteDist = 0.001;
    hessianMinDisplacement = 0.25;
    hessianWithinRadius = 5.0;

    // [Nudged Elastic Band] //
    nebImages = 5;
    nebSpring = 5.0;
    nebClimbingImageMethod = true;
    nebOldTangent = false;
    nebOptMethod = NudgedElasticBand::OPT_CG;

    // [Parallel Replica] //
    mdTimeStep = 1;
    mdSteps = 1000;
    mdMaxMovedDist = 2.0;
    mdRefine = true;
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

    // [Thermostat] //
    thermostat = Dynamics::ANDERSEN;
    thermoAndersenAlpha = 0.2; // collision strength
    thermoAndersenTcol = 10; // collision frequency in unit of dt
    thermoNoseMass = 1.0;
    thermoLangvinFriction = 0.005;

    // [Basin Hopping] //
    basinHoppingMaxDisplacement = 0.5;
    basinHoppingSteps = 10000;
    basinHoppingQuenchingSteps = 0;
    basinHoppingSingleAtomDisplace = false;
    basinHoppingSignificantStructure = false;
    basinHoppingMaxDisplacementAlgorithm = "standard";
    basinHoppingDisplacementDistribution = "uniform";
    basinHoppingSwapProbability = 0.0;
    basinHoppingJumpMax = 10;
    basinHoppingJumpSteps = 0;
    basinHoppingMDFirst = false;
    basinHoppingMDTemp = 300.0;

    // MPI
    MPIPotentialRank = -1;
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
        fprintf(stderr, "error: %s\n", strerror(errno));
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

    if(ini.ReadFile(file))
    {
        // [Main] //

        job = toLowerCase(ini.GetValue("Main", "job"));
        temperature = ini.GetValueF("Main", "temperature", temperature);
        potential = toLowerCase(ini.GetValue("Main", "potential"));
        randomSeed = ini.GetValueL("Main", "random_seed", randomSeed);
        checkpointForceCalls = ini.GetValueL("Main", "checkpoint_force_calls", checkpointForceCalls);
        quiet = ini.GetValueB("Main", "quiet", quiet);

        // Initialize random generator
        if(randomSeed < 0){
            unsigned i = time(NULL);
            randomSeed = i;
            helper_functions::random(i);
        }else{
            helper_functions::random(randomSeed);
        }

        // [Debug] //

        writeMovies= ini.GetValueB("Debug", "write_movies", writeMovies);

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

        saddleMinmodeMethod = toLowerCase(ini.GetValue("Saddle Search", "min_mode_method", "dimer"));
        saddleDisplaceMagnitude = ini.GetValueF("Saddle Search", "displace_magnitude", saddleDisplaceMagnitude);
        saddleDisplaceRadius = ini.GetValueF("Saddle Search", "displace_radius", saddleDisplaceRadius);
        saddleMaxEnergy = ini.GetValueF("Saddle Search", "max_energy", saddleMaxEnergy);
        saddleMaxStepSize = ini.GetValueF("Saddle Search", "max_step_size", saddleMaxStepSize);
        saddleMaxIterations = ini.GetValueL("Saddle Search", "max_iterations", saddleMaxIterations);
        saddleMaxSingleDisplace = ini.GetValueF("Saddle Search", "max_single_displace", saddleMaxSingleDisplace);
        saddleConvergedForce = ini.GetValueF("Saddle Search", "converged_force", saddleConvergedForce);

        saddlePerpForceRatio = ini.GetValueF("Saddle Search", "perp_force_ratio", saddlePerpForceRatio); //undocumented
        saddleDisplaceType = toLowerCase(ini.GetValue("Saddle Search", "client_displace_type", SaddleSearch::DISP_LOAD));
        // XXX: This is a result of mixing our server/client config files.
        if(saddleDisplaceType != SaddleSearch::DISP_NOT_FCC_OR_HCP &&
           saddleDisplaceType != SaddleSearch::DISP_MIN_COORDINATED && 
           saddleDisplaceType != SaddleSearch::DISP_LAST_ATOM && 
           saddleDisplaceType != SaddleSearch::DISP_RANDOM){
                saddleDisplaceType = SaddleSearch::DISP_LOAD;
        }

        // [Optimizers] //

        optMethod = toLowerCase(ini.GetValue("Optimizers", "opt_method", optMethod));
        optConvergedForce = ini.GetValueF("Optimizers", "converged_force", optConvergedForce);
        optMaxIterations = ini.GetValueL("Optimizers", "max_iterations", optMaxIterations);
        optMaxMove = ini.GetValueF("Optimizers","max_move", optMaxMove);
        optFiniteDist = ini.GetValueF("Optimizers","finite_dist", optFiniteDist);
        optTimeStep = ini.GetValueF("Optimizers","time_step", optTimeStep);

        // [Dimer] //

        dimerSeparation = ini.GetValueF("Dimer", "separation", dimerSeparation);
        dimerRotationAngle = ini.GetValueF("Dimer", "finite_angle", dimerRotationAngle);
        dimerImproved = ini.GetValueB("Dimer", "improved", dimerImproved);
        dimerConvergedAngle = ini.GetValueF("Dimer", "converged_angle", dimerConvergedAngle);
        dimerMaxIterations = ini.GetValueL("Dimer", "max_iterations", dimerMaxIterations);
        dimerOptMethod = toLowerCase(ini.GetValue("Dimer", "opt_method", "sd"));
        dimerRotationsMin = ini.GetValueL("Dimer", "rotations_min", dimerRotationsMin); // old
        dimerRotationsMax = ini.GetValueL("Dimer", "rotations_max", dimerRotationsMax); // old & new
        dimerTorqueMin = ini.GetValueF("Dimer", "torque_min", dimerTorqueMin); // old
        dimerTorqueMax = ini.GetValueF("Dimer", "torque_max", dimerTorqueMax); // old

        // [Lanczos] //

        lanczosFiniteDist = ini.GetValueF("Lanczos", "finite_dist", lanczosFiniteDist);
        lanczosTolerance = ini.GetValueF("Lanczos", "tolerance", lanczosTolerance);
        lanczosMaxIterations = ini.GetValueL("Lanczos", "max_iterations", lanczosMaxIterations);

        // [Hessian] //

        hessianType = toLowerCase(ini.GetValue("Hessian", "type", "reactant"));
        hessianFiniteDist = ini.GetValueF("Hessian", "finite_dist", hessianFiniteDist);
        hessianWithinRadius = ini.GetValueF("Hessian", "within_radius", hessianWithinRadius);
        hessianMinDisplacement = ini.GetValueF("Hessian", "min_displacement", hessianMinDisplacement);

        // [Nudged Elastic Band]

        nebImages = ini.GetValueL("NEB", "images", nebImages);
        nebSpring = ini.GetValueF("NEB", "spring", nebSpring);
        nebClimbingImageMethod = ini.GetValueB("NEB", "climbing_image_method", nebClimbingImageMethod);
        nebOldTangent = ini.GetValueB("NEB", "old_tangent", nebOldTangent);
        nebOptMethod = toLowerCase(ini.GetValue("NEB", "opt_method", nebOptMethod));

        // [Molecular Dynamics] //

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

        thermostat = toLowerCase(ini.GetValue("Dynamics", "thermostat", "andersen"));
        thermoAndersenAlpha = ini.GetValueF("Dynamics","andersen_alpha",thermoAndersenAlpha);
        thermoAndersenTcol = ini.GetValueF("Dynamics","andersen_collision_steps",thermoAndersenTcol);
        thermoNoseMass = ini.GetValueF("Dynamics","nose_mass",thermoNoseMass);
        thermoLangvinFriction = ini.GetValueF("Dynamics","langevin_friction",thermoLangvinFriction);

        // [Hyperdynamics] //

        bondBoostRMDS = ini.GetValueL("Hyperdynamics","bb_rmd_steps",bondBoostRMDS);
        bondBoostDVMAX = ini.GetValueF("Hyperdynamics","bb_dvmax",bondBoostDVMAX);
        bondBoostQRR = ini.GetValueF("Hyperdynamics","bb_stretch_threshold",bondBoostQRR );
        bondBoostPRR = ini.GetValueF("Hyperdynamics","bb_ds_curvature",bondBoostPRR );
        bondBoostQcut= ini.GetValueF("Hyperdynamics","bb_rcut",bondBoostQcut);
        biasPotential = toLowerCase(ini.GetValue("Hyperdynamics","bias_potential","none"));

        // [Basin Hopping] //        

        basinHoppingMaxDisplacement = ini.GetValueF("Basin Hopping", "max_displacement", basinHoppingMaxDisplacement);
        basinHoppingSteps = ini.GetValueL("Basin Hopping", "steps", basinHoppingSteps);
	    basinHoppingQuenchingSteps = ini.GetValueL("Basin Hopping", "quenching_steps", basinHoppingQuenchingSteps);
        basinHoppingSingleAtomDisplace = ini.GetValueB("Basin Hopping", "single_atom_displace", basinHoppingSingleAtomDisplace);
        basinHoppingSignificantStructure = ini.GetValueB("Basin Hopping", "significant_structure", basinHoppingSignificantStructure);
	    basinHoppingMaxDisplacementAlgorithm = toLowerCase(ini.GetValue("Basin Hopping", "max_displacement_algorithm", basinHoppingMaxDisplacementAlgorithm));
        basinHoppingDisplacementDistribution = toLowerCase(ini.GetValue("Basin Hopping", "displacement_distribution", basinHoppingDisplacementDistribution));
        basinHoppingSwapProbability = ini.GetValueF("Basin Hopping", "swap_probability", basinHoppingSwapProbability);
        basinHoppingJumpMax = ini.GetValueL("Basin Hopping", "jump_max", basinHoppingJumpMax);
        basinHoppingJumpSteps = ini.GetValueL("Basin Hopping", "jump_steps", basinHoppingJumpSteps);
        basinHoppingMDFirst = ini.GetValueB("Basin Hopping", "md_first", basinHoppingMDFirst);
        basinHoppingMDTemp = ini.GetValueF("Basin Hopping", "md_temp", ini.GetValueF("Main", "temperature", temperature));
	

    }
    else
    {
        fprintf(stderr, "Couldn't parse the ini file. Perhaps you are "
                        "using the old style config?\n");
        error = 1;
    }
    return error;
}

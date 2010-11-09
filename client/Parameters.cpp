/*
 *===============================================
 *  EON Parameters.cpp
 *===============================================
 */

#include <errno.h>
#include "Parameters.h"
#include "Hessian.h"
#include "INIFile.h"

Parameters::Parameters(){

    // Default values

    jobType = PROCESS_SEARCH;
    randomSeed = -1;
    reactantStateTag = 0;
    potentialTag = 1;
    potentialNoTranslation = 0;
    getPrefactorsTag = 0;
    minimizeOnly = 0;
    minimizeBox = 0;
    maxDifferencePos = 0.1;

    /*
    This is now defined in Epicenters, but it should be moved back as a parameter
    neighborCutoff = 3.3; 
    */
    
    // default parameter for relaxation   
    convergedRelax = 0.005;

    // default parameters for process search
    processSearchMinimizeFirst = 0;

    // default parameters for saddle point determination   
    saddleTypePerturbation = 1;
    saddleRefine = false;
    saddleLowestEigenmodeDetermination = 1;
    saddleConverged = 0.025;
    saddleMaxJumpAttempts = 0;
    saddleMaxStepSize = 0.2;
    saddleMaxEnergy = 20.0;
    saddleNormPerturbation = 0.1;
    saddleWithinRadiusPerturbated = 4.0;
    saddleMaxSinglePerturbation = 0.1;
    saddleMaxIterations = 512;
    saddleMaxIterationsConcave = 256;
    saddlePerpendicularForceRatio = 0.0;

    // default parameters for Hessian determination   
    hessianKind = 0;
    hessianMaxSize = 0;
    hessianMinDisplacement = 0.25;
    hessianWithinRadiusDisplaced = 5.0;
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

    //default parameters used by Parallel Repica Dynamics
    mdTimeStep = 0.1;
    mdTemperature = 300.0;
    mdSteps = 1000;
    PRD_MaxMovedDist = 2.0;
    mdRefine = false;
    mdAutoStop = false;
    RefineAccuracy = 20;
    CheckFreq = 500;
    NewRelaxSteps = 500;
    
    //default parameters used by Hyperdynamics
    BondBoost = false ;
    BBDVMAX = 0.0;
    BBQRR = 0.0001; // Can not be set to be 0.0;
    BBPRR = 0.95;
    BBQcut = 3.0;
    BBRMDS = 0;

    //default parameters used by Thermostat
    Andersen_Alpha=0.2; //collision strength
    Andersen_Tcol=10; //collision frequency in unit of dt
    return;

    // Default parameters for the DisplacementSamplingJob.
    displaceNSamples = 32;              // The number of samples to take.
    displaceIterMax = 32;               // The maximum number of rotations to perform on the dimer.
    displaceTorqueConvergence = 0.01;   // The convergence criteria of the dimer rotation.
    displaceMaxCurvature = -0.1;        // The maximum curvature for which a sample is considered good. Used to avoid shallow but negative curvatures.
    displaceMaxDE = 10.0;               // The maximum dE for which a sample is considered good. XXX: Should use saddleMaxEnergy?
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
        fprintf(stderr, "problem loading parameters file:%s\n", strerror(errno));
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
        // if we succesfully read the file, then parse it as an INI
        randomSeed = ini.GetValueL("Default", "RANDOM_SEED", randomSeed);
        potentialTag = ini.GetValueL("Default", "POTENTIAL_TAG", potentialTag);
        potentialNoTranslation = ini.GetValueL("Default", 
                                               "POTENTIAL_NO_TRANSLATION", 
                                               potentialNoTranslation);
        minimizeOnly = ini.GetValueL("Default", "MINIMIZE_ONLY", minimizeOnly);
        minimizeBox = ini.GetValueL("Default", "MINIMIZE_BOX", minimizeBox);
        convergedRelax = ini.GetValueF("Default", "CONVERGED_RELAX", 
                                       convergedRelax);
        maximumIterations = ini.GetValueL("Default", "MAXIMUM_ITERATIONS",
                                          maximumIterations);
        maxDifferencePos = ini.GetValueF("Default", "MAX_DIFFERENCE_POS", maxDifferencePos);

        string jobTypeString;
        jobTypeString = ini.GetValue("Default", "JOB_TYPE", "processsearch");
        jobTypeString = toLowerCase(jobTypeString);

        if (jobTypeString == "processsearch") {
            jobType = PROCESS_SEARCH;
        }else if (jobTypeString == "saddlesearch") {
            jobType = SADDLE_SEARCH;
        }else if (jobTypeString == "minimization") {
            jobType = MINIMIZATION;
        }else if (jobTypeString == "hessian") {
            jobType = HESSIAN;
		}else if (jobTypeString == "parallelreplica"){
			jobType = PARALLEL_REPLICA;
        }else if (jobTypeString == "replicaexchange"){
            jobType = REPLICA_EXCHANGE;         
		}else if (jobTypeString == "dimerdr"){
			jobType = DIMER_DR;
		}else if (jobTypeString == "dimerrotation"){
			jobType = DIMER_ROTATION;
		}else if (jobTypeString == "displacementsampling"){
			jobType = DISPLACEMENT_SAMPLING;
		}else if (jobTypeString == "test"){
				jobType = TEST;
        }
        else{
            fprintf(stderr, "Unknown JOB_TYPE: %s\n", jobTypeString.c_str());
            error = 1;
        }

        processSearchMinimizeFirst = ini.GetValueL("ProcessSearch", "minimize_first",
                                                   processSearchMinimizeFirst);
        
        saddleTypePerturbation = ini.GetValueL("Saddle_Point", "TYPE_PERTURBATION",
                                               saddleTypePerturbation);
        saddleLowestEigenmodeDetermination = 
            ini.GetValueL("Saddle_Point", "LOWEST_EIGENMODE_DETERMINATION", 
                          saddleLowestEigenmodeDetermination);
        saddleRefine = ini.GetValueB("Saddle_Point", "REFINE", saddleRefine);
        saddleConverged = ini.GetValueF("Saddle_Point", "CONVERGED",
                                        saddleConverged);
        saddleMaxJumpAttempts = ini.GetValueL("Saddle_Point", "MAX_JUMP_ATTEMPTS",
                                              saddleMaxJumpAttempts);
        saddleMaxStepSize = ini.GetValueF("Saddle_Point", "MAX_STEP_SIZE",
                                          saddleMaxStepSize);
        saddleMaxEnergy = ini.GetValueF("Saddle_Point", "MAX_ENERGY", 
                                        saddleMaxEnergy);
        saddleNormPerturbation = ini.GetValueF("Saddle_Point", "NORM_PERTURBATION",
                                               saddleNormPerturbation);
        saddleMaxSinglePerturbation = ini.GetValueF("Saddle_Point",
                                                    "MAX_SINGLE_PERTURBATION",
                                                    saddleMaxSinglePerturbation);
        saddleWithinRadiusPerturbated = ini.GetValueF("Saddle_Point",
                                                      "WITHIN_RADIUS_PERTURBATED",
                                                      saddleWithinRadiusPerturbated);
        saddleMaxIterations = ini.GetValueL("Saddle_Point", 
                                            "MAX_ITERATIONS", 
                                            saddleMaxIterations);
        saddlePerpendicularForceRatio = ini.GetValueF("Default",
                                                      "PERPENDICULAR_FORCE_RATIO",
                                                      saddlePerpendicularForceRatio);

        string hessianType = ini.GetValue("Hessian", "Type", "reactant");
        hessianType = toLowerCase(hessianType);
        if(hessianType == "reactant")
        {
            hessianKind = Hessian::REACTANT;
        }
        else if(hessianType == "saddle")
        {
            hessianKind = Hessian::SADDLE;
        }
        else if(hessianType == "product")
        {
            hessianKind = Hessian::PRODUCT;
        }



        hessianMaxSize = ini.GetValueL("Hessian", "MAX_SIZE", hessianMaxSize);
        hessianWithinRadiusDisplaced = ini.GetValueF("Hessian", 
                                                     "WITHIN_RADIUS_DISPLACED",
                                                     hessianWithinRadiusDisplaced);
        hessianMinDisplacement = ini.GetValueF("Hessian", 
                                               "MIN_DISPLACEMENT",
                                               hessianMinDisplacement);
        
        dimerRotationsHigh = ini.GetValueL("Dimer", "ROTATIONS_HIGH", dimerRotationsHigh);
        dimerRotationsLow = ini.GetValueL("Dimer", "ROTATIONS_LOW", dimerRotationsLow);
        dimerWindowHigh = ini.GetValueF("Dimer", "WINDOW_HIGH", dimerWindowHigh);
        dimerWindowLow = ini.GetValueF("Dimer", "WINDOW_LOW", dimerWindowLow);
        dimerSeparation = ini.GetValueF("Dimer", "SEPARATION", dimerSeparation);
        dimerRotationAngle = ini.GetValueF("Dimer", "ANGLE", dimerRotationAngle);

        displaceNSamples = ini.GetValueL("DisplacementSampling", "NSAMPLES", displaceNSamples);
        displaceIterMax = ini.GetValueL("DisplacementSampling", "ITERMAX", displaceIterMax);
        displaceTorqueConvergence = ini.GetValueF("DisplacementSampling", "TORQUE_CONVERGENCE", displaceTorqueConvergence);
        displaceMaxCurvature = ini.GetValueF("DisplacementSampling", "MAX_CURVATURE", displaceMaxCurvature);
        displaceMaxDE = ini.GetValueF("DisplacementSampling", "MAX_DE", displaceMaxDE);
        displaceCutoffs = ini.GetValue("DisplacementSampling", "CUTOFFS", displaceCutoffs);
        displaceMagnitudes = ini.GetValue("DisplacementSampling", "MAGNITUDES", displaceMagnitudes);

		mdTimeStep = ini.GetValueF("Dynamics","TIMESTEP",mdTimeStep);
 		mdTemperature = ini.GetValueF("Dynamics","TEMPERATURE",mdTemperature);
		mdSteps = ini.GetValueL("Dynamics","STEPS",mdSteps);
		PRD_MaxMovedDist = ini.GetValueF("Dynamics","PRD_MaxMovedDist",PRD_MaxMovedDist);  
		mdRefine = ini.GetValueB("Dynamics","mdRefine",mdRefine);
        mdAutoStop = ini.GetValueB("Dynamics","mdAutoStop",mdAutoStop);
		RefineAccuracy = ini.GetValueL("Dynamics","RefineAccuracy",RefineAccuracy);
		CheckFreq = ini.GetValueL("Dynamics","CheckFreq",CheckFreq);

        lanczosConvergence = ini.GetValueF("Lanczos", "CONVERGENCE", lanczosConvergence);
        lanczosIteration = ini.GetValueL("Lanczos", "ITERATION", lanczosIteration);

        NewRelaxSteps = ini.GetValueL("Dynamics","NewRelaxStep",NewRelaxSteps);

        BondBoost = ini.GetValueB("Hyper","BondBoost",BondBoost);
        BBRMDS = ini.GetValueL("Hyper","RMDS",BBRMDS);
        BBDVMAX = ini.GetValueF("Hyper","DVMAX",BBDVMAX);
        BBQRR = ini.GetValueF("Hyper","QRR",BBQRR );
        BBPRR = ini.GetValueF("Hyper","PRR",BBPRR );
        BBQcut= ini.GetValueF("Hyper","Qcut",BBQcut);
        
        
        cgCurvatureStep = ini.GetValueF("CG","CURVATURE_STEP", cgCurvatureStep);


        Andersen_Alpha = ini.GetValueF("Thermo","ANDERSEN_ALPHA",Andersen_Alpha);
        Andersen_Tcol = ini.GetValueF("Thermo","ANDERSEN_TCOL",Andersen_Tcol);

    }
    else
    {
        fprintf(stderr, "Couldn't parse the ini file. Perhaps you are "
                        "using the old style config?\n");
        error = 1;
    }
    return error;
}

//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "FiniteDifferenceJob.h"
#include "Matter.h"
#include "Constants.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"

using namespace helper_functions;


FiniteDifferenceJob::FiniteDifferenceJob(Parameters *params)
{
    parameters = params;
}

FiniteDifferenceJob::~FiniteDifferenceJob(){ }

std::vector<std::string> FiniteDifferenceJob::run(void)
{
    // No bundling for this job, so bundleNumber is ignored.
    
    // Load the displacement con file and get the position.
    Matter *reactant = new Matter(parameters);
    reactant->con2matter("reactant_passed.con");
    Matrix<double, Eigen::Dynamic, 3> posA;    
    posA = reactant->getPositions();        

    double dRs[] = {0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 1.0, 1.5, 2.0};
    int ndRs = 16;

    Matrix<double, Eigen::Dynamic, 3> forceA;    
    forceA = reactant->getForces();
    
    // Create a random displacement.
    long epicenter = EpiCenters::minCoordinatedEpiCenter(reactant,parameters->neighborCutoff);
    Matrix<double, Eigen::Dynamic, 3> displacement;    
    displacement.resize(reactant->numberOfAtoms(), 3);
    displacement.setZero();
    for(int i = 0; i < reactant->numberOfAtoms(); i++)
    {
        if(reactant->distance(epicenter, i) <= 3.3)
        {
            for(int j = 0; j < 3; j++)
            {
                if(!reactant->getFixed(i))
                {
                    displacement(i, j) = randomDouble(1.0);
                }
            }
        }
    }
    displacement.normalize(); 
    
    // Loop over values of dimer dR and print the output to results.dat.
    FILE *results = fopen("results.dat", "w");
    fprintf(results, "%14s    %14s\n", "dR", "curvature");
    Matrix<double, Eigen::Dynamic, 3> posB;    
    Matrix<double, Eigen::Dynamic, 3> forceB;  
    double curvature = 0.0;  
    for (int dRi = 0; dRi < ndRs; dRi++)
    {
        posB = posA + displacement * dRs[dRi];
        reactant->setPositions(posB);
        forceB = reactant->getForces();
        curvature = (forceB.sum() - forceA.sum()) / dRs[dRi];        
        fprintf(results, "%14.8f    %14.8f\n", dRs[dRi], curvature);
    }
    fclose(results);

    std::vector<std::string> empty;
    return empty;
}



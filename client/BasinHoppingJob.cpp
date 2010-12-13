//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "BasinHoppingJob.h"
#include "Constants.h"
#include "ConjugateGradients.h"
#include "false_boinc.h"
#include "Potentials.h"
#include "HelperFunctions.h"

#include <stdio.h>
#include <string>

using namespace std;
using namespace helper_functions;

BasinHoppingJob::BasinHoppingJob (Parameters *params)
{
    parameters = params;
}

BasinHoppingJob::~BasinHoppingJob()
{}

Matrix<double, Eigen::Dynamic, 3> BasinHoppingJob::displaceRandom()
{
    // Create a random displacement.
    Matrix<double, Eigen::Dynamic, 3> displacement;        
    displacement.resize(trial->numberOfAtoms(), 3);
    displacement.setZero();

    for(int i = 0; i < trial->numberOfAtoms(); i++)
    {
        for(int j = 0; j < 3; j++)
        {
            if(!trial->getFixed(i))
            {
                displacement(i, j) = gaussRandom(0.0, parameters->basinHoppingStepSize);
            }
        }
    }
    return displacement;
}

Matrix<double, Eigen::Dynamic, 3> BasinHoppingJob::displaceSingle()
{
    // Create a random displacement.
    Matrix<double, Eigen::Dynamic, 3> displacement;        
    displacement.resize(trial->numberOfAtoms(), 3);
    displacement.setZero();
    
    long ra = (long)randomDouble(trial->numberOfAtoms());
    while(trial->getFixed(ra))
    {
        ra = (long)random(trial->numberOfAtoms());
    }

    for(int j = 0; j < 3; j++)
    {
        displacement(ra, j) = gaussRandom(0.0, parameters->basinHoppingStepSize);
    }
    
    return displacement;
}

void BasinHoppingJob::run(int bundleNumber)
{
    current = new Matter(parameters);
    trial = new Matter(parameters);
    Matter *tmpMatter = new Matter(parameters);

    current->con2matter("reactant_passed.con");
    *trial = *current;
    *tmpMatter = *current;

    ConjugateGradients cgMin(tmpMatter, parameters);
    cgMin.setOutput(0);
    cgMin.fullRelax();

    double currentEnergy = tmpMatter->getPotentialEnergy();

    double minimumEnergy = currentEnergy;
    Matter *minimumEnergyStructure = new Matter(parameters);
    *minimumEnergyStructure = *current;

    for (int step=0; step<parameters->basinHoppingSteps; step++)
    {

        Matrix<double, Eigen::Dynamic, 3> displacement;
        if(parameters->basinHoppingSingleAtomDisplace)
        {
            displacement = displaceSingle();
        }
        else
        {
            displacement = displaceRandom();
        }
        
        trial->setPositions(current->getPositions() + displacement);

        *tmpMatter = *trial;
        ConjugateGradients cgMin(tmpMatter, parameters);
        cgMin.setOutput(0);
        cgMin.fullRelax();

        double deltaE = tmpMatter->getPotentialEnergy()-currentEnergy;

        double p = exp(-deltaE / (parameters->temperature*8.617343e-5));

        if (randomDouble(1.0)<min(1.0, p)) 
        {
            *current = *trial;
            if(parameters->basinHoppingStayMinimized)
            {
                *current = *tmpMatter;
            }
            currentEnergy = tmpMatter->getPotentialEnergy();
            if (abs(deltaE)>parameters->structureComparisonEnergyDifference) 
            {
                *current = *tmpMatter;
                if (currentEnergy < minimumEnergy) 
                {
                    minimumEnergy = currentEnergy;
                    *minimumEnergyStructure = *current;
                }
            }
            tmpMatter->matter2xyz("movie", true);
        }
        printf("step: %10i energy: %10.8f c_energy: %10.8f min_fc: %ld\n",
               step, currentEnergy, current->getPotentialEnergy(),
               cgMin.totalForceCalls);
    }

    /* Save Results */

    FILE *fileResults, *fileProduct;

    fileResults = fopen("results_0.dat", "wb");

    fprintf(fileResults, "%d termination_reason\n", 0);
    fprintf(fileResults, "%e minimum_energy\n", minimumEnergy);
    fclose(fileResults);

    fileProduct = fopen("product_0.con", "wb");
    minimumEnergyStructure->matter2con(fileProduct);
    fclose(fileProduct);
  
    return;
}

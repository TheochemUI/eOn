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

    for (int step=0; step<parameters->basinHoppingSteps; step++)
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
        
        trial->setPositions(current->getPositions() + displacement);

        *tmpMatter = *trial;
        ConjugateGradients cgMin(tmpMatter, parameters);
        cgMin.setOutput(0);
        cgMin.fullRelax();

        double deltaE = tmpMatter->getPotentialEnergy()-currentEnergy;

        double p = exp(-deltaE / (parameters->temperature*8.617343e-5));

        if (randomDouble(1.0)<min(1.0, p)) {
            *current = *trial;
            currentEnergy = tmpMatter->getPotentialEnergy();
            tmpMatter->matter2xyz("movie", true);
        }
        printf("step: %10i energy: %10.8f c_energy: %10.8f min_fc: %ld\n",
               step, currentEnergy, current->getPotentialEnergy(),
               cgMin.totalForceCalls);
            
    }
}

void BasinHoppingJob::saveData(int status, int bundleNumber)
{

}

//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "MinimizationJob.h"
#include "ConjugateGradients.h"
#include "Matter.h"
#include "Constants.h"

MinimizationJob::MinimizationJob(Parameters *params)
{
    parameters = params;
}

MinimizationJob::~MinimizationJob(){ }

void MinimizationJob::run(int bundleNumber)
{
    char buff[STRING_SIZE];
    string reactant_passed("reactant_passed");
    string reactant_output("reactant");

    if (bundleNumber < 0) {
        reactant_passed += ".con";
        reactant_output += ".con";
    }else{
        snprintf(buff, STRING_SIZE, "_%i.con", bundleNumber);
        reactant_passed += buff;
        reactant_output += buff;
    }

    Matter *reactant = new Matter(parameters);
    reactant->con2matter(reactant_passed);

    printf("\nBeginning minimization of %s\n", reactant_passed.c_str());

    ConjugateGradients cgMin(reactant, parameters);
    cgMin.setOutput(1);
    cgMin.fullRelax();

    if (cgMin.isItConverged(parameters->convergedRelax)) {
        printf("Minimization converged within tolerence\n");
    }else{
        printf("Minimization did not converge to tolerence!\n"
               "Maybe try to increase maximum_iterations?\n");
    }

    printf("Saving result to %s\n", reactant_output.c_str());
    reactant->matter2con(reactant_output);
    printf("Final Energy: %f\n", reactant->getPotentialEnergy());
}

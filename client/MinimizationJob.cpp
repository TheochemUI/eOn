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
#include "Minimizer.h"
#include "ConjugateGradients.h"
#include "Quickmin.h"
#include "QuickminBox.h"
#include "Matter.h"
#include "Constants.h"

MinimizationJob::MinimizationJob(Parameters *params)
{
    parameters = params;
}

MinimizationJob::~MinimizationJob(){ }

std::vector<std::string> MinimizationJob::run(void)
{
    string reactant_passed("reactant_passed.con");
    string reactant_output("reactant.con");

    std::vector<std::string> returnFiles;
    returnFiles.push_back(reactant_output);

    Matter *reactant = new Matter(parameters);
    reactant->con2matter(reactant_passed);

    printf("\nBeginning minimization of %s\n", reactant_passed.c_str());

    Minimizer* mizer=NULL;
    if(parameters->optMethod == "cg")
    {
        mizer = new ConjugateGradients(reactant, parameters);
    }
    else if(parameters->optMethod == "qm")
    {
        mizer = new Quickmin(reactant, parameters);
    }
    else if(parameters->optMethod == "box")
    {
        mizer = new QuickminBox(reactant, parameters);
    }else{
        printf("Unknown optMethod: %s\n", parameters->optMethod.c_str());
    }

    mizer->setOutput(1);
    mizer->fullRelax();
    
    if (mizer->isItConverged(parameters->optConvergedForce)) {
        printf("Minimization converged within tolerence\n");
    }else{
        printf("Minimization did not converge to tolerence!\n"
               "Maybe try to increase maximum_iterations?\n");
    }

    printf("Saving result to %s\n", reactant_output.c_str());
    reactant->matter2con(reactant_output);
    printf("Final Energy: %f\n", reactant->getPotentialEnergy());

    return returnFiles;
}

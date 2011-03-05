//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "PointJob.h"
#include "Matter.h"
#include "Constants.h"

PointJob::PointJob(Parameters *params)
{
    parameters = params;
}

PointJob::~PointJob(){ }

std::vector<std::string> PointJob::run(void)
{
    std::vector<std::string> returnFiles;
    string reactant_passed("reactant_passed.con");
    string data_out("results.dat");
    returnFiles.push_back(data_out);

    Matter *reactant = new Matter(parameters);
    reactant->con2matter(reactant_passed);

    printf("Energy: %f\n", reactant->getPotentialEnergy());

    FILE *fileResults = fopen(data_out.c_str(), "wb");
    fprintf(fileResults, "%f Energy\n", reactant->getPotentialEnergy());
    fclose(fileResults);

    return returnFiles;
}

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

void PointJob::run(int bundleNumber)
{
    char buff[STRING_SIZE];
    string reactant_passed("reactant_passed");
    string data_out("result");

    if (bundleNumber < 0) {
        reactant_passed += ".con";
        data_out += ".dat";
    }else{
        snprintf(buff, STRING_SIZE, "_%i.con", bundleNumber);
        reactant_passed += buff;
        snprintf(buff, STRING_SIZE, "_%i.dat", bundleNumber);
        data_out += buff;
    }

    Matter *reactant = new Matter(parameters);
    reactant->con2matter(reactant_passed);

    printf("Energy: %f\n", reactant->getPotentialEnergy());

    FILE *fileResults = fopen(data_out.c_str(), "wb");
    fprintf(fileResults, "Energy: %f\n", reactant->getPotentialEnergy());
    fclose(fileResults);
}

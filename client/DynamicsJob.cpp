//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <stdio.h>
#include <string>

#include "Dynamics.h"
#include "DynamicsJob.h"
#include "Parameters.h"
#include "Potential.h"
#include "HelperFunctions.h"

#ifdef BOINC
    #include <boinc/boinc_api.h>
    #include <boinc/diagnostics.h>
    #include <boinc/filesys.h>
#ifdef WIN32
    #include <boinc/boinc_win.h>
    #include <boinc/win_util.h>
#endif
#else
    #include "false_boinc.h"
#endif

using namespace std;
using namespace helper_functions;

DynamicsJob::DynamicsJob(Parameters *params)
{
    parameters = params;
}

DynamicsJob::~DynamicsJob() {}

std::vector<std::string> DynamicsJob::run(void)
{
    Matter *R = new Matter(parameters);
    Matter *F = new Matter(parameters);
    R->con2matter("pos.con");
    *F = *R;

    Dynamics *d = new Dynamics(R, parameters);
    d->run();

    *F = *R;
    FILE *fileProduct;
    std::string productFilename("final.con");
    returnFiles.push_back(productFilename);

    fileProduct = fopen(productFilename.c_str(), "wb");
    F->matter2con(fileProduct);
    fclose(fileProduct);

    delete R;
    delete F;
    delete d;

    std::vector<std::string> returnFiles;
    return returnFiles;
}


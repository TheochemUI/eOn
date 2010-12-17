//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <cstdlib>

#include "Matter.h"
#include "Dynamics.h"
#include "DistributedReplicaJob.h"
#include "ConjugateGradients.h"

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

DistributedReplicaJob::DistributedReplicaJob(Parameters *params)
{
    parameters = params;
}

DistributedReplicaJob::~DistributedReplicaJob()
{ 

}

void DistributedReplicaJob::run(int bundleNumber)
{
    char buff[STRING_SIZE];
    string reactant_passed("reactant_passed");

    if (bundleNumber < 0) {
        reactant_passed += ".con";
    }else{
        snprintf(buff, STRING_SIZE, "_%i.con", bundleNumber);
        reactant_passed += buff;
    }

    reactant = new Matter(parameters);
    reactant->con2matter(reactant_passed);


    printf("Now running Distributed Replica Simulation\n\n");

//    dynamics();

//    saveData(newstate,bundleNumber);
    
    delete reactant;
}


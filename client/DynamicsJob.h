//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef DYNAMICSJOB_H
#define DYNAMICSJOB_H

#include "Matter.h"
#include "Parameters.h"
#include "Job.h"

class DynamicsJob : public Job 
{

    public:

        DynamicsJob(Parameters *params);
        ~DynamicsJob(void);
        std::vector<std::string> run(void);
        Parameters *parameters;
        std::vector<std::string> returnFiles;
};

#endif

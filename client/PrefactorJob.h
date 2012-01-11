//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef PREFACTORJOB_H
#define PREFACTORJOB_H

#include "Job.h"
#include "Parameters.h"

class PrefactorJob : public Job
{
public:
    PrefactorJob(Parameters *params);
    ~PrefactorJob();
    std::vector<std::string> run(void);
    static const char PREFACTOR_REACTANT[];
    static const char PREFACTOR_SADDLE[];
    static const char PREFACTOR_PRODUCT[];
private:
    Parameters *parameters;
};

#endif

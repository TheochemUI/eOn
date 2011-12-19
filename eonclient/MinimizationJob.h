//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef MINIMAZATIONJOB_H
#define MINIMAZATIONJOB_H

#include "Job.h"
#include "Parameters.h"

class MinimizationJob: public Job {
    public:
        MinimizationJob(Parameters *params);
        ~MinimizationJob(void);
        std::vector<std::string> run(void);
    private:
        Parameters *parameters;
        int fcalls;
        enum {
            STATUS_GOOD,
            STATUS_MAX_ITERATIONS,
            STATUS_POTENTIAL_FAILED
        };
};

#endif

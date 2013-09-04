//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef PARALLELREPLICAJOB_H
#define PARALLELREPLICAJOB_H

#include "Job.h"
#include "Parameters.h"

class ParallelReplicaJob: public Job
{
    public:

        ParallelReplicaJob(Parameters *params);
        ~ParallelReplicaJob(void);
        std::vector<std::string> run(void);

    private:
        Parameters *parameters;
        std::vector<std::string> returnFiles;
        Matter *reactant;

        void dephase(Matter *trajectory);
        int refineTransition(std::vector<Matter*> MDSnapshots, bool fake=false);
};

#endif

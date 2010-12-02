//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "Job.h"
#include "Parameters.h"

class ParallelReplicaJob: public Job
{

    public:
        ParallelReplicaJob(Parameters *params);
        ~ParallelReplicaJob(void);
        void run(int bundleNumber);

    private:
        void dynamics();
        bool checkState(Matter *matter);
        bool checkState_nq(Matter *matter);
        void Refine(Matter *mdbuff[]);
        void saveData(int status,int bundleNumber);
        void dephase();
        Parameters *parameters;
        Matter *reactant;
        Matter *min1;
        Matter *min2;
        Matter *transition;
        double SPtime; 
        double RLtime;
        double *SPtimebuff;
        long min_fcalls;
        long md_fcalls;
        long nsteps;
        long nsteps_refined;
        long check_steps;
        long relax_steps;
        bool newstate;
};

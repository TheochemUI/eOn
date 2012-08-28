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

        enum{// DONT CHANGE THE ORDER OF THIS LIST
            STATUS_NEWSTATE, //0
            STATUS_NEWSTATE_CORR, //1
            STATUS_NEWSTATE_OVERFC, //2
            STATUS_TRAN_NOTIME, //3
            STATUS_TRAN_RECROSS, //4
            STATUS_NOTRAN, //5
            STATUS_BAD_RELAXFAILED, //6
            STATUS_BAD_REFINEFAILED, //7
            STATUS_BAD_INFTEMP, //8
        };

        int dynamics();
        long refine(Matter *mdBuffer[], long length, Matter *reactant);
        bool checkState(Matter *current, Matter *reactant);
        void saveData(int state);
        void dephase();
        void printEndStatus();

        Parameters *parameters;

        Matter *current;
        Matter *reactant;
        Matter *transition;
        Matter *transition_relaxed;
        Matter *product;
        Matter *product_relaxed;

        bool newStateFlag;
        bool relaxStatus;
        int jobStatus;

        long minimizeFCalls;
        long mdFCalls;
        long dephaseFCalls;
        long refineFCalls;
        long corrFCalls;

        long transitionStep;

        double time;
        double transitionTime;
        double corrTime;
        double *timeBuffer;

        std::vector<std::string> returnFiles;
};

#endif

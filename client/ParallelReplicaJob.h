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

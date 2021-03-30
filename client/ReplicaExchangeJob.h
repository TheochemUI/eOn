
#ifndef REPLICAEXCHANGEJOB_H
#define REPLICAEXCHANGEJOB_H

#include "Job.h"
#include "Parameters.h"
#include "Dynamics.h"
#include "Eigen.h"

class ReplicaExchangeJob: public Job
{
    public:

        ReplicaExchangeJob(Parameters *params);
        ~ReplicaExchangeJob(void);
        std::vector<std::string> run(void);

    private:

        void saveData();

        Parameters *parameters;
        long forceCalls;
//        Matter **replica;
        Matter *pos;
//        Dynamics **replicaDynamics;
//        double *replicaTemperature;
        std::vector<std::string> returnFiles;
};

#endif

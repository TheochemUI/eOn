//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

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
        Matter **replica;
        Matter *pos;
        Dynamics **replicaDynamics;
        double *replicaTemperature;
        std::vector<std::string> returnFiles;
};

#endif

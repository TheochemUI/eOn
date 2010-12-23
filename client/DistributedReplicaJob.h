//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef DISTRIBUTEDREPLICAJOB_H
#define DISTRIBUTEDREPLICAJOB_H

#include "Job.h"
#include "Parameters.h"
#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

class DistributedReplicaJob: public Job
{

    public:
        DistributedReplicaJob(Parameters *params);
        ~DistributedReplicaJob(void);
        void run(int bundleNumber);
        void balanceStep();
        void samplingStep();
        void saveData(int bundleNumber);
        bool checkState(Matter *matter);

    private:
        Parameters *parameters;
        Matter *reactant;
        Matter *final;
        Matter *min1;
        Matter *min2;
        long bl_fcalls;
        long sp_fcalls;
        long min_fcalls;
        long rf_fcalls;
        bool save_refine;
        double temperature;
};

#endif

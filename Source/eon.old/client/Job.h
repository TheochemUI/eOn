//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
//-----------------------------------------------------------------------------------
#ifndef JOB_H
#define JOB_H
class Job { 
    public:
        enum{
            PROCESS_SEARCH,
            SADDLE_SEARCH,
            MINIMIZATION,
            PARALLEL_REPLICA,
            DISTRIBUTED_REPLICA,
            BASIN_HOPPING,
            HESSIAN,
            FINITE_DIFFERENCE,
            DIMER_ROTATION,
            DISPLACEMENT_SAMPLING,
            TEST
        };
        virtual ~Job() {}
        virtual void run(int bundleNumber)=0;
};
#endif

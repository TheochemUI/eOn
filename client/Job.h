//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef JOB_H
#define JOB_H
#include "Parameters.h"
#include <vector>
#include <string>

class Job { 
    public:

        virtual ~Job() {}
        virtual std::vector<std::string> run()=0;

        static const char PROCESS_SEARCH[];
        static const char SADDLE_SEARCH[];
        static const char MINIMIZATION[];
        static const char POINT[];
        static const char PARALLEL_REPLICA[];
        static const char DISTRIBUTED_REPLICA[];
        static const char BASIN_HOPPING[];
        static const char HESSIAN[];
        static const char FINITE_DIFFERENCE[];
        static const char NUDGED_ELASTIC_BAND[];
        static const char DYNAMICS[];
        static const char TEST[];

        static Job *getJob(Parameters *parameters);
};
#endif

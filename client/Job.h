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

        static Job *getJob(Parameters *parameters);
};
#endif

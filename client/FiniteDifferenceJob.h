//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef FINITEDIFFERENCE_H
#define FINITEDIFFERENCE_H


#include "Job.h"
#include "Parameters.h"
#include "Eigen.h"

class FiniteDifferenceJob: public Job {
    public:
        FiniteDifferenceJob(Parameters *params);
        ~FiniteDifferenceJob(void);
        std::vector<std::string> run(void);
    private:
        Parameters *parameters;
};

#endif

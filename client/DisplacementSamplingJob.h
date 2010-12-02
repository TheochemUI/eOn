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
#include "Eigen/Eigen"

USING_PART_OF_NAMESPACE_EIGEN

class DisplacementSamplingJob: public Job {
    public:
        DisplacementSamplingJob(Parameters *params);
        ~DisplacementSamplingJob(void);
        void run(int bundleNumber);
    private:
        Parameters *parameters;
};

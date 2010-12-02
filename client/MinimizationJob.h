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

class MinimizationJob: public Job {
    public:
        MinimizationJob(Parameters *params);
        ~MinimizationJob(void);
        void run(int bundleNumber);
    private:
        Parameters *parameters;
};

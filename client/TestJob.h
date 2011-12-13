//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef TESTJOB_H
#define TESTJOB_H

#include "Job.h"
#include "Parameters.h"
#include "ConjugateGradients.h"


class TestJob: public Job {
    public:
        TestJob(Parameters *params);
        ~TestJob(void);
        std::vector<std::string> run(void);
    private:
        Parameters *parameters;
        double tolerance;
        void checkFullSearch(void);
        void checkMode(void);
        void checkPotentials(void);
        double getEnergyDiff(string potTag, double refEnergy);
        double getForceDiff(string potTag, double refForce);
};

#endif

//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef SCALEDPOTJOB_H
#define SCALEDPOTJOB_H

#include "Job.h"
#include "Parameters.h"

class ScaledPotJob: public Job
{
    public:

        ScaledPotJob(Parameters *params);
        ~ScaledPotJob(void);
        std::vector<std::string> run(void);

    private:

        int dynamics();
        long refine(Matter *mdBuffer[], long length, Matter *reactant);
        bool checkState(Matter *current, Matter *reactant);
        void saveData(int status);
        void dephase();

        Parameters *parameters;

        Matter *current;
        Matter *reactant;
        Matter *saddle;
        Matter *final;
        Matter *final_tmp;
        Matter *product;

        bool metaStateFlag;
        bool newStateFlag;

        long minimizeFCalls;
        long mdFCalls;
        long dephaseFCalls;
        long refineFCalls;

        long transitionStep;

        double time;
        double scale;
        double minCorrectedTime;
        double transitionTime;
        double transitionPot;
        double *timeBuffer;
        double *potBuffer;

        std::vector<std::string> returnFiles;
};

#endif

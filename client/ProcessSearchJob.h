//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef PROCESSSEARCHJOB_H
#define PROCESSSEARCHJOB_H

#include "Matter.h"
#include "Parameters.h"
#include "Job.h"
#include "MinModeSaddleSearch.h"

class ProcessSearchJob : public Job {
    public:
        ProcessSearchJob(Parameters *params);
        ~ProcessSearchJob(void);
        std::vector<std::string> run(void);

    private:
        int  doProcessSearch(void);
        VectorXi movedAtoms(void);
        void printEndState(int status);
        void saveData(int status);

        std::vector<std::string> returnFiles;

        Parameters *parameters;
        SaddleSearchMethod *saddleSearch; 
        Matter *initial;      // initial configuration.
        Matter *saddle;       // configuration used during the saddle point search.
        Matter *displacement; // configuration used during the saddle point search.
        Matter *min1;         // first minimum from the saddle
        Matter *min2;         // second minimum from the saddle

        double barriersValues[2];
        double prefactorsValues[2];

        int fCallsSaddle;
        long fCallsMin;
        int fCallsPrefactors;
};

#endif

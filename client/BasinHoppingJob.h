//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef BASINHOPPINGJOB_H
#define BASINHOPPINGJOB_H

#include "Matter.h"
#include "Parameters.h"
#include "Job.h"

class BasinHoppingJob : public Job {
    public:
        BasinHoppingJob(Parameters *params);
        ~BasinHoppingJob(void);
        std::vector<std::string> run(void);

    private:
        VectorXd calculateDistanceFromCenter(Matter *matter);
        AtomMatrix displaceRandom(double maxDisplacement);
        void randomSwap(Matter *matter);
        Parameters *parameters;
        Matter *current;
        Matter *trial; // initial configuration
        vector<long> getElements(Matter *matter);
        std::vector<std::string> returnFiles;
        int jump_count; // count of jump moves
        int disp_count; // count of displacement moves
        int swap_count; // count of swap moves
        int fcalls;

        std::vector<Matter *> uniqueStructures;
        std::vector<double>   uniqueEnergies;
};

#endif

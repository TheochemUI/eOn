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
        AtomMatrix displaceRandom();
        AtomMatrix displaceSingle();

        Parameters *parameters;
        Matter *current;
        Matter *trial;  // initial configuration.

        std::vector<std::string> returnFiles;
};

#endif

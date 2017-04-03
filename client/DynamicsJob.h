#ifndef DYNAMICSJOB_H
#define DYNAMICSJOB_H

#include "Matter.h"
#include "Parameters.h"
#include "Job.h"

class DynamicsJob : public Job 
{

    public:

        DynamicsJob(Parameters *params);
        ~DynamicsJob(void);
        std::vector<std::string> run(void);
        Parameters *parameters;
        std::vector<std::string> returnFiles;
};

#endif

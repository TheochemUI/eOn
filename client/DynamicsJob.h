#ifndef DYNAMICSJOB_H
#define DYNAMICSJOB_H

#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

class DynamicsJob : public Job {
public:
    DynamicsJob(Parameters *params);
    ~DynamicsJob(void);
    Parameters *parameters;
    std::vector<std::string> run(void);
    std::vector<std::string> returnFiles;
};

#endif

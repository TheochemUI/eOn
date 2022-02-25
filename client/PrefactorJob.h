
#ifndef PREFACTORJOB_H
#define PREFACTORJOB_H

#include "Job.h"
#include "Parameters.h"

class PrefactorJob : public Job
{
public:
    PrefactorJob(Parameters *params);
    ~PrefactorJob();
    std::vector<std::string> run(void);
private:
    Parameters *parameters;
};

#endif

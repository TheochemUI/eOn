
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
    static const char PREFACTOR_REACTANT[];
    static const char PREFACTOR_SADDLE[];
    static const char PREFACTOR_PRODUCT[];
private:
    Parameters *parameters;
};

#endif

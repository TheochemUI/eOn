
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
    static const std::string PREFACTOR_REACTANT;
    static const std::string PREFACTOR_SADDLE;
    static const std::string PREFACTOR_PRODUCT;
private:
    Parameters *parameters;
};

#endif

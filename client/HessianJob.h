#ifndef HESSIANJOB_H
#define HESSIANJOB_H

#include "Job.h"
#include "Parameters.h"

class HessianJob : public Job
{
public:
    HessianJob(Parameters *params);
    ~HessianJob();
    std::vector<std::string> run(void);
private:
    Parameters *parameters;
};

#endif

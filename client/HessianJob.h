#ifndef HESSIANJOB_H
#define HESSIANJOB_H

#include "Job.h"
#include "Parameters.h"

class HessianJob : public Job
{
public:
    HessianJob(Parameters *params);
    ~HessianJob();
    void run(int bundleNumber);
private:
    Parameters *parameters;
};


#endif

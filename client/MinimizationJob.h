#ifndef MINIMAZATIONJOB_H
#define MINIMAZATIONJOB_H

#include "Job.h"
#include "Parameters.h"

class MinimizationJob: public Job {
    public:
        MinimizationJob(Parameters *params);
        ~MinimizationJob(void);
        std::vector<std::string> run(void);
    private:
        Parameters *parameters;
        int fcalls;
        enum {
            STATUS_GOOD,
            STATUS_MAX_ITERATIONS,
            STATUS_POTENTIAL_FAILED
        };
};

#endif

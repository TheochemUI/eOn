#ifndef MONTECARLOJOB_H
#define MONTECARLOJOB_H

#include "Job.h"
#include "Parameters.h"

class MonteCarloJob: public Job {
    public:
        MonteCarloJob(Parameters *params);
        ~MonteCarloJob(void);
        std::vector<std::string> run(void);
    private:
        Parameters *parameters;
};

#endif

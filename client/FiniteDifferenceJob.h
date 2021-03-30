#ifndef FINITEDIFFERENCE_H
#define FINITEDIFFERENCE_H


#include "Job.h"
#include "Parameters.h"
#include "Eigen.h"

class FiniteDifferenceJob: public Job {
    public:
        FiniteDifferenceJob(Parameters *params);
        ~FiniteDifferenceJob(void);
        std::vector<std::string> run(void);
    private:
        Parameters *parameters;
};

#endif

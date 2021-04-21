#ifndef POINTJOB_H
#define POINTJOB_H

#include "Job.h"
#include "Parameters.h"

class PointJob: public Job {
    public:
        PointJob(Parameters *params);
        ~PointJob(void);
        std::vector<std::string> run(void);
    private:
        Parameters *parameters;
};

#endif

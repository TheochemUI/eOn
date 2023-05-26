#ifndef POINTJOB_H
#define POINTJOB_H

#include "Job.h"
#include "Parameters.h"

class PointJob: public Job {
    public:
        PointJob(std::unique_ptr<Parameters> parameters)
            : Job(std::move(parameters)) {}
        ~PointJob(void) = default;
        std::vector<std::string> run(void);
};

#endif

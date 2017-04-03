
#ifndef STRUCTURECOMPARISONJOB_H
#define STRUCTURECOMPARISONJOB_H

#include "Job.h"
#include "Parameters.h"

class StructureComparisonJob: public Job {
    public:
        StructureComparisonJob(Parameters *params);
        ~StructureComparisonJob(void);
        std::vector<std::string> run(void);
    private:
        Parameters *parameters;
};

#endif

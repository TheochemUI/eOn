#include "Job.h"
#include "Parameters.h"

class ParallelReplica: public Job {
    public:
        ParallelReplica(Parameters *params);
        ~ParallelReplica(void);
        void run(int bundleNumber);
    private:
        Parameters *parameters;
};

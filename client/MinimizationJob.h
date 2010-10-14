#include "Job.h"
#include "Parameters.h"

class MinimizationJob: public Job {
    public:
        MinimizationJob(Parameters *params);
        ~MinimizationJob(void);
        void run(int bundleNumber);
    private:
        Parameters *parameters;
};

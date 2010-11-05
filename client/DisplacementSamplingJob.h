#include "Job.h"
#include "Parameters.h"
#include "Eigen/Eigen"

USING_PART_OF_NAMESPACE_EIGEN

class DisplacementSamplingJob: public Job {
    public:
        DisplacementSamplingJob(Parameters *params);
        ~DisplacementSamplingJob(void);
        void run(int bundleNumber);
    private:
        Parameters *parameters;
};

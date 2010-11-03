#include "Job.h"
#include "Parameters.h"
#include "Eigen/Eigen"

USING_PART_OF_NAMESPACE_EIGEN

class DimerDisplaceJob: public Job {
    public:
        DimerDisplaceJob(Parameters *params);
        ~DimerDisplaceJob(void);
        void run(int bundleNumber);
    private:
        Parameters *parameters;
};

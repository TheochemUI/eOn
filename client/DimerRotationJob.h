#include "Job.h"
#include "Parameters.h"
#include "Eigen/Eigen"

USING_PART_OF_NAMESPACE_EIGEN

class DimerRotationJob: public Job {
    public:
        DimerRotationJob(Parameters *params);
        ~DimerRotationJob(void);
        void run(int bundleNumber);
    private:
        Parameters *parameters;
};

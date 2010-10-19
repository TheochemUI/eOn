#include "Job.h"
#include "Parameters.h"
#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

class DimerDrJob: public Job {
    public:
        DimerDrJob(Parameters *params);
        ~DimerDrJob(void);
        void run(int bundleNumber);
    private:
        Parameters *parameters;
};

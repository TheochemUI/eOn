#include "Job.h"
#include "Parameters.h"

class TestJob: public Job {
    public:
        TestJob(Parameters *params);
        ~TestJob(void);
        void run(int bundleNumber);
    private:
        Parameters *parameters;
        double tolerance;
        void checkFullSearch(void);
        void checkMode(void);        
        void checkPotentials(void);
        double getEnergyDiff(long potTag, double refEnergy);
        double getForceDiff(long potTag, double refForce);
};

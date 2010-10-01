#include "Job.h"
#include "Parameters.h"

class Minimization: public Job {
    public:
        Minimization(Parameters *params);
        ~Minimization(void);
        void run(int bundleNumber);
    private:
        Parameters *parameters;
};

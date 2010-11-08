#include "Job.h"
#include "Parameters.h"

class ReplicaExchangeJob: public Job {
    public:
        ReplicaExchangeJob(Parameters *params);
        ~ReplicaExchangeJob(void);
        void run(int bundleNumber);
    private:
    Parameters *parameters;
	Matter *reactant;
    Matter *ReplicaArray[];
    double *ReplicaTemp;
};

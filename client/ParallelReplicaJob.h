#include "Job.h"
#include "Parameters.h"

class ParallelReplicaJob: public Job {
    public:
        ParallelReplicaJob(Parameters *params);
        ~ParallelReplicaJob(void);
        void run(int bundleNumber);
    private:
	void dynamics();
	bool CheckState(Matter *matter);
	void Refine(Matter *mdbuff[]);
    void saveData(int status,int bundleNumber);
    Parameters *parameters;
	Matter *reactant;
   	Matter *min1;
	Matter *min2;
    long *stepsbuff;
    double SPtime;
    double RLtime;
    double *SPtimebuff;
	long min_fcalls;
	long md_fcalls;
	long nsteps;
    long nsteps_refined;
	long check_steps;
	long relax_steps;
	bool newstate;
};

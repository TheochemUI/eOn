#include "Job.h"
#include "Parameters.h"

class ParallelReplica: public Job {
    public:
        ParallelReplica(Parameters *params);
        ~ParallelReplica(void);
        void run(int bundleNumber);
    private:
		void dynamics();
		bool firstArchieve(Matter *matter);
		bool IsNewState();
	//	long Refine(Matter &matter[]);
        void saveData(int status,int bundleNumber);
        Parameters *parameters;
		Matter *reactant;
   		Matter *min1;
		Matter *min2;
		long min_fcalls;
		long md_fcalls;
		long nsteps;
		long ncheck;
		long nexam;
		bool newstate;
		bool stoped;
		bool remember;
};

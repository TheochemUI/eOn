#include "Job.h"
#include "Parameters.h"

class ParallelReplica: public Job {
    public:
        ParallelReplica(Parameters *params);
        ~ParallelReplica(void);
        void run(int bundleNumber);
    private:
		void dynamics();
		bool firstAchieve();
		bool IsNewState();
        void saveData(int status,int bundleNumber);
        Parameters *parameters;
		Matter *reactant;
   		Matter *min1;
		Matter *min2;
		long min_fcalls;
		long md_fcalls;
		long nsteps;
		bool newstate;
		bool stoped;
};

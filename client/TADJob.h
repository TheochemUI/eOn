
#ifndef TADJOB_H
#define TADJOB_H

#include "Job.h"
#include "Parameters.h"
#include "MinModeSaddleSearch.h"

class TADJob: public Job
{
    public:

        TADJob(Parameters *params);
        ~TADJob(void);
        std::vector<std::string> run(void);

    private:

        int dynamics();
        long refine(Matter *mdBuffer[], long length, Matter *reactant);
        bool checkState(Matter *current, Matter *reactant);
        void saveData(int status);
        void dephase();
        bool saddleSearch(Matter *tran);
        Parameters *parameters;

        Matter *current;
        Matter *reactant;
        Matter *saddle;
        Matter *crossing;
        Matter *final;
        Matter *final_tmp;
        Matter *product;

        bool metaStateFlag;
        bool newStateFlag;

        long minimizeFCalls;
        long mdFCalls;
        long dephaseFCalls;
        long refineFCalls;

        long transitionStep;

        double time;
        double minCorrectedTime;
        double transitionTime;
        double barrier;
        double *timeBuffer;
        MinModeSaddleSearch *dimerSearch;

        std::vector<std::string> returnFiles;
};

#endif

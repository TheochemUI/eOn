
#ifndef SAFEHYPERJOB_H
#define SAFEHYPERJOB_H

#include "Job.h"
#include "Parameters.h"

class SafeHyperJob: public Job
{
    public:

        SafeHyperJob(Parameters *params);
        ~SafeHyperJob(void);
        std::vector<std::string> run(void);

    private:

        int dynamics();
        long refine(Matter *mdBuffer[], long length, Matter *reactant);
        bool checkState(Matter *current, Matter *reactant);
        void saveData(int status);
        void dephase();

        Parameters *parameters;

        Matter *current;
        Matter *reactant;
        Matter *saddle;
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
        double transitionPot;
        double *timeBuffer;
        double *biasBuffer;

        std::vector<std::string> returnFiles;
};

#endif

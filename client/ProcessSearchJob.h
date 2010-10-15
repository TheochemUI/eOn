/*
 *===============================================
 *  EON ProcessSearchJob.h
 *===============================================
 */
#ifndef PROCESSSEARCHJOB_H
#define PROCESSSEARCHJOB_H

#include "Matter.h"
#include "Parameters.h"
#include "SaddlePoint.h"
#include "Hessian.h"
#include "Job.h"

class ProcessSearchJob : public Job {
    public:
        ProcessSearchJob(Parameters *params);
        ~ProcessSearchJob(void);
        void run(int bundleNumber);

    private:
        int  doProcessSearch(void);
        void printEndState(int status);
        void saveData(int status, int bundleNumber);

        Hessian *hessian;
        Parameters *parameters;
        SaddlePoint *saddlePoint; 
        Matter *initial;      // initial configuration.
        Matter *saddle;       // configuration used during the saddle point search.
        Matter *displacement; // configuration used during the saddle point search.
        Matter *min1;         // first minimum from the saddle
        Matter *min2;         // second minimum from the saddle

        double barriersValues[2];
        double prefactorsValues[2];

        int fCallsSaddle;
        int fCallsMin;
        int fCallsPrefactors;
};

#endif

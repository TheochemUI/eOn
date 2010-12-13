#ifndef BASINHOPPINGJOB_H
#define BASINHOPPINGJOB_H

#include "Matter.h"
#include "Parameters.h"
#include "Job.h"

class BasinHoppingJob : public Job {
    public:
        BasinHoppingJob(Parameters *params);
        ~BasinHoppingJob(void);
        void run(int bundleNumber);

    private:
        void saveData(int status, int bundleNumber);
        Matrix<double, Eigen::Dynamic, 3> displaceRandom();
        Matrix<double, Eigen::Dynamic, 3> displaceSingle();
        Parameters *parameters;
        Matter *current;
        Matter *trial;  // initial configuration.
};

#endif

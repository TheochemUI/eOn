#ifndef FIRE_H
#define FIRE_H

#include "Optimizer.h"
#include "Matter.h"
#include "Parameters.h"

class FIRE : public Optimizer
{

    public:
        double dt;

        FIRE(ObjectiveFunction *objf, Parameters *parameters);
        ~FIRE();

        int step(double maxMove);
        // Variant for early stopping / trust region
        int step(const double maxMove,
                 const std::vector<Matter> ppoints,
                 const double max_dist,
                 bool& isClose);
        int run(int maxIterations, double maxMove);

    private:
        ObjectiveFunction *objf;
        Parameters *parameters;
        VectorXd v;
        double alpha;
        double alpha_start;
        int N, N_min;
        double f_inc;
        double f_dec;
        double f_a;
        double dt_max;
        int iteration;
};

#endif

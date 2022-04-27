
#ifndef QUICKMIN_H
#define QUICKMIN_H

#include "Optimizer.h"
#include "Matter.h"
#include "Parameters.h"

class Quickmin : public Optimizer
{

    public:
        double dt;

        Quickmin(ObjectiveFunction *objf, Parameters *parameters);
        ~Quickmin();

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
        VectorXd velocity;
        int iteration;
};

#endif

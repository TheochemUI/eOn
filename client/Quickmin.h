
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
        int run(int maxIterations, double maxMove);

    private:
        ObjectiveFunction *objf;
        Parameters *parameters;
        VectorXd velocity;
        int iteration;
};

#endif

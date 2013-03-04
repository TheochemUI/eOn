//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

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

        bool step(double maxMove);
        bool run(int maxIterations, double maxMove);

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

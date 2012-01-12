//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

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

        bool step(double maxMove);
        bool run(int maxIterations, double maxMove);

    private:
        ObjectiveFunction *objf;
        Parameters *parameters;
        VectorXd velocity;
        int iteration;
};

#endif

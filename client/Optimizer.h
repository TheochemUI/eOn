//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "Parameters.h"
#include "ObjectiveFunction.h"
#include "Eigen.h"

class Optimizer
{
    
    public:
        virtual ~Optimizer(){};
        virtual bool step(double maxMove) = 0;
        virtual bool run(int maxIterations, double maxMove) = 0;
        virtual VectorXd getStep() = 0;
        static Optimizer *getOptimizer(ObjectiveFunction *objf, Parameters *parameters);
};

#endif

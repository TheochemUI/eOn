//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef CG_H
#define CG_H

#include "Eigen.h"
#include "Matter.h"
#include "Optimizer.h"
#include "Parameters.h"
#include "HelperFunctions.h"

class ConjugateGradients : public Optimizer
{
    public:

        ConjugateGradients(ObjectiveFunction *objf, Parameters *parameters);
        ~ConjugateGradients();

        bool step(double maxMove);
        bool run(int maxIterations, double maxMove);
        VectorXd getStep();

    private:

        ObjectiveFunction *objf;
        Parameters *parameters;

        VectorXd direction;
        VectorXd directionOld;
        VectorXd directionNorm;
        VectorXd force;
        VectorXd forceOld;
};

#endif

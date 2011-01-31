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

#include "Minimizer.h"
#include "Matter.h"
#include "Parameters.h"

class Quickmin : public Minimizer
{

    public:
        Matter *matter;
        Parameters *parameters;
        double dt;

        Quickmin(Matter *matter, Parameters *parameters);
        ~Quickmin();

        void oneStep();
        void fullRelax();
        bool isItConverged(double convergeCriterion);
        void setOutput(int level);
        
    private:
        int outputLevel;
        Matrix<double, Eigen::Dynamic, 3> velocity; 
};

#endif

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

class Quickmin : public Minimizer
{

    public:
        Quickmin(Matter *matter, Parameters *parameters);

        ~Quickmin();

        void oneStep();
        void oneStepPart1(Matrix<double, Eigen::Dynamic, 3> forces);
        void oneStepPart2(Matrix<double, Eigen::Dynamic, 3> forces);
        
        void fullRelax();
        bool isItConverged(double convergeCriterion);
        void setOutput(int level);
        
    private:
        long nAtoms;
        int outputLevel;
        Matter *matter;
        Parameters *parameters;

        Matrix<double, Eigen::Dynamic, 3> forces;
        double dtScale;

        Matrix<double, Eigen::Dynamic, 3> velocity; 
};

#endif


//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef QUICKMINBOX_H
#define QUICKMINBOX_H

#include "Minimizer.h"
#include "Matter.h"
#include "Parameters.h"

class QuickminBox : public Minimizer
{

    public:
        Matter *matter;
        Parameters *parameters;
        double dt;
        Matrix<double, Eigen::Dynamic, 3> velocity; 
        double boxXv, boxYv, boxZv, boxXf, boxYf, boxZf;

        QuickminBox(Matter *matter, Parameters *parameters);
        ~QuickminBox();

        void oneStep();
        void fullRelax();
        bool isItConverged(double convergeCriterion);
        void setOutput(int level);
        
    private:
        int outputLevel;
};

#endif

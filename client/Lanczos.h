//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef LANCZOS_H
#define LANCZOS_H

#include "Eigen.h"
#include "Matter.h"
#include "Parameters.h"
#include "LowestEigenmode.h"

// Lanczos method to find the lowest curvature mode
class Lanczos : public LowestEigenmode
{

    public:
        Lanczos(Matter *matter, Parameters *parameters);
        ~Lanczos();

        void compute(Matter *matter, AtomMatrix direction);
        double getEigenvalue();
        AtomMatrix getEigenvector();

    private:
        Parameters *parameters;
        AtomMatrix lowestEv;
        double lowestEw;
};

#endif



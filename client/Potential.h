//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Parameters.h"
#include "Eigen.h"

class Potential
{

    public:

        virtual ~Potential();

        static Potential* getPotential(Parameters *parameters);

        static int fcalls;

        AtomMatrix force(long nAtoms, AtomMatrix positions,
                         VectorXi atomicNrs, double *energy, Matrix3d box);

        void virtual initialize() = 0;
        void virtual cleanMemory(){};
        void virtual force(long nAtoms, const double *positions, 
                           const int *atomicNrs, double *forces, double *energy, 
                           const double *box) = 0;

};

#endif

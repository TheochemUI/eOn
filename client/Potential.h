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

        Potential(){}
        virtual ~Potential(){}

/*
        static const string POT_LJ;
        static const string POT_EAM_AL;
        static const string POT_MORSE_PT;
        static const string POT_EMT;
        static const string POT_QSC;
        static const string POT_ZPICE;
        static const string POT_TIP4P;
        static const string POT_LENOSKY_SI;
        static const string POT_SW_SI;
        static const string POT_TERSOFF_SI;
        static const string POT_EDIP;
        static const string POT_VASP;
        static const string POT_BOPFOX;
        static const string POT_BOP;
        static const string POT_LAMMPS;
        static const string POT_GPAW;
*/
        static const char POT_LJ[];
        static const char POT_EAM_AL[];
        static const char POT_MORSE_PT[];
        static const char POT_EMT[];
        static const char POT_QSC[];
        static const char POT_ZPICE[];
        static const char POT_TIP4P[];
        static const char POT_LENOSKY_SI[];
        static const char POT_SW_SI[];
        static const char POT_TERSOFF_SI[];
        static const char POT_EDIP[];
        static const char POT_VASP[];
        static const char POT_BOPFOX[];
        static const char POT_BOP[];
        static const char POT_LAMMPS[];
        static const char POT_GPAW[];

        static Potential* getPotential(Parameters *parameters);

        static int fcalls;

        AtomMatrix force(long nAtoms, AtomMatrix positions,
                         VectorXi atomicNrs, double *energy, Matrix3d box);

        void virtual initialize() = 0;
        void virtual force(long nAtoms, const double *positions, 
                           const int *atomicNrs, double *forces, double *energy, 
                           const double *box) = 0;

};

#endif

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

        static const char POT_LJ[];
        static const char POT_LJCLUSTER[];
        static const char POT_IMD[];
        static const char POT_EAM_AL[];
        static const char POT_MORSE_PT[];
        static const char POT_EMT[];
        static const char POT_QSC[];
        static const char POT_TIP4P[];
        static const char POT_TIP4P_PT[];
        static const char POT_TIP4P_H[];
        static const char POT_SPCE[];
        static const char POT_LENOSKY_SI[];
        static const char POT_SW_SI[];
        static const char POT_TERSOFF_SI[];
        static const char POT_EDIP[];
        static const char POT_FEHE[];
        static const char POT_VASP[];
        static const char POT_BOPFOX[];
        static const char POT_BOP[];
        static const char POT_LAMMPS[];
        static const char POT_MPI[];
        static const char POT_TERMINAL[];
        static const char POT_NEW[];

        static Potential* getPotential(Parameters *parameters);

        static int fcalls;
        static int fcallsTotal;
        static int wu_fcallsTotal;
        static double totalUserTime;
        
        Parameters *params;

        AtomMatrix force(long nAtoms, AtomMatrix positions,
                         VectorXi atomicNrs, double *energy, Matrix3d box);

        void virtual initialize() = 0;
        void virtual force(long nAtoms, const double *positions,
                           const int *atomicNrs, double *forces, double *energy,
                           const double *box) = 0;

        static Potential* pot;

};

#endif

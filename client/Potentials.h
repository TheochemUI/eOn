//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef POTENTIALS_H
#define POTENTIALS_H

#include "Parameters.h"
#include "PotentialsInterface.h"
#include "Eigen.h"

/* Class serving as a wrapper between the force calculator and the Matter object */
class Potential{

public:

    enum{
        POT_LJ,
        POT_MORSE_PT,
        POT_EMT,
        POT_EDIP,
        POT_VASP,
        POT_TERSOFF_SI,
        POT_SW_SI,
        POT_LENOSKY_SI,
        POT_LJBINARY,
        POT_EAM_AL,
        POT_QSC,
        POT_ZPICE,
        POT_TIP4P,
        POT_BOPFOX,
        POT_BOP,
        POT_LAMMPS,
        N_POTS
    };

    Potential(Parameters *parameters);
 
    ~Potential();
 
    /* this function must be provided by each force calculator
    @param[in]    nAtoms      the number of atoms
    @param[in]    *positions  pointer to the array of 3N atoms positions
    @param[in]    *atomicNrs  pointer to the array of N atomic numbers
    @param[out]   *forces     pointer to the array of 3N forces
    @param[out]   *energy     pointer to the total energy
    @param[in]    *box        pointer to the array containing the 3 lengths of the supercell */
    Matrix<double, Eigen::Dynamic, 3> force(long nAtoms, Matrix<double, Eigen::Dynamic, 3> positions, Matrix<int, Eigen::Dynamic, 1> atomicNrs, double *energy, Matrix<double, 3, 3> box);
 
    static int fcalls;

private:

    PotentialsInterface *interface_;
    Parameters *parameters_;

};
#endif

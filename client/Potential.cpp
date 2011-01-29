//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "Constants.h"
#include "Potential.h"
#include "potentials/NewPotential/NewPotential.h"
#include "potentials/EDIP/EDIP.h"
//#include "potentials/EMT/EffectiveMediumTheory.h"
#include "potentials/Morse/Morse.h"
#include "potentials/LennardJones/LJ.h"
#include "potentials/SW/SW.h"
#include "potentials/Tersoff/Tersoff.h"
#include "potentials/Aluminum/Aluminum.h"
#include "potentials/EAM/EAM.h"
#include "potentials/Lenosky/Lenosky.h"
#include "potentials/QSC/QSC.h"
#include "potentials/platinum-water/zhu_philpott_for_eon.hpp"
#ifndef WIN32
    #include "potentials/VASP/VASP.h"
#endif
#include "potentials/bopfox/bopfox.h"
#ifdef BOPFOX
    #include "potentials/bop/bop.h"
#endif
#ifdef LAMMPS_POT
    #include "potentials/LAMMPS/LAMMPS_EON.h"
#endif

#include <cstdlib>

Potential* Potential::getPotential(Parameters *parameters)
{
    Potential* pot;
    if(parameters->potential == "lj")
        pot = new LJ();
    else if(parameters->potential == "morse_pt")
        pot = new Morse();
//    else if(parameters->potential == "emt")
//        pot = new EffectiveMediumTheory();
    else if(parameters->potential == "qsc")
        pot = new QSC();
    else if(parameters->potential == "zpice")
        pot = new ZpIce();
    else if(parameters->potential == "tip4p")
        pot = new Tip4p();
    #ifndef NO_FORTRAN
    else if(parameters->potential == "eam_al")
        pot = new Aluminum();
    else if(parameters->potential == "lenosky_si")
        pot = new Lenosky();
    else if(parameters->potential == "sw_si")
        pot = new SW();
    else if(parameters->potential == "tersoff_si")
        pot = new Tersoff();
    else if(parameters->potential == "edip")
        pot = new EDIP();
    #ifndef WIN32
    else if(parameters->potential == "vasp")
        pot = new VASP();
    #endif
    else if(parameters->potential == "bopfox")
        pot = new bopfox();
    #endif
    #ifdef BOPFOX
    else if(parameters->potential == "bop")
        pot = new bop();
    #endif
    #ifdef LAMMPS_POT
    else if(parameters->potential == "lammps")
        pot = new lammps_eon();
    #endif
    else
    {
        printf("Potential tag not recognized: %ld\n", parameters->potential);
        std::exit(1);
    }	
    pot->initialize();
    return pot;
};

int Potential::fcalls = 0;

Potential::~Potential()
{
    //cleanMemory();
};

AtomMatrix Potential::force(long nAtoms, AtomMatrix positions, Matrix<int, Eigen::Dynamic, 1> atomicNrs, double *energy, Matrix<double, 3, 3> box) 
{
    // Eigen stores data in column-major format but we want row-major
    AtomMatrix forces(nAtoms,3);

    force(nAtoms, positions.data(), atomicNrs.data(), forces.data(), energy, box.data());

    fcalls+=1;

    return forces;
};

















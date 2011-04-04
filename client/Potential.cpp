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
#include "Parameters.h"
#include "potentials/NewPotential/NewPotential.h"
#include "potentials/EDIP/EDIP.h"
#include "potentials/EMT/EffectiveMediumTheory.h"
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
#ifdef MPIGPAW
    #include "potentials/GPAW/GPAW.h"
#endif

#include <cstdlib>

Potential *Potential::getPotential(Parameters *parameters)
{
    Potential* pot;
    if(parameters->potential == Parameters::LJ)
        pot = new LJ();
    else if(parameters->potential == Parameters::MORSE_PT)
        pot = new Morse();
    else if(parameters->potential == Parameters::EMT)
        pot = new EffectiveMediumTheory();
    else if(parameters->potential == Parameters::QSC)
        pot = new QSC();
    else if(parameters->potential == Parameters::ZPICE)
        pot = new ZpIce();
    else if(parameters->potential == Parameters::TIP4P)
        pot = new Tip4p();
    #ifndef NO_FORTRAN
    else if(parameters->potential == Parameters::EAM_AL)
        pot = new Aluminum();
    else if(parameters->potential == Parameters::LENOSKY_SI)
        pot = new Lenosky();
    else if(parameters->potential == Parameters::SW_SI)
        pot = new SW();
    else if(parameters->potential == Parameters::TERSOFF_SI)
        pot = new Tersoff();
    else if(parameters->potential == Parameters::EDIP)
        pot = new EDIP();
    #ifndef WIN32
    else if(parameters->potential == Parameters::VASP)
        pot = new VASP();
    #endif
    else if(parameters->potential == Parameters::BOPFOX)
        pot = new bopfox();
    #endif
    #ifdef BOPFOX
    else if(parameters->potential == Parameters::BOP)
        pot = new bop();
    #endif
    #ifdef LAMMPS_POT
    else if(parameters->potential == Parameters::LAMMPS)
        pot = new lammps_eon();
    #endif
    #ifdef MPIGPAW
    else if(parameters->potential == Parameters::GPAW)
        pot = new GPAW();
    #endif
    else
    {
        printf("Unknown Potential: %s\n", parameters->potential.c_str());
        std::exit(1);
    }
    pot->initialize();
    return pot;
};

int Potential::fcalls = 0;

AtomMatrix Potential::force(long nAtoms, AtomMatrix positions, VectorXi atomicNrs, double *energy, Matrix3d box) 
{
    AtomMatrix forces(nAtoms,3);

    force(nAtoms, positions.data(), atomicNrs.data(), forces.data(), energy, box.data());

    fcalls+=1;

    return forces;
};

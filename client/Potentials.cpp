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
#include "Potentials.h"
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

#include <cstdlib>

//int Potential::fcalls = 0;
int Potential::fcalls = 0;
// which potential to use is decided at preprocessor level
Potential::Potential(Parameters *parameters){
    parameters_ = parameters;

    if(parameters_->potential == "lj"){
        interface_ = new LJ();
        interface_->initialize();
    }
    else if(parameters_->potential == "morse_pt"){
        interface_ = new Morse();
        interface_->initialize();
    }
    else if(parameters_->potential == "emt"){
        interface_ = new EffectiveMediumTheory();
        interface_->initialize();
    }
    else if(parameters_->potential == "qsc"){
        interface_ = new QSC();
        interface_->initialize();
    }
    else if(parameters_->potential == "zpice"){
        interface_ = new ZpIce();
        interface_->initialize();
    }
    else if(parameters_->potential == "tip4p"){
        interface_ = new Tip4p();
        interface_->initialize();
    }
#ifndef NO_FORTRAN
    else if(parameters_->potential == "eam_al"){
        interface_ = new Aluminum();
        interface_->initialize();
    }
    else if(parameters_->potential == "lenosky_si"){
        interface_ = new Lenosky();
        interface_->initialize();
    }
    else if(parameters_->potential == "sw_si"){
        interface_ = new SW();
        interface_->initialize();
    }
    else if(parameters_->potential == "tersoff_si"){
        interface_ = new Tersoff();
        interface_->initialize();
    }
    else if(parameters_->potential == "edip"){
        interface_ = new EDIP();
        interface_->initialize();
    }
#ifndef WIN32
    else if(parameters_->potential == "vasp"){
        interface_ = new VASP();
        interface_->initialize();
    }		
#endif
    else if(parameters_->potential == "bopfox"){
        interface_ = new bopfox();
        interface_->initialize();
    }
#endif
#ifdef BOPFOX
    else if(parameters_->potential == "bop"){
        interface_ = new bop();
        interface_->initialize();
    }
#endif
#ifdef LAMMPS_POT
    else if(parameters_->potential == "lammps"){
        interface_ = new lammps_eon();
        interface_->initialize();
    }
#endif
    else{
        printf("Potential tag not recognized: %ld\n", parameters_->potential);
        std::exit(1);
    }	
};

Potential::~Potential(){
    interface_->cleanMemory();
};

// An alike function should be provided by the force calculator.
AtomMatrix Potential::force(long nAtoms, AtomMatrix positions, Matrix<int, Eigen::Dynamic, 1> atomicNrs, double *energy, Matrix<double, 3, 3> box) 
{
    //XXX: For now, this just serves as a wrapper for the potentials
    //     and converts from Matrix to double[]s. Later, the potentials
    //     should also use Eigen

    // Eigen stores data in column-major format but we want row-major
    AtomMatrix forces(nAtoms,3);

    interface_->force(nAtoms, positions.data(), atomicNrs.data(), forces.data(), energy, box.data());

    fcalls+=1;
    return forces;
};

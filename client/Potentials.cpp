/*
 *===============================================
 *  EON Potentials.cpp
 *===============================================
 */

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
#include "potentials/VASP/VASP.h"
#include "potentials/bopfox/bopfox.h"

#include <cstdlib>

// which potential to use is decided at preprocessor level
Potentials::Potentials(Parameters *parameters){
    parameters_ = parameters;
//_______________________    
// To use a new potential.
// An interface should be created in the file NewPotential_interface.cpp. 
// Code will stop at runtime if used and no new potential has been defined!
    if(parameters_->potentialTag == POT_USER){
        //interface_ = new NewPotential();
        //interface_->initialize();
        printf("The new potential must be commented in Potentials.cpp.\n");
        std::exit(1);
    }
//_______________________  
    else if(parameters_->potentialTag == POT_LJ){
        interface_ = new LJ();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_MORSE){
        interface_ = new Morse();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_EMT){
        interface_ = new EffectiveMediumTheory();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_EAM){
        interface_ = new EAM();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_QSC){
        interface_ = new QSC();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_ZPICE){
        interface_ = new ZpIce();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_TIP4P){
        interface_ = new Tip4p();
        interface_->initialize();
    }
#ifndef NO_FORTRAN
    else if(parameters_->potentialTag == POT_ALUMINUM){
        interface_ = new Aluminum();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_LENOSKY){
        interface_ = new Lenosky();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_SW){
        interface_ = new SW();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_TERSOFF){
        interface_ = new Tersoff();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_EDIP){
        interface_ = new EDIP();
        interface_->initialize();
    }
    else if(parameters_->potentialTag == POT_VASP){
        interface_ = new VASP();
        interface_->initialize();
    }		
    else if(parameters_->potentialTag == POT_BOPFOX){
        interface_ = new bopfox();
        interface_->initialize();
    }		
#endif
    else{
        printf("Potential tag not recognized: %ld\n", parameters_->potentialTag);
        std::exit(1);
    }	
};

Potentials::~Potentials(){
    interface_->cleanMemory();
};

// An alike function should be provided by the force calculator.
void Potentials::force(long nAtoms, const double *positions, const long *atomicNrs, 
                          double *forces, double *energy, const double *box){
    // The call is passed to the specified force calculator.
    interface_->force(nAtoms, positions, atomicNrs, forces, energy, box);
    
    if(parameters_->potentialNoTranslation){
        double tempForceX = 0;
        double tempForceY = 0;
        double tempForceZ = 0;
        
        for(long int i=0; i<nAtoms; i++) {
            tempForceX = tempForceX+forces[ 3*i ];
            tempForceY = tempForceY+forces[3*i+1];
            tempForceZ = tempForceZ+forces[3*i+2];
        }
        tempForceX = tempForceX/nAtoms;
        tempForceY = tempForceY/nAtoms;
        tempForceZ = tempForceZ/nAtoms;
        
        for(long int i=0; i<nAtoms; i++) {
            forces[ 3*i ] = forces[ 3*i ]-tempForceX;
            forces[3*i+1] = forces[3*i+1]-tempForceY;
            forces[3*i+2] = forces[3*i+2]-tempForceZ;
        }
    }
    return;
};

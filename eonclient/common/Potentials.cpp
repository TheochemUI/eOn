/*
 *===============================================
 *  Created by Andreas Pedersen on 10/4/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#include "Potentials.h"

using namespace constants;

// which potential to use is decided at preprocessor level
Potentials::Potentials(Parameters *parameters){
    parameters_ = parameters;
// _______________________    
// To use a new potential.
// An interface should be created in the file NewPotential_interface.cpp. Code will crash at runtime if used and no new potential has been defined!
    if(parameters_->getPotentialTag() == getPotentialNewPotential()){
        //interface_ = new NewPotential();
        //interface_->initialize();
    }
// _______________________  
#if 0
    else if(parameters_->getPotentialTag() == getPotentialLJ()){
        interface_ = new LJ();
        interface_->initialize();
    }
#endif
#ifdef MORSE
    else if(parameters_->getPotentialTag() == getPotentialMorse()){
        interface_ = new Morse();
        interface_->initialize();
    }
#endif
#ifdef EFFECTIVE_MEDIUM_THEORY
    else if(parameters_->getPotentialTag() == getPotentialEMT()){
        interface_ = new EffectiveMediumTheory();
        interface_->initialize();
    }
#endif
#if 0
    else if(parameters_->getPotentialTag() == getPotentialEDIP()){
        interface_ = new EDIP();
        interface_->initialize();
    }
    else if(parameters_->getpotentialtag() == getpotentialvasp()){
        interface_ = new vasp();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialTersoff()){
        interface_ = new Tersoff();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialSW()){
        interface_ = new SW();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialLenosky()){
        interface_ = new Lenosky();
        interface_->initialize();
    }
    else if(parameters_->getPotentialTag() == getPotentialLJBinary()){
        interface_ = new LJBinary();
        interface_->initialize();
    }
#endif
#ifdef ALUMINUM_POTENTIAL
    else if(parameters_->getPotentialTag() == getPotentialAluminum()){
        interface_ = new Aluminum();
        interface_->initialize();
    }
#endif
};

Potentials::~Potentials(){
    interface_->cleanMemory();
};

// An alike function should be provided by the force calculator.
void Potentials::force(long nAtoms, const double *positions, const long *atomicNrs, 
                          double *forces, double *energy, const double *box){
    // The call is passed to the specified force calculator.
    interface_->force(nAtoms, positions, atomicNrs, forces, energy, box);
    
    if(parameters_->getPotentialNoTranslation()){
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

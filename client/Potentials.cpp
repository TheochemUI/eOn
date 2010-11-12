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
#ifdef BOPFOX
    #include "potentials/bop/bop.h"
#endif

#include <cstdlib>

//int Potentials::fcalls = 0;
int Potentials::fcalls = 0;
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
    //else if(parameters_->potentialTag == POT_EAM){
    //    interface_ = new EAM();
    //    interface_->initialize();
    //}
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
#ifdef BOPFOX
    else if(parameters_->potentialTag == POT_BOP){
        interface_ = new bop();
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
Matrix<double, Eigen::Dynamic, 3> Potentials::force(long nAtoms, Matrix<double, Eigen::Dynamic, 3> positions, Matrix<int, Eigen::Dynamic, 1> atomicNrs, double *energy, Matrix<double, 3, 3> box) 
{
    //XXX: For now, this just serves as a wrapper for the potentials
    //     and converts from Matrix to double[]s. Later, the potentials
    //     should also use Eigen

    // Eigen stores data in column-major format but we want row-major
    Matrix<double, 3, 3> boxT= box.transpose();
    Matrix<double, 3, Eigen::Dynamic> positionsT = positions.transpose();
    Matrix<double, 3, Eigen::Dynamic> forcesT(3, (int)nAtoms);

    interface_->force(nAtoms, positionsT.data(), atomicNrs.data(), forcesT.data(), energy, boxT.data());

    Matrix<double, Eigen::Dynamic, 3> forces = forcesT.transpose();  

    if(parameters_->potentialNoTranslation){
        //XXX: What does this do?
        Vector3d tempForce(3);
        tempForce = forces.colwise().sum()/nAtoms;

        for(long int i=0; i<nAtoms; i++) 
        {
            forces.row(i) -= tempForce.transpose();
        }
    }

    fcalls+=1;
    return forces;
};


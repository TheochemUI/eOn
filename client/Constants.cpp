/*
 *===============================================
 *  Created by Andreas Pedersen on 10/27/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */

#include "Constants.h"

// Constants used to determine which potential to use.
// _______________________
// To use a new potential an interface should be created in the file NewPotential_interface.cpp.
// Code will crash at runtime if a potential is used and no new potential has been defined.
long constants::getPotentialNewPotential(){return 0;};
// _______________________    
long constants::getPotentialLJ(){return 1;};
long constants::getPotentialMorse(){return 2;};
long constants::getPotentialEMT(){return 3;};
long constants::getPotentialEDIP(){return 4;};
long constants::getPotentialVASP(){return 5;};
long constants::getPotentialTersoff(){return 6;};
long constants::getPotentialSW(){return 7;};
long constants::getPotentialLenosky(){return 8;};
//long constants::getPotentialLJBinary(){return 9;};
long constants::getPotentialAluminum(){return 10;};
long constants::getPotentialEAM(){return 11;};
long constants::getPotentialQSC(){return 12;};
long constants::getPotentialZpIce(){return 13;};
long constants::getPotentialTip4p(){return 14;};

// Values passed between server and clients
long constants::getStateGood(){return 0;};
long constants::getStateInit(){return 1;};
long constants::getStateSaddleSearchNoConvexRegion(){return 2;};
long constants::getStateSaddleSearchTerminatedBarrier(){return 3;};
long constants::getStateSaddleSearchTerminatedTotalIterations(){return 4;};
long constants::getStateSaddleSearchTerminatedConcaveIterations(){return 5;};
long constants::getStateNotConnected(){return 6;};
long constants::getStateBadPrefactor(){return 7;};
long constants::getStateBadBarrier(){return 8;};
long constants::getStateMinimumNotConverged(){return 9;};

// Constants used in the client
double constants::getMaxDifferencePos(){return 0.1;};

// Constants used in prefactor determination
double constants::getPrefactorMax(){return 10e20;};
double constants::getPrefactorMin(){return 10e8;};

// Constants used in conjugate gradients
double constants::getCurvatureStep(){return 0.001;};
double constants::getMaxMoveFullRelax(){return 0.20;};
double constants::getTimeStep(){return 0.1;};

// Constants used in displacement prior saddle point search
long constants::getNoPerturbation(){return 0;};
long constants::getPerturbateNotFccOrHcp(){return 1;};
long constants::getPerturbateMinimalCoordinated(){return 2;};
long constants::getPerturbateLastAtom(){return 3;};
double constants::getNeighborCutoff(){return 3.3;};

// Constants used in saddle point determination to set which algorithm that is going to be used for the lowest eigenmode determination
long constants::getLowestEigenmodeDimer(){return 1;};
long constants::getLowestEigenmodeLanczos(){return 2;};

// Constant used in saddle point search
//long constants::getMaximumJumpAttemptsSP(){return 50;};
long constants::getMaximumIterationsConcaveSeriesSP(){return 256;};

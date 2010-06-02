/*
 *===============================================
 *  DimerConstants.h
 *  eon2
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 04/02/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef DIMER_CONSTANTS_H
#define DIMER_CONSTANTS_H

#include "Constants.h"

/** Collection of dimer constants.*/

namespace constants {
    // Constants used in dimer  
    double getDimerSize();///< Distance between the configurations defining the dimer.
    double getDimerRotationAngle();///< The infinitesimal rotation step to determine the optimal rotational angle.
}
#endif

/*
 *===============================================
 *  EpiCenters.h
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 10/24/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef EPI_CENTERS_H
#define EPI_CENTERS_H

#include <vector>
#include <limits>
#include <cassert>
#include <climits>

#include "Matter.h"
#include "HelperFunctions.h"
#include "Constants.h"

namespace epi_centers {
    long cnaEpiCenter(const Matter *matter);///< The value returned is the index of a random atom being both free and neither FCC nor HCP coordinated.
    long minimalCoordinatedEpiCenter(const Matter *matter);///< The value returned is the index of a random atom being both free and minimally coordinated.
    long lastAtom(const Matter *matter);///< The value returned is the index of a last atom which is assumed to be free.
    long randomFreeAtomEpiCenter(const Matter *matter);///< The value returned is the index of a random atom being free.

    void cna(long *cna, const Matter *matter);///< Perform a CNA, the result is stored in \a cna that should be of the same length as the number of matter. 
    void coordination(long *coordinationVal, const Matter *matter);///< Determirne the coordination for the individual atoms, the result is stored in \a coordinationVal that should be of the same length as the number of matter. 
    void coordinationEqualOrBellow(bool *result, long coordinationMaxVal, const Matter *matter);///< Determine which atoms that is equal or less coordinationed than \a coordinationMaxVal, the result is stored in \a result that should be of the same length as the number of matter. 
    long minCoordination(const Matter *matter);///< Determine the \b first minimally coordinated atom. 
}
#endif

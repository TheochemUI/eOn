//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
//-----------------------------------------------------------------------------------
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

#include "Matter.h"

// radius used in the local atomic structure analysis: MOVE to paramters!
#define neighborCutoff  3.3

namespace EpiCenters {
    long cnaEpiCenter(const Matter *matter); // index of a random atom being both free and neither FCC nor HCP coordinated
    long minimalCoordinatedEpiCenter(const Matter *matter); // index of a random atom being both free and minimally coordinated
    long lastAtom(const Matter *matter); // index of a last atom which is assumed to be free
    long randomFreeAtomEpiCenter(const Matter *matter); // index of a random atom being free

    void cna(long *cna, const Matter *matter); // do a common neighbor analysis
    void coordination(long *coordinationVal, const Matter *matter); // determine the coordination for the individual atoms
    void coordinationEqualOrBellow(bool *result, long coordinationMaxVal, const Matter *matter); // determine which atoms that is equal or less coordinationed than coordinationMaxVal
    long minCoordination(const Matter *matter); // Determine the first minimally coordinated atom
}
#endif

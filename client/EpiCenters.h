//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef EPI_CENTERS_H
#define EPI_CENTERS_H

#include "Matter.h"

namespace EpiCenters
{
    long cnaEpiCenter(const Matter *matter, double neighborCutoff); // index of random atom that is free and neither FCC nor HCP coordinated
    long minCoordinatedEpiCenter(const Matter *matter, double neighborCutoff); // index of a random atom that is free and minimally coordinated
    long lastAtom(const Matter *matter); // index of last atom -- assumed to be free
    long randomFreeAtomEpiCenter(const Matter *matter); // index of a random atom that is free

    void cna(long *cna, const Matter *matter, double neighborCutoff); // do a common neighbor analysis
    void coordination(long *coordinationVal, const Matter *matter, double neighborCutoff); // determine the coordination for individual atoms
    void coordinationLessOrEqual(bool *result, long coordinationMaxVal, const Matter *matter, double neighborCutoff); // determine atoms with coordination <= coordinationMaxVal
    long minCoordination(const Matter *matter, double neighborCutoff); // determine first minimally coordinated atom
}
#endif

#ifndef EPI_CENTERS_H
#define EPI_CENTERS_H

#include "Matter.h"

namespace EpiCenters {

const char DISP_LOAD[] = "load";
const char DISP_NOT_FCC_OR_HCP[] = "not_fcc_hcp_coordinated";
const char DISP_MIN_COORDINATED[] = "least_coordinated";
const char DISP_LAST_ATOM[] = "last_atom";
const char DISP_RANDOM[] = "random";

// index of random atom that is free and neither FCC nor HCP coordinated
long cnaEpiCenter(const Matter *matter, double neighborCutoff);

// index of a random atom that is free and minimally coordinated
long minCoordinatedEpiCenter(const Matter *matter, double neighborCutoff);

// index of last atom -- assumed to be free
long lastAtom(const Matter *matter);

// index of a random atom that is free
long randomFreeAtomEpiCenter(const Matter *matter);

// do a common neighbor analysis
void cna(long *cna, const Matter *matter, double neighborCutoff);

// determine the coordination for individual atoms
void coordination(long *coordinationVal, const Matter *matter,
                  double neighborCutoff);

// determine atoms with coordination <= coordinationMaxVal
void coordinationLessOrEqual(bool *result, long coordinationMaxVal,
                             const Matter *matter, double neighborCutoff);

// determine first minimally coordinated atom
long minCoordination(const Matter *matter, double neighborCutoff);
} // namespace EpiCenters
#endif

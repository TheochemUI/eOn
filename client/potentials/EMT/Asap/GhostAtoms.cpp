#include "GhostAtoms.h"
#include "Vec.h"

GhostAtoms::GhostAtoms(const Vec *pos, int nAtoms, SuperCell *superCell)
    : Atoms(pos, nAtoms, superCell), ghostPotential(NULL), nGhosts(0) {}

void GhostAtoms::SetNumberOfGhosts(int g) {
  nGhosts = g;
  nImages = 0;
  positions.resize(nAtoms + nGhosts);
  types.resize(nAtoms + nGhosts);
}

#ifndef GHOSTATOMS_H
#define GHOSTATOMS_H

#include "Atoms.h"

class SuperCell;
class Vec;
class GhostPotential;

class GhostAtoms : public Atoms {
public:
  GhostAtoms(const Vec *pos, int nAtoms, SuperCell *superCell);
  virtual ~GhostAtoms() {}
  GhostPotential *GetGhostPotential() { return ghostPotential; }
  void SetGhostPotential(GhostPotential *g) { ghostPotential = g; }
  int GetNumberOfGhosts() const { return nGhosts; }
  void SetNumberOfGhosts(int g);

private:
  GhostPotential *ghostPotential;
  int nGhosts;
};

#endif // GHOSTATOMS_H

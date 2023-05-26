// Emacs: This is -*- C++ -*-

#include "Vec.h"
#include <vector>
using std::vector;

class Atoms;
class SuperCell;
struct Image;

/// CellList is a helper class for NeighborList and QCPotential.

/// CellList(atoms, rMax) will put the atoms into cells of size minimum
/// rMax (when used by the NeighborList class, rMax will be rCut + 2 * drift).
/// MakeNeighborLists(...) will generate neighbor lists and periodic images.
/// GetNeighbors(...) will return all neighbors of point within a distance
/// of rMax (useful for decoration in QC).

class CellList {
public:
  CellList(Atoms *a, double rMax);
  void MakeNeighborList(int neighborList[], int indices[], int &nMax,
                        vector<Image> &images);
  int GetNeighbors(Vec point, int *neighbors, int size) const;
  bool AnyNeighbors(Vec point) const;

private:
  int nAtoms;
  int nAllAtoms;
  double rMax;
  double rMax2;
  Vec minimum, size;
  int nCells[3];
  int nTotalCells[4];
  Vec *positions;
  const SuperCell *superCell;
  const double *superCellHeights;
  const bool *periodic;
  vector<vector<int>> cells;
  vector<int> cellIndices;
};

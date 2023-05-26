#include "NeighborList.h"
#include "Atoms.h"
#include "CellList.h"
#include "GhostAtoms.h"
#include "GhostPotential.h"
#include "SuperCell.h"
// #include "Timing.h"
#include <assert.h>
#include <math.h>
#include <string.h>

using std::cerr;
using std::endl;
using std::flush;
extern int verbose;

NeighborList::NeighborList(Atoms *a, double rCut, double driftfactor,
                           int defaultNeighborEstimate)
    : atoms(a) {
  // Driftfactor defaults to 0.05, but is set to 0.0 by CNA and friends
  if (verbose > 2)
    cerr << "NeighborList::NeighborList beginning" << endl;
  indices = 0;
  Allocate(); // Initializes nAtoms and nSize
  neighborList = 0;
  rCut2 = rCut * rCut;
  double drift = driftfactor * rCut;
  drift2 = drift * drift;
  rMax = rCut + 2.0 * drift;
  positions = atoms->GetPositionsPtr();
  // following is estimated number of neighbors. The lists will only store
  // half-neighbors, but a safety factor of 2 is a good thing. We make sure to
  // have at least 50 neighbors per atom though; for example a surface
  // involving a vacuum region would lower the apparent density
  double density = nSize / atoms->GetSuperCell()->GetVolume();
  double PI = 4.0 * atan(1.0);
  int estNbrsPerAtom = (int)(density * 4 * PI * rCut * rCut * rCut / 3.);
  if (estNbrsPerAtom < defaultNeighborEstimate)
    estNbrsPerAtom = defaultNeighborEstimate;

  nMax = estNbrsPerAtom * nSize;

  if (verbose > 2)
    cerr << "NeighborList::NeighborList calling MakeList()" << endl;
  atoms->NormalizePositions();
  MakeList();
  if (verbose > 2)
    cerr << "NeighborList::NeighborList calling UpdateImagePositions()" << endl;
  UpdateImagePositions();
  if (verbose > 2)
    cerr << "NeighborList::NeighborList completed." << endl;
  CheckAndUpdateNeighborList();
}

NeighborList::~NeighborList() {
  if (indices != 0)
    delete[] indices;
  if (neighborList != 0)
    delete[] neighborList;
}

void NeighborList::Allocate() {
  nAtoms = atoms->GetNumberOfRealAtoms();
  nSize = nAtoms;
  // Do we have any ghosts?
  GhostAtoms *ghostAtoms;
  //  try {
  //	ghostAtoms = dynamic_cast<GhostAtoms *>(atoms); //use below line for
  // windows

  //  }catch(...)
  //  {
  ghostAtoms = NULL;
  //  }

  if (ghostAtoms != 0)
    nSize += ghostAtoms->GetNumberOfGhosts(); // Yes!
  if (indices != 0)
    delete[] indices;
  indices = new int[nAtoms + 1];
}

/// \bug Changing the supercell while scaling the atomic positions has
/// little effect on the interatomic distances, but that little effect
/// should be included in the decision to update the neighborlist.
/// That is not done, so currently this function doesn't do anything.
void NeighborList::UpdateSuperCell(const SuperCell *superCell) {}

void NeighborList::MakeList() {
  // USETIMER("NeighborList::MakeList");
  if (verbose >= 1)
    cerr << " NeighborList-Update ";

  if (verbose > 2) {
    Vec realmax;
    Vec realmin;
    Vec ghmax;
    Vec ghmin;
    const Vec *pos = atoms->GetCartesianPositions();
    realmin = realmax = pos[0];
    for (int i = 1; i < nAtoms; i++)
      for (int j = 0; j < 3; j++) {
        if (realmax[j] < pos[i][j])
          realmax[j] = pos[i][j];
        if (realmin[j] > pos[i][j])
          realmin[j] = pos[i][j];
      }
    if (nSize > nAtoms) {
      ghmin = ghmax = pos[nAtoms];
      for (int i = nAtoms; i < nSize; i++)
        for (int j = 0; j < 3; j++) {
          if (ghmax[j] < pos[i][j])
            ghmax[j] = pos[i][j];
          if (ghmin[j] > pos[i][j])
            ghmin[j] = pos[i][j];
        }
    }
    cerr << "  Real atoms [" << nAtoms << "]:  " << realmin << " - " << realmax
         << endl;
    if (nSize > nAtoms)
      cerr << "  Ghost atoms [" << nSize - nAtoms << "]: " << ghmin << " - "
           << ghmax << endl;
  }

  // put atoms into cells:
  CellList cellList(atoms, rMax);

  // delete old images and old neighbor lists:
  images.resize(0);
  if (neighborList != 0)
    delete[] neighborList;

  // make space for new neighbor lists slightly bigger than the last one:
  //  nMax += int(0.15 * nMax) + 100;
  nMax += int(0.25 * nMax) + 100;
  neighborList = new int[nMax];

  // make the new neighbor list and new image atoms:
  cellList.MakeNeighborList(neighborList, indices, nMax, images);

  atoms->SetNumberOfImages(images.size());
  // this could have changed the adress of positions. Get it again:
  positions = atoms->GetPositionsPtr();

  // update reference positions for neighbor list:
  referencePositions.resize(nAtoms);
  memcpy(&referencePositions[0], &positions[0], nAtoms * 3 * sizeof(double));

  // Now find the maximal list length
  maxlistlen = 0;
  for (int i = 0; i < nAtoms; i++)
    if (indices[i + 1] - indices[i] > maxlistlen)
      maxlistlen = indices[i + 1] - indices[i];
  if (verbose > 1)
    cerr << "NeighborList::MakeList completed\n";
}

bool NeighborList::CheckAndUpdateNeighborList() {
  // USETIMER("NeighborList::CheckAndUpdateNeighborList");
  // if (verbose > 2)
  //   cerr << "NeighborList::CheckAndUpdate" << endl;

  bool updateRequired = false;
  for (int n = 0; n < nAtoms; n++)
    if (Length2(positions[n] - referencePositions[n]) > drift2) {
      updateRequired = true;
      break;
    }

  // try to get GhostAtoms class:
  GhostAtoms *ghostAtoms;
  //	try {
  //		ghostAtoms = dynamic_cast<GhostAtoms *>(atoms);
  //
  //  }catch(...)
  //  {
  ghostAtoms = NULL;
  //  }

  // if any atom on any processor has moved more than "drift" then migrate all
  // atoms and make new ghost atoms:
  bool reallocationRequired = false;
  if (ghostAtoms != 0) {
    updateRequired =
        ghostAtoms->GetGhostPotential()->CheckAndUpdateAtoms(updateRequired);
    // We might need to reallocate some stuff, if the number of atoms or
    // the number of ghosts has changed. Also the potential might want to
    // do some reallocation (that's why we return "reallocationRequired"
    // from this method.
    if (atoms->GetNumberOfRealAtoms() != nAtoms ||
        ghostAtoms->GetNumberOfGhosts() != nSize - nAtoms) {
      reallocationRequired = true;
      Allocate();
    }
  } else if (updateRequired) {
    // ghostAtoms->GetGhostPotential()->CheckAndUpdateAtoms has called
    // atoms->NormalizePositions, but for simulations without ghosts it must
    // also be done.
    atoms->NormalizePositions();
  }
  // update ghost positions:
  if (ghostAtoms != 0)
    ghostAtoms->GetGhostPotential()->UpdateGhostPositions();

  if (updateRequired)
    MakeList();
  UpdateImagePositions();
  return reallocationRequired;
}

void NeighborList::UpdateImagePositions() {
  const Vec *translations = atoms->GetSuperCell()->translations;
  for (unsigned int n = 0; n < images.size(); n++)
    positions[nSize + n] =
        positions[images[n].number] + translations[images[n].nTranslation];
}

int NeighborList::GetNeighbors(int a1, int *neighbors, Vec *diffs,
                               double *diffs2, int &size, double r) const {
  int i0 = indices[a1];
  int i1 = indices[a1 + 1];
  if (i1 - i0 > size) {
    cerr << "NBLIST OVERRUN" << endl << flush;
    assert(i1 - i0 > maxlistlen);
    size = -1;
    return 0;
  }
  double rC2 = rCut2;
  if (r > 0.0)
    rC2 = r * r;
  Vec pos1 = positions[a1];
  int nNeighbors = 0;
  for (int i = i0; i < i1; i++) {
    int a2 = neighborList[i];
    diffs[nNeighbors] = positions[a2] - pos1;
    double d2 = Length2(diffs[nNeighbors]);
    if (d2 < rC2) {
      diffs2[nNeighbors] = d2;
      if (a2 < nSize)
        neighbors[nNeighbors] = a2;
      else
        neighbors[nNeighbors] = images[a2 - nSize].number;
      nNeighbors++;
    }
  }
  size -= nNeighbors;
  assert(size >= 0);
  return nNeighbors;
}

const int *NeighborList::GetList(int a, int &n) {
  n = indices[a + 1] - indices[a];
  return neighborList + indices[a];
}

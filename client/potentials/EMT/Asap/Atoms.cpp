#include "Atoms.h"
#include "Exception.h"
#include "Potential.h"
#include "SuperCell.h"
// #include "Timing.h"
#include <assert.h>
#include <math.h>
#include <string.h>

extern int verbose;

#if verbose
#define VERB(x)                                                                \
  if (verbose == 1)                                                            \
  std::cerr << x
#else
#define VERB(x)
#endif

Atoms::Atoms(const Vec *p, int n, SuperCell *s)
    : potential(NULL), superCell(s), nAtoms(n), counter(1) {
  positions.resize(nAtoms);
  pos_translations.resize(nAtoms);
  memset(&pos_translations[0], 0, nAtoms * 3 * sizeof(double));
  types.resize(nAtoms);
  SetCartesianPositions(p);
  nImages = 0;
  superCell->SetAtoms(this);
}

/// Change the number of atoms, and resize the arrays containing the
/// data.  If the number of atoms is increased, the end of the arrays
/// will contain garbage.  Call both SetCartesianPositions() and
/// SetAtomicNumbers() after calling this method, and remember to
/// update any relevent Python arrays.  Mainly for use by
/// ParallelAtoms when the atoms migrate between processors.
void Atoms::SetNumberOfAtoms(int n) {
  nAtoms = n;
  positions.resize(nAtoms);
  pos_translations.resize(nAtoms);
  types.resize(nAtoms);
  MarkChanged();
  VERB("\nSetNumber");
}

void Atoms::SetNumberOfImages(int i) {
  positions.resize(positions.size() - nImages + i);
  nImages = i;
}

/// The cartesian positions are copied from the array of 3-vectors
/// (Vec) to the internal array.  No length checks are done.
void Atoms::SetCartesianPositions(const Vec *p) {
  memcpy(&positions[0], p, nAtoms * 3 * sizeof(double));
  MarkChanged();
  VERB("\nSetPos");
}

void Atoms::GetCartesianPositions(Vec *p) const {
  memcpy(p, &positions[0], nAtoms * 3 * sizeof(double));
}

void Atoms::GetUnwrappedPositions(vector<Vec> &p) const {
  p.resize(nAtoms);
  GetCartesianPositions(&p[0]);
  for (int i = 0; i < nAtoms; i++)
    p[i] -= pos_translations[i];
}

void Atoms::SetUnwrappedPositions(const Vec *p) {
  SetCartesianPositions(p);
  for (int i = 0; i < nAtoms; i++)
    positions[i] += pos_translations[i];
}

const Vec *Atoms::GetUnitCell() const { return superCell->vectors; }

void Atoms::SetUnitCell(const Vec newbasis[3], bool fix) {
  if (fix) {
    // Do not move the atoms, as seen from Python.
    for (int i = 0; i < nAtoms; i++) {
      positions[i] += pos_translations[i];
      pos_translations[i][0] = pos_translations[i][1] = pos_translations[i][2] =
          0.0;
    }
    superCell->SetBasis(newbasis);
    NormalizePositions();
  } else {
    // Keep the scaled positions fixed, allowing the real positions to move.
    // This will not influence the normalization.
    superCell->TransformToScaledSpace(&positions[0], nAtoms);
    superCell->TransformToScaledSpace(&pos_translations[0], nAtoms);
    superCell->SetBasis(newbasis);
    superCell->TransformToRealSpace(&pos_translations[0], nAtoms);
    superCell->TransformToRealSpace(&positions[0], nAtoms);
  }
  MarkChanged();
  if (GetPotential())
    GetPotential()->UpdateSuperCell(superCell);
}

///
/// Normalize an array of vectors containing differences in positions
/// using periodic boundary conditions and the minimum image
/// convention.  If there are free boundary conditions in a specific
/// direction, nothing is done to that coordinate of the vectors.
void Atoms::NormalizeDifferences(Vec *diff) const {
  //  USETIMER("Atoms::NormalizeDifferences");
  superCell->TransformToScaledSpace(diff, GetNumberOfAtoms());
  for (int i = 0; i < 3; i++)
    if (superCell->periodic[i]) {
      Vec *a = diff;
      for (int j = 0; j < nAtoms; ++j, ++a)
        (*a)[i] -= floor((*a)[i] + 0.5);
    }
  superCell->TransformToRealSpace(diff, GetNumberOfAtoms());
}

void Atoms::NormalizeDifference(Vec &diff) const {
  int i;
  Vec scaled;
  for (i = 0; i < 3; i++) {
    scaled[i] = superCell->inverse[i] * diff;
    if (superCell->periodic[i])
      scaled[i] -= floor(scaled[i] + 0.5);
  }
  diff = superCell->vectors[0] * scaled[0];
  for (i = 1; i < 3; i++)
    diff += superCell->vectors[i] * scaled[i];
}

void Atoms::NormalizePositions() {
  //  USETIMER("Atoms::NormalizePositions");
  if (!superCell->periodic[0] && !superCell->periodic[1] &&
      !superCell->periodic[2])
    return;
  superCell->TransformToScaledSpace(&positions[0], nAtoms);
  superCell->TransformToScaledSpace(&pos_translations[0], nAtoms);
  for (int i = 0; i < 3; i++)
    if (superCell->periodic[i]) {
      std::vector<Vec>::iterator a = positions.begin();
      std::vector<Vec>::iterator b = pos_translations.begin();
      for (int j = 0; j < nAtoms; ++j, ++a, ++b) {
        double d = -floor((*a)[i]);
        (*a)[i] += d;
        (*b)[i] += d;
      }
    }
  superCell->TransformToRealSpace(&pos_translations[0], nAtoms);
  superCell->TransformToRealSpace(&positions[0], nAtoms);
}

void Atoms::NormalizePosition(Vec &pos) const {
  Vec dummy;
  NormalizePosition(pos, dummy);
}

void Atoms::NormalizePosition(Vec &pos, Vec &scaled_translation) const {
  Vec scaled;
  int i;
  for (i = 0; i < 3; i++) {
    scaled[i] = superCell->inverse[i] * pos;
    if (superCell->periodic[i]) {
      double d = -floor(scaled[i]);
      scaled[i] += d;
      scaled_translation[i] = d;
    } else
      scaled_translation[i] = 0.0;
  }
  pos = superCell->vectors[0] * scaled[0];
  for (i = 1; i < 3; i++)
    pos += superCell->vectors[i] * scaled[i];
}

void Atoms::ReNormalizePosition(Vec &pos, Vec &scaled_translation) const {
  for (int i = 0; i < 3; i++)
    pos += superCell->vectors[i] * scaled_translation[i];
}

void Atoms::SetAtomicNumbers(const int *z) {
  memcpy(&types[0], z, nAtoms * sizeof(int));
  MarkChanged();
  VERB("\nSetZ");
}

#if 0
void Atoms::DeformAtoms(SuperCell *newSuperCell)
{
  for (int i = 0; i < 3; i++)
    if (superCell->periodic[i] != newSuperCell->periodic[i])
      throw (Exception("Periodicity can not be changed! p[") << "i" << "]: "
	     << superCell->periodic[i] << " -> " << newSuperCell->periodic[i]);
  potential->UpdateSuperCell(newSuperCell);
  for (int a = 0; a < nAtoms; a++)
    {
      Vec scaled;
      for ( i = 0; i < 3; i++)
	scaled[i] = superCell->inverse[i] * positions[a];
      positions[a] = newSuperCell->vectors[0] * scaled[0];
      for ( i = 1; i < 3; i++)
	positions[a] += newSuperCell->vectors[i] * scaled[i];
    }
  superCell = newSuperCell;
}
#endif

void Atoms::SetCalculator(AsapPotential *potential) {
  this->potential = potential;
  potential->SetAtoms(this);
}

void Atoms::GetListOfElements(set<int> &elements) const {
  const int *atomicnumbers = GetAtomicNumbers();

  elements.clear();
  for (int i = 0; i < nAtoms; i++) {
    int z = atomicnumbers[i];
    if (elements.find(z) == elements.end())
      elements.insert(z);
  }
}

#include "SuperCell.h"
#include "Atoms.h"
#include "Exception.h"
#include "Potential.h"
#include <iostream>
#include <math.h>

SuperCell::SuperCell(const Vec v[3], const bool p[3]) {
  atoms = 0;
  for (int i = 0; i < 3; i++)
    periodic[i] = p[i];

  SetBasis(v);
}

#if 0 // Not needed
SuperCell::SuperCell(const SuperCell &original)
{
  atoms = 0;
  for (i = 0; i < 3; i++)
    periodic[i] = original.periodic[i];

  SetBasis(original.vectors);
}
#endif

SuperCell::SuperCell(const SuperCell *original) {
  atoms = 0;
  for (int i = 0; i < 3; i++)
    periodic[i] = original->periodic[i];

  SetBasis(original->vectors);
}

void SuperCell::SetBasis(const Vec v[3]) {
  int i;
  for (i = 0; i < 3; i++)
    vectors[i] = v[i];
  double determinant = Cross(vectors[0], vectors[1]) * vectors[2];
  for (i = 0; i < 3; i++) {
    inverse[i] = Cross(vectors[(i + 1) % 3], vectors[(i + 2) % 3]);
    heights[i] = fabs(determinant) / sqrt(Length2(inverse[i]));
    inverse[i] /= determinant;
  }
  int s = 0;
  int k[3];
  for (k[2] = -1; k[2] <= 1; k[2]++)
    for (k[1] = -1; k[1] <= 1; k[1]++)
      for (k[0] = -1; k[0] <= 1; k[0]++) {
        Vec vec(0.0, 0.0, 0.0);
        for (i = 0; i < 3; i++)
          vec += vectors[i] * k[i];
        translations[s++] = vec;
      }
}

void SuperCell::SetAtoms(Atoms *a) {
  if (atoms && atoms != a)
    throw Exception("A SuperCell can only belong to one Atoms object.");
  else
    atoms = a;
}

void SuperCell::TransformToScaledSpace(Vec *positions, int nAtoms) const {
  for (int a = 0; a < nAtoms; a++) {
    Vec scaled;
    for (int i = 0; i < 3; i++)
      scaled[i] = inverse[i] * positions[a];
    positions[a] = scaled;
  }
}

void SuperCell::TransformToRealSpace(Vec *positions, int nAtoms) const {
  for (int a = 0; a < nAtoms; a++) {
    Vec scaled = positions[a];
    positions[a] = vectors[0] * scaled[0];
    for (int i = 1; i < 3; i++)
      positions[a] += vectors[i] * scaled[i];
  }
}

void SuperCell::ClosestLatticeSiteIndices(int intVector[3], Vec &diff) const {
  for (int i = 0; i < 3; i++)
    intVector[i] = int(floor(diff * inverse[i] + 0.5));
}

double SuperCell::GetVolume() const {
  double det;
  det = -vectors[0][2] * vectors[1][1] * vectors[2][0] +
        vectors[0][1] * vectors[1][2] * vectors[2][0] +
        vectors[0][2] * vectors[1][0] * vectors[2][1] -
        vectors[0][0] * vectors[1][2] * vectors[2][1] -
        vectors[0][1] * vectors[1][0] * vectors[2][2] +
        vectors[0][0] * vectors[1][1] * vectors[2][2];
  return fabs(det);
}

void SuperCell::SetPeriodic(const bool newperiodic[3]) {
  for (int i = 0; i < 3; i++)
    periodic[i] = newperiodic[i];

  if (atoms->GetPotential())
    atoms->GetPotential()->UpdateSuperCell(this);
}

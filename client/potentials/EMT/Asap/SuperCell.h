// Emacs, this is -*- C++ -*-

#ifndef _SUPERCELL_H
#define _SUPERCELL_H

#include "Vec.h"
class Atoms;

// vectors: supercell basis vectors
// inverse: vectors[i] * inverse[j] = ((i = j) ? 1 : 0)
// heights: heights of cell

/// The SuperCell defines the size, shape and boundary condx of a simulation.

/// The SuperCell serves two purposes.  The most important is the
/// computational box of a simulation: it defines the size, shape and
/// boundary conditions of the simulation.  A second use is for the
/// unit cell of the lattice in the elements of a Quasicontinuum
/// simulation.
class SuperCell {
public:
  /// Create a new supercell with basis vectors Vec and periodicity p.
  SuperCell(const Vec v[3], const bool p[3]);

  /// Copy constructor.
  SuperCell(const SuperCell *original);

  /// Set the atoms associated with the supercell.
  void SetAtoms(Atoms *a);

  /// In-place transformation of positions from real to scaled space.
  void TransformToScaledSpace(Vec *positions, int nAtoms) const;

  /// In-place transformation of positions from scaled to real space.
  void TransformToRealSpace(Vec *positions, int nAtoms) const;

  /// Quasicontinuum only: convert position to closest lattice site.

  /// When the supercell is used in a QC simulation to represent the
  /// crystal lattice of an element, this function is used to convert
  /// the position of an atom into lattice coordinates.
  void ClosestLatticeSiteIndices(int intVector[3], Vec &diff) const;

  /// Get the volume of the supercell.
  double GetVolume() const;

  /// Change the basis.  The atoms are not moved.  Should only be called from
  /// the atoms.
  void SetBasis(const Vec v[3]);

  /// Get the basis vectors
  Vec GetBasis(int a) { return vectors[a]; }

  /// Set the periodicity.  USE WITH CAUTION!

  /// Change the periodicity of the supercell.  Other classes may not
  /// notice that this has happened.  This cannot be done after the
  /// atoms have been assigned a Potential.
  void SetPeriodic(const bool newperiodic[3]);

public:
  bool periodic[3];  ///< Periodic boundary conditions along the 3 axes?
  Vec vectors[3];    ///< The basis vectors
  Vec inverse[3];    ///< The inverse of the vectors matrix.
  double heights[3]; ///< The heights of the supercell.
  Vec translations[27];
  ///< The 27 possible translations due to the periodic boundaries.
private:
  Atoms *atoms;
  ///< A pointer to the atoms associated with this supercell, if any.
};

#endif // _SUPERCELL_H

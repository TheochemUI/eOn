#ifndef ATOMS_H
#define ATOMS_H

#include <set>
#include <vector>
using std::set;
using std::vector;
#include "Vec.h" // Needed because we are destroying a vector<Vec> in
                 // Atoms::~Atoms().  (at least with the alpha compiler)

class SuperCell;
// class Vec;
class AsapPotential;

// Atomic positions are stored in a vector<Vec>:
//
// positions[0, ..., nAtoms - 1] :
//   real atoms
//
// positions[nAtoms, ..., nAtoms + nGhosts - 1] :
//   ghost atoms (from QC or Parallel)
//
// positions[nAtoms + nGhosts, ..., positions.size() - 1] :
//   image atoms (from periodic boundaries)

/// The main list of atoms class.

/// The Atoms class stores the coordinates and atomic numbers of the
/// atoms.  Coordinates are stored as Cartesian coordinates.  Atoms
/// has a reference to a SuperCell object containing the unitcell of
/// the simulation, and a pointer to a Potential object responsible of
/// force calculations.  Other information about atoms (such as
/// velocities and forces) are not stored in the C++ code but live in
/// Python.

class Atoms {
public:
  /// \brief Construct Atoms from array of positions, number of elements and
  /// a SuperCell.
  Atoms(const Vec *p, int n, SuperCell *s);

  /// Deletes the atoms.  Does not delete Potential or SuperCell.
  virtual ~Atoms(){};

  /// Return a pointer to the SuperCell object.
  virtual SuperCell *GetSuperCell() const { return superCell; }

  /// Return the unitcell vectors.
  virtual const Vec *GetUnitCell() const;

  /// Set a new unit cell, rescaling the positions unless fix is true
  virtual void SetUnitCell(const Vec newbasis[3], bool fix);

  /// Set a Potential.  An eventual old Potential is not deallocated.
  virtual void SetCalculator(AsapPotential *potential);

  /// Get a pointer to the Potential.
  AsapPotential *GetPotential() const { return potential; }

  /// Update the change pointer.

  /// The change counter is updated whenever the atoms are changed.  It is used
  /// to decide if e.g. forces should be recalculated.  Methods such as
  /// SetCartesianPositions() increment the counter.
  inline void MarkChanged() { ++counter; }

  /// Get the change counter.
  inline int GetChangeCounter() const { return counter; }

  /// Set the Cartesian positions.

  /// \warning When atoms.SetCartesianPositions() is called from
  /// Python, that is translated into a call to
  /// SetUwrappedPositions(), *not* to a call to this method.
  virtual void SetCartesianPositions(const Vec *p);

  /// Return a copy of the Cartesian positions.

  /// This version returns all degrees of freedom, i.e. both atoms and
  /// nodes in a quasicontinuum simulation.  It does not return ghost
  /// atoms.
  ///
  /// \warning When atoms.GetCartesianPositions() is called from
  /// Python, that is translated into a call to
  /// GetUwrappedPositions(), *not* to a call to this method.
  virtual void GetCartesianPositions(Vec *p) const;

  /// Return a const pointer to the Cartesian positions.

  /// This version returns a pointer to the coordinates of all atoms.
  /// If ghost atoms are present, they are at the end of the array.
  /// It is up to the caller to know if ghost atoms are present, and
  /// to access them (or not) according to the need.  Non-atomistic
  /// degrees of freedom, such as nodes in QC simulations, are not
  /// returned by this version.
  ///
  /// \bug The difference in semantics between the two versions of
  /// GetCartesianPositions is most unfortunate: one of them should be
  /// renamed.
  virtual const Vec *GetCartesianPositions() const { return &positions[0]; }

  /// Return a non-const pointer to the Cartesian positions.  DANGEROUS!

  /// GetPositionsPtr returns a pointer allowing modification of the
  /// positions.  Be sure to call MarkChanged if this is done.
  Vec *GetPositionsPtr() { return &positions[0]; }

  /// \brief Return the positions, hiding that they may have been
  /// wrapped due to the boundary conditions.

  /// This is the function used by Python when GetCartesianPositions
  /// is called.
  virtual void GetUnwrappedPositions(vector<Vec> &p) const;

  /// \brief Set the positions from an array of positions not taking
  /// any wrapping into account.

  /// The already existing offset (caused by wrapping through periodic
  /// boundaries) is applied to the atoms.  This is the function used
  /// by Python when SetCartesianPositions is called.
  virtual void SetUnwrappedPositions(const Vec *p);

  /// Apply periodic boundary conditions to differences in positions.
  virtual void NormalizeDifferences(Vec *diff) const;

  /// Apply periodic boundary conditions to difference in positions.
  virtual void NormalizeDifference(Vec &diff) const;

  /// \brief Apply periodic boundary conditions to the positions.
  /// Does NOT call MarkChanged !

  /// Normalize the positions, i.e. translate the atoms back into the
  /// supercell if there are periodic boundary conditions.  This does
  /// not count as changing the atoms, so the change counter is NOT
  /// updated.  The actual work is done by the SuperCell.  This method
  /// should be called by the Potential, which may (should?) delegate
  /// it to the NeighborList.
  virtual void NormalizePositions();

  /// Normalize a position which is not owned by the atoms.

  /// If another object holds a position needing to be normalized,
  /// this function should be called.  It exists in two versions, one
  /// just normalizing the position, and one also returning the scaled
  /// space correction applied.
  virtual void NormalizePosition(Vec &pos) const;

  /// Normalize a position which is not owned by the atoms.

  /// This version of the call reports the translation in
  /// scaled_translation, so the same translation can be obtained
  /// later (by calling ReNormalizePosition) without redoing the
  /// decision on whether to wrap.  This is useful if positions need
  /// to be normalized between neighborlist updates, and is currently
  /// used by the Quasicontinuum code.
  virtual void NormalizePosition(Vec &pos, Vec &scaled_translation) const;

  /// Repeat a normalization.  See NormalizePosition(a,b).
  virtual void ReNormalizePosition(Vec &pos, Vec &scaled_translation) const;

  /// Set the atomic numbers.  Should be done right after construction.
  virtual void SetAtomicNumbers(const int *z);

  /// Return a const pointer to the atomic numbers.
  virtual const int *GetAtomicNumbers() const { return &types[0]; }

  /// Return a non-const pointer to the atomic numbers.  DANGEROUS!

  /// GetAtomicNumbersPtr returns a pointer allowing modification of the atomic
  /// numbers.  Use with care, and call MarkChanged manually!
  int *GetAtomicNumbersPtr() { return &types[0]; }

  /// Get a set of all elements present in the simulations.
  virtual void GetListOfElements(set<int> &elements) const;

  /// Set the number of images.  Used by the neighbor list.
  void SetNumberOfImages(int i);

  /// Get the number of images.  Used by the neighbor list.
  int GetNumberOfImages() const { return nImages; }

  /// Get the number of atoms (including nodes in a QC simulation).

  /// In a QuasiContinuum calculation we have real atoms and node atoms.
  /// GetNumberOfAtoms() is a virtual function: For a QCAtoms class it will
  /// return the number of real atoms plus the number of node atoms. This
  /// number is needed by SWIG to get the size of arrays. The method
  /// GetNumberOfRealAtoms() is needed by potentials and the Neighbor List
  /// object (they only deal with the real atoms).
  virtual int GetNumberOfAtoms() const { return nAtoms; }

  /// Get the number of atoms (excluding nodes in a QC simulation).

  /// See also the documentation of GetNumberOfAtoms().
  ///
  int GetNumberOfRealAtoms() const { return nAtoms; }

  /// Change the number of atoms.
  virtual void SetNumberOfAtoms(int n);

protected:
  AsapPotential *potential;     ///< A pointer to the potential.
  SuperCell *superCell;         ///< A pointer to the supercell.
  vector<Vec> positions;        ///< Contains the Cartesian positions.
  vector<Vec> pos_translations; ///< \brief Translations applied to the
                                ///  positions to conform with the periodic
                                ///  boundary conditions
  vector<int> types;            ///< Contains the atomic numbers.
  int nAtoms;                   ///< The number of atoms in the list.
  int nImages; ///< The number of image atoms used by NeighborList.
  int counter; ///< Used to track changes.
};

#endif // ATOMS_H

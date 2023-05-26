// Emacs: this is -*- C++ -*-

#ifndef ASAP_POTENTIAL_H
#define ASAP_POTENTIAL_H

#include "Exception.h"
#include <assert.h>

class Atoms;
class Vec;
class SuperCell;

typedef double symTensor[6];

/// Abstract base class of all potentials.

/// A Potential calculates forces, energies and stresses for a list of
/// atoms.  A given instance of a Potential is associated with a
/// specific Atoms object on a one-to-one bases, this is established
/// when Atoms.SetCalculator() is called.
///
/// Two types of classes are derived from Potential.
///   - Concrete implementations of potentials, such as EMT
///     and MoPotential.
///   - Potentials wrapping concrete implementations, but providing
///     special functionalily, such as ParallelPotential (implementing
///     parallel simulations) or QCPotential (implementing the
///     QuasiContinuum method).
class AsapPotential {
public:
  virtual ~AsapPotential() {}

  /// Set the atoms belonging to this potential.

  ///
  /// This is called automatically by Atoms.SetCalculator() and should
  /// not be called elsewhere.
  virtual void SetAtoms(Atoms *a) = 0;

  /// Calculate the total energy of the system.
  virtual double GetPotentialEnergy() = 0;

  /// Calculate the forces on all atoms and return the result.
  virtual const Vec *GetCartesianForces() = 0;

  /// Calculate the stress on all atoms.
  virtual const symTensor *GetStresses(const Vec *momenta = 0) = 0;

  /// Calculate the total stress of the system.

  /// Note that the output variable stress is not zeroed first.  In
  /// quasicontinuum simulations and other situations where there is a
  /// contribution to the stress which is not from the atoms, the
  /// stress parameter can initially be set to the appropriate
  /// derivatives of the energy of that part of the calculation
  /// (without the volume normalization, which will be performed by
  /// the potential).
  virtual void GetStress(double stress[6], const Vec *momenta = 0) = 0;

  /// Calculate the energy of all atoms.
  virtual const double *GetPotentialEnergies() = 0;

  virtual void CheckNeighborLists() = 0;
  /// Calculate the energy of an atom in a regular fcc(?) lattice

  /// Given the lattice vectors, this function calculates the energy
  /// of a single atom in the lattice.  Used for the QuasiContinuum method.
  virtual double CalculateLatticeEnergy(const Vec a[3]) {
    throw Exception("The QC method is not implemented for this potential");
    return 0.0;
  }

  /// Calculate derivative of the energy of an atom in a regular lattice.

  ///
  /// Used by the QuasiContinuum method.
  virtual void CalculateDerivatives(const Vec a[3], double dEdaDota[6]) {
    throw Exception("The QC method is not implemented for this potential");
  }

  /// Get data used by the QuasiContinuum method.
  virtual double GetData() const {
    throw Exception("GetData() not implemented.");
    return 0.0;
  }

  /// Return the cutoff radius used in the potential.
  virtual double GetCutoffRadius() const = 0;

  /// Return the lattice constant of the material, if well-defined.

  /// If a lattice constant of the material can be defined, return it
  /// in Angstrom, otherwise throw an exception.
  virtual double GetLatticeConstant() const = 0;

  /// Tell the Potential that a new SuperCell has been assigned to the atoms.
  virtual void UpdateSuperCell(const SuperCell *newSuperCell) = 0;

  /// Get the number of atoms.
  virtual int GetNumberOfAtoms() const = 0;
};

#endif // POTENTIAL_H

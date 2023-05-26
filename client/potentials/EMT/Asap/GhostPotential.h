// Emacs: This is -*- C++ -*-

#ifndef GHOSTPOTENTIAL_H
#define GHOSTPOTENTIAL_H

#include "Potential.h"

/// Common base class of ParallelPotential and QCPotential.

/// GhostPotential is a common base class of the functional extensions
/// of Potential for extensions of Atoms which have the concept of
/// ghost atoms (ParallelAtoms and QCAtoms).
///
/// A GhostPotential is derived from Potential since it replaces the
/// normal potential in parallel/QC simulations.  It also wraps a
/// Potential, namely a specific implementation of Potential (EMT,
/// MoPotential ...) describing the interatomic interactions.  In this
/// way both parallel and QC simulations are supported with minimal
/// intrusion in the potentials (the specific potentials occationally
/// need to check if the Atoms are really GhostAtoms and then call a
/// few helper functions).

class GhostPotential : public AsapPotential {
public:
  GhostPotential(AsapPotential *p) : potential(p) {}

  /// Check if the ghost atoms want the neighbor list updated.

  /// Called by NeihborList::CheckAndUpdateNeighborList() to
  /// post-process the decision to update the neighbor list.  Only
  /// used by ParallelPotential, where all neighbor lists on all
  /// processors should be updated when one processor does it.
  virtual bool CheckAndUpdateAtoms(bool update) = 0;

  /// Place data on virtual atoms.

  /// CommunicateData is used by potentials where exactly one double
  /// per atom needs to be updated on all atoms.  The data may
  /// originate from the corresponding arrays on other processors
  /// (this is the case in ParallelPotential) or be calculated
  /// (QCPotential).  In EMT this quantity is sigma1.  If future
  /// potentials need to communicate more date, this interface must be
  /// extended.
  virtual void CommunicateData(double data[]) = 0;
  /// Place data on virtual atoms.

  /// \bug This is a primitive extension of CommunicateData(double *),
  /// it does not work with the quasicontinuum method and must be
  /// improved.
  virtual void CommunicateData(double data[], int n) = 0;

  /// Update the position of ghost atoms.
  virtual void UpdateGhostPositions() = 0;

  /// Get the cutoff of the wrapped potential.
  double GetCutoffRadius() const { return potential->GetCutoffRadius(); }

  /// Get the lattice constant of the wrapped potential.
  double GetLatticeConstant() const { return potential->GetLatticeConstant(); }

protected:
  AsapPotential *potential; ///< Pointer to the wrapped potential.
};

#endif // GHOSTPOTENTIAL_H

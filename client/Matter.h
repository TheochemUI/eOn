/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once
#include "ConFileIO.h"
#include "Eigen.h"
#include "EonLogger.h"
#include "Parameters.h"
#include "Potential.h"
#include "SurrogatePotential.h"
#include <array>
#include <memory>
#include <string>

// This is a forward declaration of BondBoost to avoid a circular dependency.
namespace eonc {
class BondBoost;

namespace pbc {

inline AtomMatrix apply(const AtomMatrix &diff, const Matrix3d &cell,
                        const Matrix3d &cellInverse) {
  // Transform to fractional coordinates, wrap to [-0.5, 0.5), transform back.
  // Uses floor(x + 0.5) instead of double-fmod: single x86 vroundsd instruction
  // vs expensive fmod library call. Vectorized via Eigen .array() operations.
  AtomMatrix frac = diff * cellInverse;
  frac.array() -= (frac.array() + 0.5).floor();
  return frac * cell;
}

inline VectorXd applyV(const VectorXd &diffVector, const Matrix3d &cell,
                       const Matrix3d &cellInverse) {
  AtomMatrix pbcMatrix =
      apply(AtomMatrix::Map(diffVector.data(), diffVector.size() / 3, 3), cell,
            cellInverse);
  return VectorXd(VectorXd::Map(pbcMatrix.data(), diffVector.size()));
}

} // namespace pbc

/* Data describing an atomic structure. This class has been devised to handle
 * information about an atomic structure such as positions, velocities, masses,
 * etc. It also allow to associate a forcefield for the structure through a
 * pointer to function (potential()). The class can read and save data to a .con
 * file (atom2con() and con2atom()). It can also save to a .xyz file
 * (atom2xyz()).*/

class Matter {
public:
  ~Matter() = default;
  Matter(std::shared_ptr<Potential> pot, const Parameters &params)
      : potential{pot},
        usePeriodicBoundaries{true},
        recomputePotential{true},
        forceCalls{0},
        removeNetForce{params.main_options.removeNetForce},
        structComp{params.structure_comparison_options},
        parameters{&params},
        nAtoms{0},
        positions{MatrixXd::Zero(0, 3)},
        velocities{MatrixXd::Zero(0, 3)},
        forces{MatrixXd::Zero(0, 3)},
        biasForces{MatrixXd::Zero(0, 3)},
        biasPotential{nullptr},
        masses{Eigen::VectorXd::Zero(0)},
        atomicNrs{Eigen::VectorXi::Zero(0)},
        isFixed{Eigen::VectorXi::Zero(0)},
        cell{Matrix3d::Zero()},
        cellInverse{Matrix3d::Zero()},
        energyVariance{0.0},
        potentialEnergy{0.0} {} // the number of atoms shall be set later
  // using resize()
  Matter(const Matter &matter);                  // create a copy of matter
  const Matter &operator=(const Matter &matter); // copy the matter object
  bool compare(const Matter &matter, bool indistinguishable = false);

  double
  distanceTo(const Matter &matter); // the distance to the given matter object
  double perAtomNorm(const Matter &matter); // the maximum distance between two
                                            // atoms in the Matter objects
  void
  setPotential(std::shared_ptr<Potential> pot); // set potential function to use
  std::shared_ptr<Potential> getPotential();    // get potential function to use
  void resize(long int nAtoms);   // set or reset the number of atoms
  long int numberOfAtoms() const; // return the number of atoms
  Matrix3d getCell() const;
  void setCell(const Matrix3d &newCell);
  double getPosition(long int atom, int axis)
      const; // return the position of an atom along one of the axis
  void setPosition(
      long int atom, int axis,
      double position); // set the position of atom along axis to position
  void setVelocity(
      long int atom, int axis,
      double velocity); // set the velocity of atom along axis to velocity
  bool relax(bool quiet = false, bool writeMovie = false,
             bool checkpoint = false, std::string prefixMovie = std::string(),
             std::string prefixCheckpoint = std::string());

  AtomMatrix pbc(const AtomMatrix &diff) const {
    return eonc::pbc::apply(diff, cell, cellInverse);
  }
  VectorXd pbcV(const VectorXd &diff) const {
    return eonc::pbc::applyV(diff, cell, cellInverse);
  }

  size_t getPotentialCalls() const;
  const AtomMatrix &getPositions() const; // return coordinates of atoms
  AtomMatrix getPositionsCopy() const;    // return a modifiable copy
  VectorXd getPositionsV() const;
  AtomMatrix
  getPositionsFree() const; // return coordinates of free atoms in array pos
  VectorXd getPositionsFreeV() const;
  void
  setPositions(const AtomMatrix &pos); // update Matter with the new positions
                                       // of the free atoms given in array pos
  void setPositionsV(const VectorXd &pos);
  void setPositionsFree(
      const AtomMatrix &pos); // update Matter with the new positions of the
                              // free atoms given in array pos
  void setPositionsFreeV(const VectorXd &pos);

  AtomMatrix getVelocities() const;
  void setVelocities(const AtomMatrix &v);
  void setBiasForces(const AtomMatrix &bf);
  void setBiasPotential(BondBoost *bondBoost);
  void setForces(const AtomMatrix &f);
  AtomMatrix getAccelerations();

  const AtomMatrix &getForces() const;
  const AtomMatrix &getForcesRaw() const;
  AtomMatrix getBiasForces();
  VectorXd getForcesV() const;
  AtomMatrix getForcesFree();
  VectorXd getForcesFreeV();

  double getMass(long int atom) const; // return the mass of the atom specified
  void setMass(long int atom, double mass); // set the mass of an atom
  void setMasses(const VectorXd &massesIn); // set the mass of an atom
  long getAtomicNr(
      long int atom) const; // return the atomic number of the atom specified
  void setAtomicNr(long int atom,
                   long atomicNr);           // set the atomic number of an atom
  VectorXi getAtomicNrs() const;             // Get the vector of atomic numbers
  VectorXi getAtomicNrsFree() const;         // Get the vector of atomic numbers
  void setAtomicNrs(const VectorXi &atmnrs); // set the vector of atomic numbers

  int getFixed(long int atom)
      const; // return true if the atom is fixed, false if it is movable
  void setFixed(long int atom,
                int isFixed); // set the atom to fixed (true) or movable (false)
  double getEnergyVariance() const;
  double getPotentialEnergy() const;
  double getKineticEnergy() const;
  double getMechanicalEnergy() const;

  double distance(long index1, long index2)
      const; // return the distance between two atoms in same configuration
  double pdistance(long index1, long index2, int axis) const;
  double distance(const Matter &matter, long index)
      const; // the distance between the same atom in two cofigurations

  long int
  numberOfFreeAtoms() const; // return the number of free (or movable) atoms
  long int numberOfFixedAtoms() const; // return the number of fixed atoms

  long
  getForceCalls() const; // return how many force calls that have been performed
  void resetForceCalls(); // zeroing the value of force calls

  double maxForce(void) const;

  // I/O delegates to eonc::io free functions
  void writeTibble(std::string filename) { io::writeTibble(*this, filename); }
  bool con2matter(std::string filename) {
    return io::con2matter(*this, filename);
  }
  bool con2matter(FILE *file) { return io::con2matter(*this, file); }
  bool convel2matter(std::string filename) {
    return io::convel2matter(*this, filename);
  }
  bool convel2matter(FILE *file) { return io::convel2matter(*this, file); }
  bool matter2con(std::string filename, bool append = false) {
    return io::matter2con(*this, filename, append);
  }
  bool matter2con(FILE *file) { return io::matter2con(*this, file); }
  bool matter2convel(std::string filename) {
    return io::matter2convel(*this, filename);
  }
  bool matter2convel(FILE *file) { return io::matter2convel(*this, file); }
  void matter2xyz(std::string filename, bool append = false) {
    io::matter2xyz(*this, filename, append);
  }

  AtomMatrix getFree() const;
  VectorXd getFreeV() const;
  Eigen::Matrix<double, Eigen::Dynamic, 1> getMasses() const;

private:
  // Friend declarations for eonc::io free functions that need private access
  friend bool io::con2matter(Matter &, FILE *);
  friend bool io::convel2matter(Matter &, FILE *);
  friend bool io::matter2con(Matter &, FILE *);
  friend bool io::matter2convel(Matter &, FILE *);
  friend void io::matter2xyz(Matter &, std::string, bool);

  eonc::log::Scoped m_log;
  std::shared_ptr<Potential>
      potential; // pointer to function calculating the energy and forces
  bool usePeriodicBoundaries; // boolean telling periodic boundaries are used
  mutable bool recomputePotential; // boolean indicating if the potential energy
                                   // and forces need to be recalculated
  mutable long
      forceCalls; // keep track of how many force calls have been performed

  // CON file header lines (indices 0-4 map to old headerCon1,2,4,5,6)
  std::array<std::string, 5> headerCon;

  void computePotential() const;
  void applyPeriodicBoundary();
  void applyPeriodicBoundary(double &component, int axis);
  void applyPeriodicBoundary(AtomMatrix &diff);

  // Narrowed from Parameters: only fields Matter actually reads
  bool removeNetForce{true};
  Parameters::structure_comparison_options_t structComp;
  // Full Parameters pointer retained solely for relax() delegation
  const Parameters *parameters;
  long nAtoms;
  AtomMatrix positions;
  AtomMatrix velocities;
  mutable AtomMatrix forces;
  AtomMatrix biasForces;
  BondBoost *biasPotential;
  VectorXd masses;
  VectorXi atomicNrs;
  VectorXi isFixed;   // array of bool, false for movable atom, true for fixed
  VectorXi atomIndex; // original atom index from .con column 5
  mutable AtomMatrix freeMask; // cached Nx3 mask (1.0 for free, 0.0 for fixed)
  mutable AtomMatrix maskedForces;      // cached forces with fixed atoms zeroed
  mutable std::vector<int> freeIndices; // cached indices of free atoms
  mutable bool recomputeFreeMask{true};
  mutable bool recomputeMaskedForces{true};
  Matrix3d cell;
  Matrix3d cellInverse;
  mutable double energyVariance;
  mutable double potentialEnergy;
};

} // namespace eonc

using eonc::Matter;

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
#include "Eigen.h"
#include "Parameters.h"
#include "Potential.h"
#include <memory>
#include <spdlog/sinks/basic_file_sink.h>
#include <string>

namespace eonc {
// This is a forward declaration of BondBoost to avoid a circular dependency.
class BondBoost;

// XXX: Really need this to go..
constexpr size_t MAXC{100}; // maximum number of components for functions
                            // matter2con and con2matter

/* Data describing an atomic structure. This class has been devised to handle
 * information about an atomic structure such as positions, velocities, masses,
 * etc. It also allow to associate a forcefield for the structure through a
 * pointer to function (potential()). The class can read and save data to a .con
 * file (atom2con() and con2atom()). It can also save to a .xyz file
 * (atom2xyz()).*/

class Matter {
public:
  ~Matter() = default;
  Matter(std::shared_ptr<Potential> pot)
      : potential{pot} {
    m_log = spdlog::get("combi");
  }
  Matter(const Matter &matter);                  // create a copy of matter
  const Matter &operator=(const Matter &matter); // copy the matter object
  // bool operator==(const Matter& matter); // true if differences in positions
  // are below differenceDistance bool operator!=(const Matter& matter); //
  // inverse of ==
  bool compare(const Matter &matter, bool indistinguishable = false);

  double
  distanceTo(const Matter &matter); // the distance to the given matter object
  double
  perAtomNorm(const Matter &matter) const; // the maximum distance between two
                                           // atoms in the Matter objects
  void
  setPotential(std::shared_ptr<Potential> pot); // set potential function to use
  std::shared_ptr<Potential> getPotential();    // get potential function to use
  void resize(long int nAtoms); // set or reset the number of atoms
  size_t numberOfAtoms() const; // return the number of atoms
  Matrix3S getCell() const;
  void setCell(Matrix3S newCell);
  double getPosition(long int atom, int axis)
      const; // return the position of an atom along one of the axis
  void setPosition(
      long int atom, int axis,
      double position); // set the position of atom along axis to position
  void setVelocity(
      long int atom, int axis,
      double velocity); // set the velocity of atom along axis to velocity
  AtomMatrix pbc(AtomMatrix diff) const;
  VectorType pbcV(VectorType diff) const;

  size_t getPotentialCalls() const;
  AtomMatrix getPositions() const; // return coordinates of atoms in array pos
  VectorType getPositionsV() const;
  AtomMatrix
  getPositionsFree() const; // return coordinates of free atoms in array pos
  VectorType getPositionsFreeV() const;
  void
  setPositions(const AtomMatrix pos); // update Matter with the new positions of
                                      // the free atoms given in array pos
  void setPositionsV(const VectorType pos);
  void setPositionsFree(
      const AtomMatrix pos); // update Matter with the new positions of the free
                             // atoms given in array pos
  void setPositionsFreeV(const VectorType pos);

  AtomMatrix getVelocities() const;
  void setVelocities(const AtomMatrix v);
  void setBiasForces(const AtomMatrix bf);
  void setBiasPotential(BondBoost *bondBoost);
  void setForces(const AtomMatrix f);
  AtomMatrix getAccelerations();

  AtomMatrix getForces(); // return forces applied on all atoms in array force
  AtomMatrix getBiasForces();
  VectorType getForcesV();
  AtomMatrix getForcesFree();
  VectorType getForcesFreeV();

  double getMass(long int atom) const; // return the mass of the atom specified
  void setMass(long int atom, double mass); // set the mass of an atom
  void setMasses(VectorType massesIn);      // set the mass of an atom
  size_t getAtomicNr(
      size_t atom) const; // return the atomic number of the atom specified
  void setAtomicNr(long int atom,
                   long atomicNr);         // set the atomic number of an atom
  Vector<size_t> getAtomicNrs() const;     // Get the vector of atomic numbers
  Vector<size_t> getAtomicNrsFree() const; // Get the vector of atomic numbers
  void
  setAtomicNrs(const Vector<size_t> atmnrs); // Get the vector of atomic numbers

  int getFixed(long int atom)
      const; // return true if the atom is fixed, false if it is movable
  void setFixed(long int atom,
                int isFixed); // set the atom to fixed (true) or movable (false)
  // void setPotentialEnergy(double);
  double getEnergyVariance();
  // VectorType getForceVariance();
  // double getMaxVariance();
  double getPotentialEnergy();
  double getKineticEnergy() const;
  double getMechanicalEnergy(); // return the mechanical energy (i.e. kinetic
                                // plus potential energy)

  double distance(long index1, long index2)
      const; // return the distance between two atoms in same configuration
  double pdistance(long index1, long index2, int axis) const;
  double distance(const Matter &matter, long index)
      const; // the distance between the same atom in two cofigurations

  long int
  numberOfFreeAtoms() const; // return the number of free (or movable) atoms
  long int numberOfFixedAtoms() const; // return the number of fixed atoms

  size_t
  getForceCalls() const; // return how many force calls that have been performed
  void resetForceCalls(); // zeroing the value of force calls

  double maxForce(void);

  void writeTibble(std::string filename);
  bool con2matter(std::string filename); // read con file into Matter, return
                                         // true if successful
  bool con2matter(FILE *file); // read con file and load data into Matter,
                               // return true if successful
  bool
  convel2matter(std::string filename); // read con file with both coordinates
                                       // and velocities into Matter
  bool convel2matter(FILE *file); // read con file with both coordinates and
                                  // velocities and load data into Matter
  bool
  matter2con(std::string filename,
             bool append = false); // print con file from data in Class Matter
  bool matter2con(FILE *file);     // print con file from data in Class Matter
  bool
  matter2convel(std::string filename); // print con file with both coordinates
                                       // and velocities  in Class Matter
  bool matter2convel(FILE *file); // print con file with both coordinates and
                                  // velocities from data in Class Matter
  void matter2xyz(std::string filename,
                  bool append = false); // print xyz file from data in Matter
  AtomMatrix getFree() const;
  VectorType getFreeV() const;
  VectorType getMasses() const;

private:
  std::shared_ptr<spdlog::logger> m_log;
  std::shared_ptr<Potential> potential;
  // --- These should be parsed in
  bool usePeriodicBoundaries{true};
  // Indicates if the potential energy and forces need to be recalculated
  bool recomputePotential{true};
  bool removeNetForce{true};
  // ----
  size_t forceCalls{0};

  // CON file header information, which is not used in the eon code
  char headerCon1[512];
  char headerCon2[512];
  char headerCon4[512];
  char headerCon5[512];
  char headerCon6[512];

  void computePotential();
  void initializeDataMembers(std::shared_ptr<Parameters> parameters);
  void applyPeriodicBoundary();
  void applyPeriodicBoundary(double &component, int axis);
  void applyPeriodicBoundary(AtomMatrix &diff);

  // Stuff which used to be in MatterPrivateData
  size_t nAtoms{0};
  AtomMatrix positions{MatrixType::Zero(0, 3)};
  AtomMatrix velocities{MatrixType::Zero(0, 3)};
  AtomMatrix forces{MatrixType::Zero(0, 3)};
  AtomMatrix biasForces{MatrixType::Zero(0, 3)};
  BondBoost *biasPotential{nullptr};
  VectorType masses{VectorType::Zero(0)};
  Vector<size_t> atomicNrs{Vector<size_t>::Zero(0)};
  Vector<int> isFixed{Vector<int>::Zero(
      0)}; // array of bool, false for movable atom, true for fixed
  Matrix3S cell{Matrix3S::Zero()};
  Matrix3S cellInverse{Matrix3S::Zero()};
  double energyVariance{0.0};
  double potentialEnergy{0.0};
};

} // namespace eonc

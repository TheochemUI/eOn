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
#include "client/Eigen.h"
#include "client/Parameters.h"
#include "client/Potential.h"
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
  Matter(std::shared_ptr<PotBase> pot)
      : potential{pot} {
    m_log = spdlog::get("combi");
  }
  Matter(const Matter &matter);                  // create a copy of matter
  const Matter &operator=(const Matter &matter); // copy the matter object
  // bool operator==(const Matter& matter); // true if differences in positions
  // are below differenceDistance bool operator!=(const Matter& matter); //
  // inverse of ==
  bool compare(const Matter &matter, bool indistinguishable = false);
  double distanceTo(const Matter &matter);
  double perAtomNorm(const Matter &matter) const;
  void setPotential(std::shared_ptr<PotBase> pot);
  std::shared_ptr<PotBase> getPotential();
  void resize(long int nAtoms);
  size_t numberOfAtoms() const;
  double getPosition(long int atom, int axis) const;
  void setPosition(long int atom, int axis, double position);
  void setVelocity(long int atom, int axis, double velocity);
  AtomMatrix pbc(AtomMatrix diff) const;
  VectorType pbcV(VectorType diff) const;

  size_t getPotentialCalls() const;
  AtomMatrix getPositions() const;
  VectorType getPositionsV() const;
  AtomMatrix getPositionsFree() const;
  VectorType getPositionsFreeV() const;
  void setPositions(const AtomMatrix pos);
  void setPositionsV(const VectorType pos);
  void setPositionsFree(const AtomMatrix pos);
  void setPositionsFreeV(const VectorType pos);

  AtomMatrix getVelocities() const;
  void setVelocities(const AtomMatrix v);
  void setBiasForces(const AtomMatrix bf);
  void setBiasPotential(BondBoost *bondBoost);
  void setForces(const AtomMatrix f);
  AtomMatrix getAccelerations();

  AtomMatrix getForces();
  AtomMatrix getBiasForces();
  VectorType getForcesV();
  AtomMatrix getForcesFree();
  VectorType getForcesFreeV();

  double getMass(long int atom) const;
  void setMass(long int atom, double mass);
  void setMasses(VectorType massesIn);
  size_t getAtomicNr(size_t atom) const;
  void setAtomicNr(long int atom, long atomicNr);
  Vector<size_t> getAtomicNrs() const;
  Vector<size_t> getAtomicNrsFree() const;
  void setAtomicNrs(const Vector<size_t> atmnrs);

  int getFixed(long int atom) const;
  // set the atom to fixed (true) or movable (false)
  void setFixed(long int atom, int isFixed);
  // void setPotentialEnergy(double);
  double getEnergyVariance();
  // VectorType getForceVariance();
  // double getMaxVariance();
  double getPotentialEnergy();
  double getKineticEnergy() const;
  // return the mechanical energy (i.e. kinetic
  // plus potential energy)
  double getMechanicalEnergy();

  double distance(long index1, long index2) const;
  double pdistance(long index1, long index2, int axis) const;
  double distance(const Matter &matter, long index) const;

  long int numberOfFreeAtoms() const;
  long int numberOfFixedAtoms() const;

  size_t getForceCalls() const;
  void resetForceCalls();

  double maxForce(void);

  bool con2matter(std::string filename);
  bool con2matter(FILE *file);
  bool convel2matter(std::string filename);
  bool convel2matter(FILE *file);
  AtomMatrix getFree() const;
  VectorType getFreeV() const;
  VectorType getMasses() const;

  Matrix3S cell{Matrix3S::Zero()};

  // --- These should be parsed in
  bool usePeriodicBoundaries{true};
  bool removeNetForce{true};
  // ----

private:
  std::shared_ptr<spdlog::logger> m_log;
  std::shared_ptr<PotBase> potential;
  // Indicates if the potential energy and forces need to be recalculated
  bool recomputePotential{true};
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
  size_t npotcalls{0};
  AtomMatrix positions{MatrixType::Zero(0, 3)};
  AtomMatrix velocities{MatrixType::Zero(0, 3)};
  AtomMatrix forces{MatrixType::Zero(0, 3)};
  AtomMatrix biasForces{MatrixType::Zero(0, 3)};
  BondBoost *biasPotential{nullptr};
  VectorType masses{VectorType::Zero(0)};
  Vector<size_t> atomicNrs{Vector<size_t>::Zero(0)};
  // array of bool, false for movable atom, true for fixed
  Vector<int> isFixed{Vector<int>::Zero(0)};
  Matrix3S cellInverse{Matrix3S::Zero()};
  double energyVariance{0.0};
  double potentialEnergy{0.0};
};

} // namespace eonc

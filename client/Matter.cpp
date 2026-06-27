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
#include "Matter.h"
#include "BaseStructures.h"
#include "BondBoost.h"
#include "GeometryAnalysis.h"
#include "HelperFunctions.h"
#include "SurrogatePotential.h"

#include "EonLogger.h"
#include <memory>
#include <stdexcept>

Matter::Matter(const Matter &matter) { operator=(matter); }

const Matter &Matter::operator=(const Matter &matter) {
  nAtoms = matter.nAtoms;
  resize(nAtoms);

  positions = matter.positions;
  forces = matter.forces;
  masses = matter.masses;
  atomicNrs = matter.atomicNrs;
  isFixed = matter.isFixed;
  atomIndex = matter.atomIndex;
  cell = matter.cell;
  cellInverse = matter.cellInverse;
  velocities = matter.velocities;

  removeNetForce = matter.removeNetForce;
  structComp = matter.structComp;
  parameters = matter.parameters;

  usePeriodicBoundaries = matter.usePeriodicBoundaries;
  pbcConvention = matter.pbcConvention;

  potential = matter.potential;
  potentialEnergy = matter.potentialEnergy;
  recomputePotential = matter.recomputePotential;
  recomputeFreeMask = true; // Force cache rebuild after assignment

  headerCon = matter.headerCon;

  return *this;
}

// The == comparison considers identity. This is crucial for process search.
// bool Matter::operator==(const Matter& matter) {
//     if(structComp.check_rotation) {
//         return eonc::helpers::rotationMatch(this, &matter,
//         structComp.distance_difference);
//     }else{
//         return (structComp.distance_difference)
//         > perAtomNorm(matter);
//     }
// }

bool Matter::compare(const Matter &matter, bool indistinguishable) {
  if (nAtoms != matter.numberOfAtoms())
    return false;
  if (structComp.check_rotation && indistinguishable) {
    return eonc::helpers::sortedR(*this, matter,
                                  structComp.distance_difference);
  } else if (indistinguishable) {
    if (this->numberOfFixedAtoms() == 0 and structComp.remove_translation)
      eonc::helpers::translationRemove(*this, matter);
    return eonc::helpers::identical(*this, matter,
                                    structComp.distance_difference);
  } else if (structComp.check_rotation) {
    return eonc::helpers::rotationMatch(*this, matter,
                                        structComp.distance_difference);
  } else {
    if (this->numberOfFixedAtoms() == 0 and structComp.remove_translation)
      eonc::helpers::translationRemove(*this, matter);
    return (structComp.distance_difference) > perAtomNorm(matter);
  }
}

// bool Matter::operator!=(const Matter& matter) {
//     return !operator==(matter);
// }

// Returns the distance to the given matter object.
double Matter::distanceTo(const Matter &matter) {
  return pbc(positions - matter.positions).norm();
}

// Returns the maximum distance between two atoms in the Matter objects.
double Matter::perAtomNorm(const Matter &matter) {
  long i = 0;
  double max_distance = 0.0;

  if (matter.numberOfAtoms() == nAtoms) {
    AtomMatrix diff = pbc(positions - matter.positions);
    for (i = 0; i < nAtoms; i++) {
      max_distance = std::max(diff.row(i).norm(), max_distance);
    }
  }
  return max_distance;
}

void Matter::resize(const long int length) {
  if (length > 0) {
    nAtoms = length;
    positions.resize(length, 3);
    positions.setZero();

    velocities.resize(length, 3);
    velocities.setZero();

    biasForces.resize(length, 3);
    biasForces.setZero();

    forces.resize(length, 3);
    forces.setZero();

    masses.resize(length);
    masses.setZero();

    atomicNrs.resize(length);
    atomicNrs.setZero();

    isFixed.resize(length);
    isFixed.setZero();

    atomIndex.resize(length);
    for (long i = 0; i < length; i++)
      atomIndex(i) = static_cast<int>(i); // default: sequential
  }
  recomputePotential = true;
  recomputeMaskedForces = true;
  recomputeFreeMask = true;
}

long int Matter::numberOfAtoms() const { return (nAtoms); }

Matrix3d Matter::getCell() const { return cell; }

void Matter::setCell(const Matrix3d &newCell) {
  cell = newCell;
  cellInverse = cell.inverse();
}

double Matter::getPosition(long int indexAtom, int axis) const {
  return positions(indexAtom, axis);
}

void Matter::setPosition(long int indexAtom, int axis, double position) {
  positions(indexAtom, axis) = position;
  if (usePeriodicBoundaries) {
    applyPeriodicBoundary();
  }
  recomputePotential = true;
  recomputeMaskedForces = true;
}

void Matter::setVelocity(long int indexAtom, int axis, double vel) {
  velocities(indexAtom, axis) = vel;
}

// return coordinates of atoms by const reference (zero-copy)
const AtomMatrix &Matter::getPositions() const { return positions; }
// return a modifiable copy of positions
AtomMatrix Matter::getPositionsCopy() const { return positions; }

VectorXd Matter::getPositionsV() const {
  return VectorXd::Map(positions.data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getPositionsFree() const {
  getFree(); // ensure freeIndices is up to date
  AtomMatrix ret(static_cast<long>(freeIndices.size()), 3);
  for (size_t j = 0; j < freeIndices.size(); j++) {
    ret.row(static_cast<long>(j)) = positions.row(freeIndices[j]);
  }
  return ret;
}

VectorXi Matter::getAtomicNrsFree() const {
  return this->atomicNrs.array() * getFreeV().cast<int>().array();
}

bool Matter::relax(bool quiet, bool writeMovie, bool checkpoint,
                   std::string prefixMovie, std::string prefixCheckpoint) {
  return eonc::helpers::relaxMatter(*this, *parameters, quiet, writeMovie,
                                    checkpoint, prefixMovie, prefixCheckpoint);
}

VectorXd Matter::getPositionsFreeV() const {
  return VectorXd::Map(getPositionsFree().data(), 3 * numberOfFreeAtoms());
}

// update Matter with the new positions of the free atoms given in array 'pos'
void Matter::setPositions(const AtomMatrix &pos) {
  positions = pos;
  if (usePeriodicBoundaries) {
    applyPeriodicBoundary();
  }
  recomputePotential = true;
  recomputeMaskedForces = true;
}

// Same but takes vector instead of n x 3 matrix
void Matter::setPositionsV(const VectorXd &pos) {
  setPositions(AtomMatrix::Map(pos.data(), numberOfAtoms(), 3));
}

void Matter::setPositionsFree(const AtomMatrix &pos) {
  getFree(); // ensure freeIndices is up to date
  for (size_t j = 0; j < freeIndices.size(); j++) {
    positions.row(freeIndices[j]) = pos.row(static_cast<long>(j));
  }
  recomputePotential = true;
  recomputeMaskedForces = true;
}

void Matter::setPositionsFreeV(const VectorXd &pos) {
  setPositionsFree(AtomMatrix::Map(pos.data(), numberOfFreeAtoms(), 3));
}

AtomMatrix Matter::getBiasForces() {
  if (biasPotential != nullptr) {
    biasPotential->boost();
  }
  return biasForces.array() * getFree().array();
}

void Matter::setBiasPotential(BondBoost *bondBoost) {
  biasPotential = bondBoost;
}

void Matter::setBiasForces(const AtomMatrix &bf) {
  biasForces = bf.array() * getFree().array();
}
// Return forces with fixed atoms zeroed (cached).
// Note: not thread-safe. Concurrent reads on the same Matter instance
// may race on the mutable maskedForces/recomputeMaskedForces members.
const AtomMatrix &Matter::getForces() const {
  computePotential();
  if (recomputeMaskedForces) {
    // Use the cached freeMask (Nx3, 1.0 for free / 0.0 for fixed) to zero
    // fixed-atom forces in a single vectorized Eigen operation.
    maskedForces = forces.array() * getFree().array();
    recomputeMaskedForces = false;
  }
  return maskedForces;
}

const AtomMatrix &Matter::getForcesRaw() const {
  computePotential();
  return forces;
}

VectorXd Matter::getForcesV() const {
  return VectorXd::Map(getForces().data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getForcesFree() const {
  AtomMatrix allForces = getForces();
  getFree(); // ensure freeIndices is up to date (mutable cache)
  AtomMatrix ret(static_cast<long>(freeIndices.size()), 3);
  for (size_t j = 0; j < freeIndices.size(); j++) {
    ret.row(static_cast<long>(j)) = allForces.row(freeIndices[j]);
  }
  return ret;
}

VectorXd Matter::getForcesFreeV() const {
  AtomMatrix freeForces = getForcesFree();
  return VectorXd::Map(freeForces.data(), 3 * numberOfFreeAtoms());
}

// return distance between the atoms with index1 and index2
double Matter::distance(long index1, long index2) const {
  return pbc(positions.row(index1) - positions.row(index2)).norm();
}

// return projected distance between the atoms with index1 and index2 on asix
// (0-x,1-y,2-z)
double Matter::pdistance(long index1, long index2, int axis) const {
  Matrix<double, 1, 3> ret;
  ret.setZero();
  ret(0, axis) = positions(index1, axis) - positions(index2, axis);
  ret = pbc(ret);
  return ret(0, axis);
}

// return the distance atom with index has moved between the current Matter
// object and the Matter object passed as argument
double Matter::distance(const Matter &matter, long index) const {
  return pbc(positions.row(index) - matter.getPositions().row(index)).norm();
}

double Matter::getMass(long int indexAtom) const { return (masses[indexAtom]); }

void Matter::setMass(long int indexAtom, double mass) {
  masses[indexAtom] = mass;
}

void Matter::setMasses(const VectorXd &massesIn) {
  for (int i = 0; i < nAtoms; i++) {
    masses[i] = massesIn[i];
  }
}

long Matter::getAtomicNr(long int indexAtom) const {
  return (atomicNrs[indexAtom]);
}

void Matter::setAtomicNr(long int indexAtom, long atomicNr) {
  atomicNrs[indexAtom] = atomicNr;
  recomputePotential = true;
  recomputeMaskedForces = true;
}

int Matter::getFixed(long int indexAtom) const { return (isFixed[indexAtom]); }

void Matter::setFixed(long int indexAtom, int isFixed_passed) {
  isFixed[indexAtom] = isFixed_passed;
  recomputeFreeMask = true;
  recomputeMaskedForces = true;
}

// void Matter::setPotentialEnergy(double epot_input)
//{
//	potentialEnergy=epot_input;
// }

double Matter::getPotentialEnergy() const {
  if (nAtoms > 0) {
    computePotential();
    return potentialEnergy;
  } else
    return 0.0;
}

double Matter::getKineticEnergy() const {
  // Vectorized: 0.5 * sum(mass_i * |v_i|^2) for free atoms
  // Use free mask (1-isFixed) to zero out fixed atom contributions
  Eigen::VectorXd speed2(nAtoms);
  for (long i = 0; i < nAtoms; i++) {
    speed2[i] = velocities.row(i).squaredNorm();
  }
  auto freeMask = (1 - isFixed.array()).cast<double>();
  return 0.5 * (masses.array() * freeMask * speed2.array()).sum();
}

double Matter::getMechanicalEnergy() const {
  return getPotentialEnergy() + getKineticEnergy();
}

long int Matter::numberOfFreeAtoms() const { return nAtoms - isFixed.sum(); }

long int Matter::numberOfFixedAtoms() const { return isFixed.sum(); }

long Matter::getForceCalls() const { return (forceCalls); }

void Matter::resetForceCalls() {
  forceCalls = 0;
  return;
}

void Matter::computePotential() const {
  if (recomputePotential) {
    if (!potential) {
      throw std::runtime_error(
          "Matter::computePotential called without a potential");
    }
    if (potential->isSurrogate()) {
      // Surrogate potential case: uses free-atom subset interface
      auto surrogatePotential =
          static_cast<SurrogatePotential *>(potential.get());
      auto [freePE, freeForces, vari] = surrogatePotential->get_ef_var(
          this->getPositionsFree(), this->getAtomicNrsFree(), cell);
      this->potentialEnergy = freePE;
      this->energyVariance = vari;
      for (long idx{0}, jdx{0}; idx < nAtoms; idx++) {
        if (!isFixed(idx)) {
          forces.row(idx) = freeForces.row(jdx);
          jdx++;
        }
      }
    } else {
      // Hot path: call force() directly into member storage.
      // No intermediate allocation, no tuple, no copy.
      double var{0};
      potential->force(nAtoms, positions.data(), atomicNrs.data(),
                       forces.data(), &potentialEnergy, &var, cell.data());
      potential->forceCallCounter++;
      PotRegistry::get().on_force_call(potential->getType());
    }
    forceCalls = forceCalls + 1;
    recomputePotential = false;

    if (isFixed.sum() == 0 && removeNetForce) {
      Vector3d tempForce = forces.colwise().sum() / nAtoms;
      for (long int i = 0; i < nAtoms; i++) {
        forces.row(i) -= tempForce.transpose();
      }
    }
  }
}

// Transform coordinates into the cell using the selected PBC convention (#176).
// Legacy: fractional [0,1) via fmod (historical). MinimumImage: fractional
// [-0.5,0.5) via floor (same MIC as eonc::pbc::apply for differences).
void Matter::applyPeriodicBoundary() {
  positions =
      eonc::pbc::applyPositions(positions, cell, cellInverse, pbcConvention);
}

double Matter::maxForce() const {
  // Ensures that the forces are up to date
  computePotential();

  // I think this can be done in one line with the rowwise method
  double maxForce = 0.0;
  for (int i = 0; i < nAtoms; i++) {
    if (getFixed(i)) {
      continue;
    }
    maxForce = std::max(forces.row(i).norm(), maxForce);
  }
  return maxForce;
}

VectorXi Matter::getAtomicNrs() const { return this->atomicNrs; }

void Matter::setAtomicNrs(const VectorXi &atmnrs) {
  if (atmnrs.size() != this->nAtoms) {
    throw std::invalid_argument(
        "Vector of atomic numbers not equal to the number of atoms");
  } else {
    this->atomicNrs = atmnrs;
  }
}

AtomMatrix Matter::getFree() const {
  if (recomputeFreeMask) {
    freeMask.resize(nAtoms, 3);
    freeIndices.clear();
    freeIndices.reserve(nAtoms);
    for (long i = 0; i < nAtoms; i++) {
      double val = isFixed(i) ? 0.0 : 1.0;
      freeMask(i, 0) = val;
      freeMask(i, 1) = val;
      freeMask(i, 2) = val;
      if (!isFixed(i))
        freeIndices.push_back(static_cast<int>(i));
    }
    recomputeFreeMask = false;
  }
  return freeMask;
}

VectorXd Matter::getFreeV() const {
  return VectorXd::Map(getFree().data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getVelocities() const {
  return velocities.array() * getFree().array();
}

void Matter::setVelocities(const AtomMatrix &v) {
  velocities = v.array() * getFree().array();
}

void Matter::setForces(const AtomMatrix &f) {
  forces = f.array() * getFree().array();
}

AtomMatrix Matter::getAccelerations() {
  AtomMatrix totF = getForces() + getBiasForces();
  AtomMatrix ret = totF.array() * getFree().array();
  // Single reciprocal mass computation, replicated across 3 columns
  auto invMass = masses.array().inverse();
  ret.col(0).array() *= invMass;
  ret.col(1).array() *= invMass;
  ret.col(2).array() *= invMass;
  return ret;
}

Matrix<double, Eigen::Dynamic, 1> Matter::getMasses() const { return masses; }

void Matter::setPotential(std::shared_ptr<Potential> pot) {
  this->potential = pot;
  recomputePotential = true;
  recomputeMaskedForces = true;
}

void Matter::setComputedPotential(double energy, double variance) {
  potentialEnergy = energy;
  energyVariance = variance;
  recomputePotential = false;
  recomputeMaskedForces = true;
  forceCalls++;

  // Apply the same net force removal as computePotential()
  if (isFixed.sum() == 0 && removeNetForce) {
    Vector3d tempForce = forces.colwise().sum() / nAtoms;
    for (long int i = 0; i < nAtoms; i++) {
      forces.row(i) -= tempForce.transpose();
    }
  }
}

size_t Matter::getPotentialCalls() const {
  return this->potential->forceCallCounter;
}

double Matter::getEnergyVariance() const { return this->energyVariance; }

// Eigen::VectorXd Matter::getForceVariance() {
//   return this->variance.segment(1, numberOfFreeAtoms() * 3);
// }

// double Matter::getMaxVariance() { return this->variance.maxCoeff(); }

std::shared_ptr<Potential> Matter::getPotential() { return this->potential; }

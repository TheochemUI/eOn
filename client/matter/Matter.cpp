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
#include "client/matter/Matter.h"
#include <cachelot/hash_fnv1a.h>
#include <stdexcept>
// #include "BondBoost.h"

namespace eonc {

Matter::Matter(const Matter &matter) { operator=(matter); }

const Matter &Matter::operator=(const Matter &matter) {
  nAtoms = matter.nAtoms;
  resize(nAtoms);

  positions = matter.positions;
  masses = matter.masses;
  atomicNrs = matter.atomicNrs;
  isFixed = matter.isFixed;
  cell = matter.cell;
  cellInverse = matter.cellInverse;
  velocities = matter.velocities;
  atmindices = matter.atmindices;

  removeNetForce = matter.removeNetForce;
  usePeriodicBoundaries = matter.usePeriodicBoundaries;

  potential = matter.potential;

  return *this;
}

// The == comparison considers identity. This is crucial for process search.
// bool Matter::operator==(const Matter& matter) {
//     if(parameters->checkRotation) {
//         return helper_functions::rotationMatch(this, &matter,
//         parameters->distanceDifference);
//     }else{
//         return (parameters->distanceDifference) > perAtomNorm(matter);
//     }
// }

// bool Matter::operator!=(const Matter& matter) {
//     return !operator==(matter);
// }

// Returns the distance to the given matter object.
double Matter::distanceTo(const Matter &matter) {
  return pbc(positions - matter.positions).norm();
}

AtomMatrix Matter::pbc(AtomMatrix diff) const {
  AtomMatrix ddiff = diff * cellInverse;

  int i, j;
  for (i = 0; i < diff.rows(); i++) {
    for (j = 0; j < 3; j++) {
      ddiff(i, j) = fmod(fmod(ddiff(i, j), 1.0) + 1.5, 1.0) - .5;
    }
  }

  return ddiff * cell;
}

VectorType Matter::pbcV(VectorType diffVector) const {
  AtomMatrix pbcMatrix =
      pbc(AtomMatrix::Map(diffVector.data(), diffVector.size() / 3, 3));
  return VectorType::Map(pbcMatrix.data(), diffVector.size());
}

// Returns the maximum distance between two atoms in the Matter objects.
double Matter::perAtomNorm(const Matter &matter) const {
  size_t i{0};
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

    masses.resize(length);
    masses.setZero();

    atomicNrs.resize(length);
    atomicNrs.setZero();

    isFixed.resize(length);
    isFixed.setZero();

    atmindices.resize(length);
    atmindices.setZero();
  }
}

size_t Matter::numberOfAtoms() const { return (nAtoms); }

double Matter::getPosition(long int indexAtom, int axis) const {
  return positions(indexAtom, axis);
}

void Matter::setPosition(long int indexAtom, int axis, double position) {
  positions(indexAtom, axis) = position;
  if (usePeriodicBoundaries) {
    applyPeriodicBoundary();
  }
}

void Matter::setVelocity(long int indexAtom, int axis, double vel) {
  velocities(indexAtom, axis) = vel;
}

// return coordinates of atoms in array 'pos'
AtomMatrix Matter::getPositions() const { return positions; }

VectorType Matter::getPositionsV() const {
  return VectorType::Map(positions.data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getPositionsFree() const {
  AtomMatrix ret(numberOfFreeAtoms(), 3);
  size_t i{0}, j{0};
  for (i = 0; i < nAtoms; i++) {
    if (!isFixed(i)) {
      ret.row(j) = positions.row(i);
      j++;
    }
  }
  return ret;
}

Vector<size_t> Matter::getAtomicNrsFree() const {
  return this->atomicNrs.array() * getFreeV().cast<size_t>().array();
}

// bool Matter::relax(bool quiet, bool writeMovie, bool checkpoint,
//                    string prefixMovie, string prefixCheckpoint) {
//   auto objf =
//       std::make_shared<MatterObjectiveFunction>(Matter(*this), parameters);
//   auto optim =
//       helpers::create::mkOptim(objf, parameters->optim.method, parameters);

//   ostringstream min;
//   min << prefixMovie;
//   if (writeMovie) {
//     matter2con(min.str(), false);
//   }

//   int iteration = 0;
//   if (!quiet) {
//     SPDLOG_LOGGER_DEBUG(m_log, "{} {:10s}  {:14s}  {:18s}  {:13s}\n",
//                         "[Matter]", "Iter", "Step size",
//                         parameters->optim.convergenceMetricLabel, "Energy");
//     SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}\n",
//                         "[Matter]", iteration, 0.0, objf->getConvergence(),
//                         getPotentialEnergy());
//   }

//   while (!objf->isConverged() && iteration < parameters->optim.maxIterations)
//   {

//     AtomMatrix pos = getPositions();

//     optim->step(parameters->optim.maxMove);
//     iteration++;
//     setPositionsFreeV(objf->getPositions());

//     double stepSize =
//         helper_functions::maxAtomMotion(pbc(getPositions() - pos));

//     if (!quiet) {
//       SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
//                           "[Matter]", iteration, stepSize,
//                           objf->getConvergence(), getPotentialEnergy());
//     }

//     if (writeMovie) {
//       matter2con(min.str(), true);
//     }

//     if (checkpoint) {
//       ostringstream chk;
//       chk << prefixCheckpoint << "_cp";
//       matter2con(chk.str(), false);
//     }
//   }

//   if (iteration == 0) {
//     if (!quiet) {
//       SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
//                           "[Matter]", iteration, 0.0, objf->getConvergence(),
//                           getPotentialEnergy());
//     }
//   }
//   //    bool converged = optimizer->run(parameters->optMaxIterations,
//   //    parameters->optMaxMove);
//   return objf->isConverged();
// }

VectorType Matter::getPositionsFreeV() const {
  return VectorType::Map(getPositionsFree().data(), 3 * numberOfFreeAtoms());
}

// update Matter with the new positions of the free atoms given in array 'pos'
void Matter::setPositions(const AtomMatrix pos) {
  positions = pos;
  if (usePeriodicBoundaries) {
    applyPeriodicBoundary();
  }
}

// Same but takes vector instead of n x 3 matrix
void Matter::setPositionsV(const VectorType pos) {
  setPositions(AtomMatrix::Map(pos.data(), numberOfAtoms(), 3));
}

void Matter::setPositionsFree(const AtomMatrix pos) {
  // FIXME: Ensure pos and existing data are in the same form with atom ids
  size_t i{0}, j{0};
  for (i = 0; i < nAtoms; i++) {
    if (!isFixed(i)) {
      positions.row(i) = pos.row(j);
      j += 1;
    }
  }
}

void Matter::setPositionsFreeV(const VectorType pos) {
  setPositionsFree(AtomMatrix::Map(pos.data(), numberOfFreeAtoms(), 3));
}

// AtomMatrix Matter::getBiasForces() {
//   if (biasPotential != NULL) {
//     biasPotential->boost();
//   }
//   return biasForces.array() * getFree().array();
// }

// void Matter::setBiasForces(const AtomMatrix bf) {
//   biasForces = bf.array() * getFree().array();
// }

// void Matter::setBiasPotential(BondBoost *bondBoost) {
//   biasPotential = bondBoost;
// }

AtomMatrix Matter::getForces() const {
  if (!potential) {
    throw(std::runtime_error("Whoops, you need a potential.."));
  }

  auto &&[_, frcs] = potential->get_ef(positions, atomicNrs, cell);

  if (isFixed.sum() == 0 && /* params->removeNetForce */ false) {
    Eigen::Vector3d tempForce(3);
    tempForce = frcs.colwise().sum() / nAtoms;

    for (long int i = 0; i < nAtoms; i++) {
      frcs.row(i) -= tempForce.transpose();
    }
  }

  size_t i{0};
  for (i = 0; i < nAtoms; i++) {
    if (isFixed[i]) {
      frcs.row(i).setZero();
    }
  }
  return frcs;
}

VectorType Matter::getForcesV() const {
  return VectorType::Map(getForces().data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getForcesFree() const {
  AtomMatrix allForces = getForces();
  AtomMatrix ret(numberOfFreeAtoms(), 3);
  size_t i{0}, j{0};
  for (i = 0; i < nAtoms; i++) {
    if (!isFixed[i]) {
      ret.row(j) = allForces.row(i);
      j++;
    }
  }
  return ret;
}

VectorType Matter::getForcesFreeV() const {
  return VectorType::Map(getForcesFree().data(), 3 * numberOfFreeAtoms());
}

// return distance between the atoms with index1 and index2
double Matter::distance(long index1, long index2) const {
  return pbc(positions.row(index1) - positions.row(index2)).norm();
}

// return projected distance between the atoms with index1 and index2 on axis
// (0-x,1-y,2-z)
double Matter::pdistance(long index1, long index2, int axis) const {
  MatrixTRC<double, 1, 3> ret;
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

void Matter::setMasses(VectorType massesIn) {
  for (size_t i = 0; i < nAtoms; i++) {
    masses[i] = massesIn[i];
  }
}

size_t Matter::getAtomicNr(size_t indexAtom) const {
  return (atomicNrs[indexAtom]);
}

void Matter::setAtomicNr(long int indexAtom, long atomicNr) {
  atomicNrs[indexAtom] = atomicNr;
}

int Matter::getFixed(long int indexAtom) const { return (isFixed[indexAtom]); }

void Matter::setFixed(long int indexAtom, int isFixed_passed) {
  isFixed[indexAtom] = isFixed_passed;
}

void Matter::setFixedMask(const Vector<int> &fixmask) { isFixed = fixmask; }
// void Matter::setPotentialEnergy(double epot_input)
//{
//	potentialEnergy=epot_input;
// }

double Matter::getPotentialEnergy() const {
#ifdef EON_CHECKS
  if (nAtoms <= 0) {
    throw std::runtime_error("Called for a potential with 0 atoms");
  }
#endif
  if (!potential) {
    throw(std::runtime_error("Whoops, you need a potential.."));
  }
  auto &&[potE, _] = potential->get_ef(positions, atomicNrs, cell);
  return potE;
}

double Matter::getKineticEnergy() const {
  double K = 0;
  for (size_t i = 0; i < nAtoms; i++) {
    if (!isFixed[i])
      K += masses[i] * 0.5 * velocities.row(i).squaredNorm();
  }
  return K;
}

double Matter::getMechanicalEnergy() {
  return getPotentialEnergy() + getKineticEnergy();
}

long int Matter::numberOfFreeAtoms() const { return nAtoms - isFixed.sum(); }

long int Matter::numberOfFixedAtoms() const { return isFixed.sum(); }

size_t Matter::getForceCalls() const { return this->forceCalls; }

void Matter::resetForceCalls() {
  forceCalls = 0;
  return;
}

// Transform the coordinate to use the minimum image convention.
void Matter::applyPeriodicBoundary() {
  AtomMatrix ddiff = positions * cellInverse;

  int i, j;
  for (i = 0; i < ddiff.rows(); i++) {
    for (j = 0; j < 3; j++) {
      ddiff(i, j) = fmod(ddiff(i, j) + 1.0, 1.0);
    }
  }
  positions = ddiff * cell;
}

double Matter::maxForce(void) {
  // Ensures that the forces are up to date
  return (getForcesFree()).rowwise().norm().maxCoeff();
}

Vector<size_t> Matter::getAtomicNrs() const { return this->atomicNrs; }

void Matter::setAtomicNrs(const Vector<size_t> atmnrs) {
  if (static_cast<size_t>(atmnrs.size()) != this->nAtoms) {
    throw std::invalid_argument(
        "Vector of atomic numbers not equal to the number of atoms");
  } else {
    this->atomicNrs = atmnrs;
  }
}

AtomMatrix Matter::getFree() const {
  AtomMatrix ret(nAtoms, 3);
  size_t i{0}, j{0};
  for (i = 0; i < nAtoms; i++) {
    for (j = 0; j < 3; j++) {
      ret(i, j) = double(!bool(isFixed(i)));
    }
  }
  return ret;
}

VectorType Matter::getFreeV() const {
  return VectorType::Map(getFree().data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getVelocities() const {
  return velocities.array() * getFree().array();
}

void Matter::setVelocities(const AtomMatrix v) {
  velocities = v.array() * getFree().array();
}

// AtomMatrix Matter::getAccelerations() {
//   AtomMatrix totF = getForces() + getBiasForces();
//   AtomMatrix ret = totF.array() * getFree().array();
//   ret.col(0).array() /= masses.array();
//   ret.col(1).array() /= masses.array();
//   ret.col(2).array() /= masses.array();
//   return ret;
// }

VectorType Matter::getMasses() const { return masses; }

void Matter::setPotential(std::shared_ptr<PotBase> pot) {
  this->potential = pot;
}

size_t Matter::getPotentialCalls() const { return this->npotcalls; }

double Matter::getEnergyVariance() { return this->energyVariance; }

// VectorType Matter::getForceVariance() {
//   return this->variance.segment(1, numberOfFreeAtoms() * 3);
// }

// double Matter::getMaxVariance() { return this->variance.maxCoeff(); }

std::shared_ptr<PotBase> Matter::getPotential() { return this->potential; }

} // namespace eonc

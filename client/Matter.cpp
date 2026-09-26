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
#include "eon/Matter.h"
#include "eon/BaseStructures.h"
#include "eon/BondBoost.h"
#include "eon/ForceNorm.h"
#include "eon/GeometryAnalysis.h"
#include "eon/HelperFunctions.h"
#include "eon/Parameters.h"
#include "eon/SurrogatePotential.h"

#include "eon/EonLogger.h"
#include <cmath>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>

namespace eonc {

struct Matter::Impl {
  AtomMatrix positions{AtomMatrix::Zero(0, 3)};
  AtomMatrix velocities{AtomMatrix::Zero(0, 3)};
  mutable AtomMatrix forces{AtomMatrix::Zero(0, 3)};
  AtomMatrix biasForces{AtomMatrix::Zero(0, 3)};
  VectorXd masses{VectorXd::Zero(0)};
  VectorXi atomicNrs{VectorXi::Zero(0)};
  // Nx3; 1.0 if that axis is fixed, 0.0 if free.
  AtomMatrix isFixed{AtomMatrix::Zero(0, 3)};
  // Original atom index from .con column 5.
  Eigen::Matrix<std::int64_t, Eigen::Dynamic, 1> atomIndex;
  mutable AtomMatrix freeMask;     // Nx3; 1.0 free, 0.0 fixed
  mutable AtomMatrix maskedForces; // forces with fixed atoms zeroed
  Matrix3d cell{Matrix3d::Zero()};
  Matrix3d cellInverse{Matrix3d::Zero()};
};

Matter::~Matter() = default;

Matter::Matter(std::shared_ptr<Potential> pot, const Parameters &params)
    : potential{pot},
      usePeriodicBoundaries{!(pot && pot->requiresIsolatedMoleculeLayout())},
      pbcConvention{PbcConvention::Legacy},
      recomputePotential{true},
      forceCalls{0},
      removeNetForce{params.main_options().removeNetForce},
      structComp{params.structure_comparison_options()},
      parameters{&params},
      nAtoms{0},
      impl_{std::make_unique<Impl>()},
      biasPotential{nullptr},
      energyVariance{0.0},
      potentialEnergy{0.0} {}

bool Matter::getWriteConForces() const noexcept {
  return parameters != nullptr && parameters->main_options().writeConForces;
}

namespace {
void checkAtom(long nAtoms, long indexAtom, const char *fn) {
  if (indexAtom < 0 || indexAtom >= nAtoms) {
    throw std::out_of_range(std::string(fn) + ": atom index out of range");
  }
}
void checkAxis(int axis, const char *fn) {
  if (axis < 0 || axis > 2) {
    throw std::out_of_range(std::string(fn) + ": axis out of range");
  }
}
} // namespace

Matter::Matter(const Matter &matter)
    : impl_{std::make_unique<Impl>()} {
  operator=(matter);
}

const Matter &Matter::operator=(const Matter &matter) {
  if (this == &matter) {
    return *this;
  }
  nAtoms = matter.nAtoms;
  resize(nAtoms);

  impl_->positions = matter.impl_->positions;
  impl_->forces = matter.impl_->forces;
  impl_->masses = matter.impl_->masses;
  impl_->atomicNrs = matter.impl_->atomicNrs;
  impl_->isFixed = matter.impl_->isFixed;
  impl_->atomIndex = matter.impl_->atomIndex;
  fileToMatter = matter.fileToMatter;
  impl_->cell = matter.impl_->cell;
  impl_->cellInverse = matter.impl_->cellInverse;
  impl_->velocities = matter.impl_->velocities;

  removeNetForce = matter.removeNetForce;
  structComp = matter.structComp;
  parameters = matter.parameters;

  usePeriodicBoundaries = matter.usePeriodicBoundaries;
  pbcConvention = matter.pbcConvention;

  potential = matter.potential;
  potentialEnergy = matter.potentialEnergy;
  energyVariance = matter.energyVariance;
  forceCalls = matter.forceCalls;
  recomputePotential = matter.recomputePotential;
  // Both caches describe the forces this object held before the assignment.
  // resize() above already raises them; state it here alongside the members
  // this function owns.
  recomputeFreeMask = true;
  recomputeMaskedForces = true;

  // A BondBoost binds to one Matter (BondBoost.h:32 takes a Matter *), so a
  // copy cannot share the source's: boosting through it would drive the
  // original. The copy starts without one and re-establishes it through
  // setBiasPotential. biasForces is zeroed by the resize above, which is the
  // state that matches having no bias potential.
  biasPotential = nullptr;

  headerCon = matter.headerCon;
  // ConFrame is move-only; copy does not retain movie trajectory.
  movie_frames_.clear();

  return *this;
}

Matter::Matter(Matter &&other) noexcept
    : impl_{std::make_unique<Impl>()} {
  operator=(std::move(other));
}

Matter &Matter::operator=(Matter &&other) noexcept {
  if (this == &other) {
    return *this;
  }
  potential = std::move(other.potential);
  usePeriodicBoundaries = other.usePeriodicBoundaries;
  pbcConvention = other.pbcConvention;
  recomputePotential = other.recomputePotential;
  forceCalls = other.forceCalls;
  headerCon = std::move(other.headerCon);
  removeNetForce = other.removeNetForce;
  structComp = other.structComp;
  parameters = other.parameters;
  nAtoms = other.nAtoms;
  impl_->positions = std::move(other.impl_->positions);
  impl_->velocities = std::move(other.impl_->velocities);
  impl_->forces = std::move(other.impl_->forces);
  impl_->biasForces = std::move(other.impl_->biasForces);
  biasPotential = other.biasPotential;
  other.biasPotential = nullptr;
  impl_->masses = std::move(other.impl_->masses);
  impl_->atomicNrs = std::move(other.impl_->atomicNrs);
  impl_->isFixed = std::move(other.impl_->isFixed);
  impl_->atomIndex = std::move(other.impl_->atomIndex);
  fileToMatter = std::move(other.fileToMatter);
  impl_->freeMask = std::move(other.impl_->freeMask);
  impl_->maskedForces = std::move(other.impl_->maskedForces);
  freeIndices = std::move(other.freeIndices);
  recomputeFreeMask = other.recomputeFreeMask;
  recomputeMaskedForces = other.recomputeMaskedForces;
  impl_->cell = std::move(other.impl_->cell);
  impl_->cellInverse = std::move(other.impl_->cellInverse);
  energyVariance = other.energyVariance;
  movie_frames_ = std::move(other.movie_frames_);
  potentialEnergy = other.potentialEnergy;

  other.nAtoms = 0;
  other.recomputePotential = true;
  other.recomputeFreeMask = true;
  other.recomputeMaskedForces = true;
  return *this;
}

bool Matter::compare(const Matter &matter, bool indistinguishable) {
  if (nAtoms != matter.numberOfAtoms())
    return false;
  if (structComp.check_rotation && indistinguishable) {
    return eonc::geometry::sortedR(*this, matter,
                                   structComp.distance_difference);
  } else if (indistinguishable) {
    if (this->numberOfFixedAtoms() == 0 and structComp.remove_translation)
      eonc::geometry::translationRemove(*this, matter);
    return eonc::geometry::identical(*this, matter,
                                     structComp.distance_difference);
  } else if (structComp.check_rotation) {
    return eonc::geometry::rotationMatch(*this, matter,
                                         structComp.distance_difference);
  } else {
    if (this->numberOfFixedAtoms() == 0 and structComp.remove_translation)
      eonc::geometry::translationRemove(*this, matter);
    return (structComp.distance_difference) > perAtomNorm(matter);
  }
}

// Returns the distance to the given matter object.
double Matter::distanceTo(const Matter &matter) {
  if (matter.numberOfAtoms() != nAtoms) {
    throw std::invalid_argument("Matter::distanceTo: size mismatch");
  }
  return pbc(impl_->positions - matter.impl_->positions).norm();
}

// Returns the maximum distance between two atoms in the Matter objects.
double Matter::perAtomNorm(const Matter &matter) {
  long i = 0;
  double max_distance = 0.0;

  if (matter.numberOfAtoms() == nAtoms) {
    AtomMatrix diff = pbc(impl_->positions - matter.impl_->positions);
    for (i = 0; i < nAtoms; i++) {
      max_distance = std::max(diff.row(i).norm(), max_distance);
    }
  }
  return max_distance;
}

void Matter::resize(const long int length) {
  if (length < 0) {
    throw std::invalid_argument("Matter::resize: negative atom count");
  }
  // Same-N resize still zeros coordinates. Keep .con column-5 ids and
  // the file-order map so a later matter2con does not stamp 1..N.
  const bool keepAtomIds =
      (length == nAtoms && impl_->atomIndex.size() == length &&
       fileToMatter.size() == static_cast<size_t>(length));
  // Zero is a real size: leaving nAtoms at the old value there sends
  // setMasses and every other nAtoms loop off the end of an empty array.
  nAtoms = length;
  impl_->positions.resize(length, 3);
  impl_->positions.setZero();

  impl_->velocities.resize(length, 3);
  impl_->velocities.setZero();

  impl_->biasForces.resize(length, 3);
  impl_->biasForces.setZero();

  impl_->forces.resize(length, 3);
  impl_->forces.setZero();

  impl_->masses.resize(length);
  impl_->masses.setZero();

  impl_->atomicNrs.resize(length);
  impl_->atomicNrs.setZero();

  impl_->isFixed.resize(length, 3);
  impl_->isFixed.setZero();

  if (!keepAtomIds) {
    impl_->atomIndex.resize(length);
    fileToMatter.resize(static_cast<size_t>(length));
    for (long i = 0; i < length; i++) {
      impl_->atomIndex(i) = static_cast<std::int64_t>(i); // default: sequential
      fileToMatter[static_cast<size_t>(i)] = i;
    }
  }
  recomputePotential = true;
  recomputeMaskedForces = true;
  recomputeFreeMask = true;
}

long int Matter::numberOfAtoms() const { return (nAtoms); }

Matrix3d Matter::getCell() const { return impl_->cell; }

void Matter::setCell(const Matrix3d &newCell) {
  impl_->cell = newCell;
  impl_->cellInverse = impl_->cell.inverse();
  recomputePotential = true;
  recomputeMaskedForces = true;
}

double Matter::getPosition(long int indexAtom, int axis) const {
  checkAtom(nAtoms, indexAtom, "Matter::getPosition");
  checkAxis(axis, "Matter::getPosition");
  return impl_->positions(indexAtom, axis);
}

void Matter::setPosition(long int indexAtom, int axis, double position) {
  checkAtom(nAtoms, indexAtom, "Matter::setPosition");
  checkAxis(axis, "Matter::setPosition");
  impl_->positions(indexAtom, axis) = position;
  if (usePeriodicBoundaries) {
    applyPeriodicBoundary();
  }
  recomputePotential = true;
  recomputeMaskedForces = true;
}

void Matter::setVelocity(long int indexAtom, int axis, double vel) {
  checkAtom(nAtoms, indexAtom, "Matter::setVelocity");
  checkAxis(axis, "Matter::setVelocity");
  impl_->velocities(indexAtom, axis) = vel;
}

// return coordinates of atoms by const reference (zero-copy)
const AtomMatrix &Matter::getPositions() const { return impl_->positions; }
// return a modifiable copy of positions
AtomMatrix Matter::getPositionsCopy() const { return impl_->positions; }

VectorXd Matter::getPositionsV() const {
  return VectorXd::Map(impl_->positions.data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getPositionsFree() const {
  getFree(); // ensure freeIndices is up to date
  AtomMatrix ret(static_cast<long>(freeIndices.size()), 3);
  for (size_t j = 0; j < freeIndices.size(); j++) {
    ret.row(static_cast<long>(j)) = impl_->positions.row(freeIndices[j]);
  }
  return ret;
}

VectorXi Matter::getAtomicNrsFree() const {
  getFree();
  VectorXi ret(static_cast<Eigen::Index>(freeIndices.size()));
  for (size_t j = 0; j < freeIndices.size(); j++) {
    ret[static_cast<Eigen::Index>(j)] = impl_->atomicNrs[freeIndices[j]];
  }
  return ret;
}

bool Matter::relax(bool quiet, bool writeMovie, bool checkpoint,
                   std::string prefixMovie, std::string prefixCheckpoint,
                   bool retainMovieFrames) {
  if (retainMovieFrames) {
    movie_frames_.clear();
  }
  return eonc::helpers::relaxMatter(
      *this, *parameters, quiet, writeMovie, checkpoint, prefixMovie,
      prefixCheckpoint, retainMovieFrames ? &movie_frames_ : nullptr);
}

VectorXd Matter::getPositionsFreeV() const {
  return VectorXd::Map(getPositionsFree().data(), 3 * numberOfFreeAtoms());
}

// update Matter with the new positions of the free atoms given in array 'pos'
void Matter::setPositions(const AtomMatrix &pos) {
  if (pos.rows() != nAtoms) {
    throw std::invalid_argument("Matter::setPositions: row count mismatch");
  }
  impl_->positions = pos;
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
    impl_->positions.row(freeIndices[j]) = pos.row(static_cast<long>(j));
  }
  // Optimizers write free-atom coords only (TIP4P/SPCE and any PBC pot).
  // Match setPositions: wrap the full configuration when PBC is on (#171).
  if (usePeriodicBoundaries) {
    applyPeriodicBoundary();
  }
  recomputePotential = true;
  recomputeMaskedForces = true;
}

void Matter::setPositionsFreeV(const VectorXd &pos) {
  setPositionsFree(AtomMatrix::Map(pos.data(), numberOfFreeAtoms(), 3));
}

AtomMatrix Matter::getBiasForces() {
  if (biasPotential != nullptr) {
    // Evaluate the current bias only. The MD job advances the bond-boost
    // schedule once per step via BondBoost::advance().
    biasPotential->boost();
  }
  return impl_->biasForces.array() * getFree().array();
}

void Matter::setBiasPotential(BondBoost *bondBoost) {
  biasPotential = bondBoost;
}

void Matter::setBiasForces(const AtomMatrix &bf) {
  impl_->biasForces = bf.array() * getFree().array();
}
// Return forces with fixed atoms zeroed (cached).
// Note: not thread-safe. Concurrent reads on the same Matter instance
// may race on the mutable maskedForces/recomputeMaskedForces members.
const AtomMatrix &Matter::getForces() const {
  computePotential();
  if (recomputeMaskedForces) {
    // Use the cached freeMask (Nx3, 1.0 for free / 0.0 for fixed) to zero
    // fixed-atom forces in a single vectorized Eigen operation.
    impl_->maskedForces = impl_->forces.array() * getFree().array();
    recomputeMaskedForces = false;
  }
  return impl_->maskedForces;
}

const AtomMatrix &Matter::getForcesRaw() const {
  computePotential();
  return impl_->forces;
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
  checkAtom(nAtoms, index1, "Matter::distance");
  checkAtom(nAtoms, index2, "Matter::distance");
  return pbc(impl_->positions.row(index1) - impl_->positions.row(index2))
      .norm();
}

// return projected distance between the atoms with index1 and index2 on asix
// (0-x,1-y,2-z)
double Matter::pdistance(long index1, long index2, int axis) const {
  checkAtom(nAtoms, index1, "Matter::pdistance");
  checkAtom(nAtoms, index2, "Matter::pdistance");
  checkAxis(axis, "Matter::pdistance");
  Matrix<double, 1, 3> ret;
  ret.setZero();
  ret(0, axis) =
      impl_->positions(index1, axis) - impl_->positions(index2, axis);
  ret = pbc(ret);
  return ret(0, axis);
}

// return the distance atom with index has moved between the current Matter
// object and the Matter object passed as argument
double Matter::distance(const Matter &matter, long index) const {
  checkAtom(nAtoms, index, "Matter::distance");
  checkAtom(matter.nAtoms, index, "Matter::distance");
  return pbc(impl_->positions.row(index) - matter.getPositions().row(index))
      .norm();
}

double Matter::getMass(long int indexAtom) const {
  checkAtom(nAtoms, indexAtom, "Matter::getMass");
  return (impl_->masses[indexAtom]);
}

void Matter::setMass(long int indexAtom, double mass) {
  checkAtom(nAtoms, indexAtom, "Matter::setMass");
  impl_->masses[indexAtom] = mass;
}

void Matter::setMasses(const VectorXd &massesIn) {
  if (massesIn.size() != nAtoms) {
    throw std::invalid_argument("Matter::setMasses: size mismatch");
  }
  impl_->masses = massesIn;
}

long Matter::getAtomicNr(long int indexAtom) const {
  return (impl_->atomicNrs[indexAtom]);
}

void Matter::setAtomicNr(long int indexAtom, long atomicNr) {
  impl_->atomicNrs[indexAtom] = atomicNr;
  recomputePotential = true;
  recomputeMaskedForces = true;
}

int Matter::getFixed(long int indexAtom) const {
  checkAtom(nAtoms, indexAtom, "Matter::getFixed");
  return (impl_->isFixed(indexAtom, 0) > 0.5 &&
          impl_->isFixed(indexAtom, 1) > 0.5 &&
          impl_->isFixed(indexAtom, 2) > 0.5)
             ? 1
             : 0;
}

int Matter::getFixed(long int indexAtom, int axis) const {
  checkAtom(nAtoms, indexAtom, "Matter::getFixed");
  checkAxis(axis, "Matter::getFixed");
  return impl_->isFixed(indexAtom, axis) > 0.5 ? 1 : 0;
}

std::array<bool, 3> Matter::getFixedMask(long int indexAtom) const {
  checkAtom(nAtoms, indexAtom, "Matter::getFixedMask");
  return {impl_->isFixed(indexAtom, 0) > 0.5,
          impl_->isFixed(indexAtom, 1) > 0.5,
          impl_->isFixed(indexAtom, 2) > 0.5};
}

void Matter::setFixed(long int indexAtom, int isFixed_passed) {
  checkAtom(nAtoms, indexAtom, "Matter::setFixed");
  const double v = isFixed_passed ? 1.0 : 0.0;
  impl_->isFixed(indexAtom, 0) = v;
  impl_->isFixed(indexAtom, 1) = v;
  impl_->isFixed(indexAtom, 2) = v;
  recomputeFreeMask = true;
  recomputeMaskedForces = true;
}

void Matter::setFixed(long int indexAtom, int axis, int isFixed_passed) {
  checkAtom(nAtoms, indexAtom, "Matter::setFixed");
  checkAxis(axis, "Matter::setFixed");
  impl_->isFixed(indexAtom, axis) = isFixed_passed ? 1.0 : 0.0;
  recomputeFreeMask = true;
  recomputeMaskedForces = true;
}

void Matter::setFixedMask(long int indexAtom, std::array<bool, 3> mask) {
  checkAtom(nAtoms, indexAtom, "Matter::setFixedMask");
  impl_->isFixed(indexAtom, 0) = mask[0] ? 1.0 : 0.0;
  impl_->isFixed(indexAtom, 1) = mask[1] ? 1.0 : 0.0;
  impl_->isFixed(indexAtom, 2) = mask[2] ? 1.0 : 0.0;
  recomputeFreeMask = true;
  recomputeMaskedForces = true;
}

double Matter::getPotentialEnergy() const {
  if (nAtoms > 0) {
    computePotential();
    return potentialEnergy;
  } else
    return 0.0;
}

double Matter::getKineticEnergy() const {
  // 0.5 * sum(mass_i * |v_i,free|^2); a constrained axis does not contribute.
  AtomMatrix vfree = impl_->velocities.array() * getFree().array();
  Eigen::VectorXd speed2 = vfree.rowwise().squaredNorm();
  return 0.5 * (impl_->masses.array() * speed2.array()).sum();
}

double Matter::getMechanicalEnergy() const {
  return getPotentialEnergy() + getKineticEnergy();
}

long int Matter::numberOfFixedAtoms() const {
  long n = 0;
  for (long i = 0; i < nAtoms; ++i) {
    if (getFixed(i)) {
      ++n;
    }
  }
  return n;
}

long int Matter::numberOfFreeAtoms() const {
  return nAtoms - numberOfFixedAtoms();
}

long Matter::getForceCalls() const { return (forceCalls); }

void Matter::resetForceCalls() {
  forceCalls = 0;
  return;
}

void Matter::assertIsolatedMoleculeLayoutSafe() const {
  if (!potential || !potential->requiresIsolatedMoleculeLayout()) {
    return;
  }
  if (usePeriodicBoundaries) {
    throw std::runtime_error(
        "Potential requires an isolated (non-periodic) molecular layout "
        "(NWChem/ORCA-class). Disable periodic boundaries before optimizing; "
        "PBC wraps can tear non-centered molecules (issue #188).");
  }
}

void Matter::computePotential() const {
  if (recomputePotential) {
    if (!potential) {
      throw std::runtime_error(
          "Matter::computePotential called without a potential");
    }
    assertIsolatedMoleculeLayoutSafe();
    if (potential->isSurrogate()) {
      // Surrogate potential case: uses free-atom subset interface
      auto surrogatePotential =
          static_cast<SurrogatePotential *>(potential.get());
      auto [freePE, freeForces, vari] = surrogatePotential->get_ef_var(
          this->getPositionsFree(), this->getAtomicNrsFree(), impl_->cell);
      this->potentialEnergy = freePE;
      this->energyVariance = vari;
      for (long idx{0}, jdx{0}; idx < nAtoms; idx++) {
        if (!getFixed(idx)) {
          impl_->forces.row(idx) = freeForces.row(jdx);
          jdx++;
        }
      }
    } else {
      // Hot path: call force() directly into member storage.
      // No intermediate allocation, no tuple, no copy.
      double var{0};
      potential->setFixedMask(nAtoms, impl_->isFixed.data());
      const auto n = static_cast<size_t>(nAtoms);
      // Isolated molecules still store a box for I/O. Pots that infer PBC
      // from a non-zero cell (GFN2) must see a zero box here.
      const Matrix3d force_cell =
          usePeriodicBoundaries ? impl_->cell : Matrix3d::Zero();
      potential->force(std::span<const double>(impl_->positions.data(), n * 3),
                       std::span<const int>(impl_->atomicNrs.data(), n),
                       std::span<double>(impl_->forces.data(), n * 3),
                       &potentialEnergy, &var,
                       std::span<const double>(force_cell.data(), 9));
      potential->forceCallCounter++;
      potential->notifyForceCall();
    }
    if (!std::isfinite(potentialEnergy) || !impl_->forces.allFinite()) {
      throw std::runtime_error(
          "Potential returned non-finite energy or forces");
    }
    forceCalls = forceCalls + 1;
    recomputePotential = false;

    // One free atom: subtracting the mean force is identically zero
    // (eOn-zjri). NEB would then report immediate GOOD.
    if (impl_->isFixed.maxCoeff() < 0.5 && removeNetForce && nAtoms > 1) {
      Vector3d tempForce = impl_->forces.colwise().sum() / nAtoms;
      for (long int i = 0; i < nAtoms; i++) {
        impl_->forces.row(i) -= tempForce.transpose();
      }
    }
  }
}

// Transform coordinates into the cell using the selected PBC convention (#176).
// Legacy: fractional [0,1) via fmod (historical). MinimumImage: fractional
// [-0.5,0.5) via floor (same MIC as eonc::pbc::apply for differences).
void Matter::applyPeriodicBoundary() {
  assertIsolatedMoleculeLayoutSafe();
  // Unset / singular cell (default after Matter construct is Zero) — do not
  // wipe coordinates; callers set the cell before wrapping makes sense.
  if (std::abs(impl_->cell.determinant()) < 1e-30) {
    return;
  }
  impl_->positions = eonc::pbc::applyPositions(
      impl_->positions, impl_->cell, impl_->cellInverse, pbcConvention);
}

double Matter::maxFreeAtomForce(const AtomMatrix &rows) const {
  if (rows.rows() != nAtoms) {
    throw std::invalid_argument(
        "Matter::maxFreeAtomForce: row count does not match atom count");
  }
  return maxFreeAtomForceNorm(rows.data(), isFixed.data(), nAtoms);
}

double Matter::maxForce() const {
  // Ensures that the forces are up to date
  computePotential();
  return maxFreeAtomForce(getForces());
}

VectorXi Matter::getAtomicNrs() const { return this->impl_->atomicNrs; }

void Matter::setAtomicNrs(const VectorXi &atmnrs) {
  if (atmnrs.size() != this->nAtoms) {
    throw std::invalid_argument(
        "Vector of atomic numbers not equal to the number of atoms");
  } else {
    this->impl_->atomicNrs = atmnrs;
    recomputePotential = true;
    recomputeMaskedForces = true;
  }
}

AtomMatrix Matter::getFree() const {
  if (recomputeFreeMask) {
    impl_->freeMask.resize(nAtoms, 3);
    impl_->freeMask = 1.0 - impl_->isFixed.array();
    freeIndices.clear();
    freeIndices.reserve(static_cast<size_t>(nAtoms));
    for (long i = 0; i < nAtoms; i++) {
      if (impl_->freeMask.row(i).sum() > 0.5) {
        freeIndices.push_back(static_cast<int>(i));
      }
    }
    recomputeFreeMask = false;
  }
  return impl_->freeMask;
}

VectorXd Matter::getFreeV() const {
  return VectorXd::Map(getFree().data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getVelocities() const {
  return impl_->velocities.array() * getFree().array();
}

void Matter::setVelocities(const AtomMatrix &v) {
  impl_->velocities = v.array() * getFree().array();
}

void Matter::setForces(const AtomMatrix &f) {
  impl_->forces = f.array() * getFree().array();
  impl_->maskedForces = impl_->forces;
  recomputeMaskedForces = false;
  recomputePotential = false;
}

AtomMatrix Matter::getAccelerations() {
  AtomMatrix totF = getForces() + getBiasForces();
  AtomMatrix ret = totF.array() * getFree().array();
  // Single reciprocal mass computation, replicated across 3 columns
  auto invMass = impl_->masses.array().inverse();
  ret.col(0).array() *= invMass;
  ret.col(1).array() *= invMass;
  ret.col(2).array() *= invMass;
  return ret;
}

Matrix<double, Eigen::Dynamic, 1> Matter::getMasses() const {
  return impl_->masses;
}

void Matter::setPotential(std::shared_ptr<Potential> pot) {
  this->potential = pot;
  // Molecular QM backends (NWChem/ORCA) must not use PBC wraps (#188). Auto-off
  // with a hard fail if code later forces PBC while this pot is attached.
  if (potential && potential->requiresIsolatedMoleculeLayout() &&
      usePeriodicBoundaries) {
    usePeriodicBoundaries = false;
    // Use EONC_LOG_* (not bare QUILL_LOG_* with eonc::log::get() as arg): the
    // quill macros expand logger->… and choke on a call expression as the first
    // macro argument on some GCC versions (CI: expected primary-expression).
    EONC_LOG_WARNING(
        "Disabled PBC for isolated-molecule potential (NWChem/ORCA-class); "
        "re-enabling PBC will throw (issue #188)");
  }
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
  if (impl_->isFixed.maxCoeff() < 0.5 && removeNetForce && nAtoms > 1) {
    Vector3d tempForce = impl_->forces.colwise().sum() / nAtoms;
    for (long int i = 0; i < nAtoms; i++) {
      impl_->forces.row(i) -= tempForce.transpose();
    }
  }
}

size_t Matter::getPotentialCalls() const {
  return this->potential->forceCallCounter;
}

double Matter::getEnergyVariance() const { return this->energyVariance; }

std::shared_ptr<Potential> Matter::getPotential() { return this->potential; }

AtomMatrix Matter::pbc(const AtomMatrix &diff) const {
  if (!usePeriodicBoundaries) {
    return diff;
  }
  return eonc::pbc::apply(diff, impl_->cell, impl_->cellInverse);
}

VectorXd Matter::pbcV(const VectorXd &diff) const {
  if (!usePeriodicBoundaries) {
    return diff;
  }
  return eonc::pbc::applyV(diff, impl_->cell, impl_->cellInverse);
}

double *Matter::forcesData() { return impl_->forces.data(); }

std::int64_t Matter::getAtomIndex(long int atom) const {
  return impl_->atomIndex(atom);
}

void Matter::setAtomIndex(long int atom, std::int64_t index) {
  impl_->atomIndex(atom) = index;
}

void Matter::restoreFileForces(const AtomMatrix &fileForces, bool trustEnergy,
                               double energy) {
  impl_->forces = fileForces;
  recomputeMaskedForces = true;
  if (trustEnergy) {
    potentialEnergy = energy;
    energyVariance = 0.0;
    recomputePotential = false;
  } else {
    recomputePotential = true;
  }
}

} // namespace eonc

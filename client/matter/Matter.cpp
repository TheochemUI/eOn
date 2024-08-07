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
#include "client/Element.hpp"
// #include "BondBoost.h"

namespace eonc {

Matter::Matter(const Matter &matter) { operator=(matter); }

const Matter &Matter::operator=(const Matter &matter) {
  nAtoms = matter.nAtoms;
  resize(nAtoms);

  positions = matter.positions;
  forces = matter.forces;
  masses = matter.masses;
  atomicNrs = matter.atomicNrs;
  isFixed = matter.isFixed;
  cell = matter.cell;
  cellInverse = matter.cellInverse;
  velocities = matter.velocities;

  removeNetForce = matter.removeNetForce;
  usePeriodicBoundaries = matter.usePeriodicBoundaries;

  potential = matter.potential;
  potentialEnergy = matter.potentialEnergy;
  recomputePotential = matter.recomputePotential;

  strcpy(headerCon1, matter.headerCon1);
  strcpy(headerCon2, matter.headerCon2);
  strcpy(headerCon4, matter.headerCon4);
  strcpy(headerCon5, matter.headerCon5);
  strcpy(headerCon6, matter.headerCon6);

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

    forces.resize(length, 3);
    forces.setZero();

    masses.resize(length);
    masses.setZero();

    atomicNrs.resize(length);
    atomicNrs.setZero();

    isFixed.resize(length);
    isFixed.setZero();
  }
  recomputePotential = true;
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
  recomputePotential = true;
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
  recomputePotential = true;
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
  recomputePotential = true;
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

// void Matter::setBiasPotential(BondBoost *bondBoost) {
//   biasPotential = bondBoost;
// }

// void Matter::setBiasForces(const AtomMatrix bf) {
//   biasForces = bf.array() * getFree().array();
// }
// return forces applied on all atoms in array 'force'
AtomMatrix Matter::getForces() {
  computePotential();
  AtomMatrix ret = forces;
  size_t i{0};
  for (i = 0; i < nAtoms; i++) {
    if (isFixed[i]) {
      ret.row(i).setZero();
    }
  }
  return ret;
}

VectorType Matter::getForcesV() {
  return VectorType::Map(getForces().data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getForcesFree() {
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

VectorType Matter::getForcesFreeV() {
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
  recomputePotential = true;
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

double Matter::getPotentialEnergy() {
  if (nAtoms > 0) {
    computePotential();
    return potentialEnergy;
  } else
    return 0.0;
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

void Matter::computePotential() {
  if (recomputePotential) {
    if (!potential) {
      throw(std::runtime_error("Whoops, you need a potential.."));
    }
    // auto surrogatePotential =
    //     std::dynamic_pointer_cast<SurrogatePotential>(potential);
    // if (surrogatePotential) {
    //   // Surrogate potential case
    //   auto [freePE, freeForces, vari] = surrogatePotential->get_ef_var(
    //       this->getPositionsFree(), this->getAtomicNrsFree().cast<int>(),
    //       cell);
    //   // Now populate full structures
    //   this->potentialEnergy = freePE;
    //   this->energyVariance = vari;
    //   for (size_t idx{0}, jdx{0}; idx < nAtoms; idx++) {
    //     if (!isFixed(idx)) {
    //       forces.row(idx) = freeForces.row(jdx);
    //       jdx++;
    //     }
    //   }
    // } else {
    // Non-surrogate potential case
    npotcalls++;
    auto [pE, frcs] = potential->get_ef(positions, atomicNrs, cell);
    potentialEnergy = pE;
    forces = frcs;
    // }
    // TODO(rg) :: Something about this isn't right..
    forceCalls += 1;
    recomputePotential = false;

    if (isFixed.sum() == 0 && removeNetForce) {
      FixedVecType<3> tempForce(3);
      tempForce = forces.colwise().sum() / nAtoms;

      for (size_t i = 0; i < nAtoms; i++) {
        forces.row(i) -= tempForce.transpose();
      }
    }
  }
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
  computePotential();
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

void Matter::setForces(const AtomMatrix f) {
  forces = f.array() * getFree().array();
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

bool Matter::convel2matter(std::string filename) {
  bool state;
  FILE *file;
  // Add the .con extension to filename if it is not already there.
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 3, "con")) {
    filename += ".con";
  }
  file = fopen(filename.c_str(), "rb");
  if (!file) {
    SPDLOG_LOGGER_ERROR(m_log, "File {} was not found.", filename);
    return (false);
  }
  state = convel2matter(file);
  fclose(file);
  return (state);
}

bool Matter::convel2matter(FILE *file) {
  char line[255]; // Temporary string of character to read from the file.
  fgets(headerCon1, sizeof(line), file);

  //    if (strchr(headerCon1,'\r')) {
  //        /* Files created on Windows or on Mac with Excell have carriage
  //        returns (\r) instead of or along with the new line charater (\n). C
  //        recognises only the \n as the end of line. */ cerr << "A carriage
  //        return ('\\r') has been detected. To work correctly, new lines
  //        should be indicated by the new line character (\\n)."; return false;
  //        // return false for error
  //    }

  size_t i{0}, j{0};

  fgets(headerCon2, sizeof(line), file);

  double lengths[3];
  // The third line contains the length of the periodic cell
  fgets(line, sizeof(line), file);
  sscanf(line, "%lf %lf %lf", &lengths[0], &lengths[1], &lengths[2]);

  double angles[3];
  fgets(headerCon4, sizeof(line), file);
  // The fourth line contains the angles of the cell vectors
  sscanf(headerCon4, "%lf %lf %lf", &angles[0], &angles[1], &angles[2]);

  if (angles[0] == 90.0 && angles[1] == 90.0 && angles[2] == 90.0) {
    cell(0, 0) = lengths[0];
    cell(1, 1) = lengths[1];
    cell(2, 2) = lengths[2];
  } else {
    angles[0] *= M_PI / 180.0;
    angles[1] *= M_PI / 180.0;
    angles[2] *= M_PI / 180.0;

    cell(0, 0) = 1.0;
    cell(1, 0) = cos(angles[0]);
    cell(1, 1) = sin(angles[0]);
    cell(2, 0) = cos(angles[1]);
    cell(2, 1) = (cos(angles[2]) - cell(1, 0) * cell(2, 0)) / cell(1, 1);
    cell(2, 2) = sqrt(1.0 - pow(cell(2, 0), 2) - pow(cell(2, 1), 2));

    cell(0, 0) *= lengths[0];
    cell(1, 0) *= lengths[1];
    cell(1, 1) *= lengths[1];
    cell(2, 0) *= lengths[2];
    cell(2, 1) *= lengths[2];
    cell(2, 2) *= lengths[2];
  }
  cellInverse = cell.inverse();

  fgets(headerCon5, sizeof(line), file);
  fgets(headerCon6, sizeof(line), file);

  fgets(line, sizeof(line), file);
  size_t Ncomponent; // Number of components or different types of atoms. For
                     // instance H2O has two components (H and O).
  if (sscanf(line, "%lud", &Ncomponent) == 0) {
    SPDLOG_LOGGER_INFO(m_log, "The number of components cannot be read. One "
                              "component is assumed instead");
    Ncomponent = 1;
  }
  if ((Ncomponent > MAXC) || (Ncomponent < 1)) {
    SPDLOG_LOGGER_ERROR(
        m_log,
        "con2atoms doesn't support more that {} components or less than 1",
        MAXC);
    return false;
  }
  /* to store the position of the
      first atom of each element 'MAXC+1': the last element is used to store the
     total number of atom.*/
  size_t first[MAXC + 1];
  size_t Natoms = 0;
  first[0] = 0;

  // Now we want to know the number of atom of each type. Ex with H2O, two
  // hydrogens and one oxygen
  for (j = 0; j < Ncomponent; j++) {
    fscanf(file, "%lud", &Natoms);
    first[j + 1] = Natoms + first[j];
  }

  fgets(line, sizeof(line), file); // Discard the rest of the line
  // Set the total number of atoms, and allocates memory
  resize(first[Ncomponent]);
  double mass[MAXC];
  for (j = 0; j < Ncomponent; j++) {
    // Now we want to know the number of atom of each type. Ex with
    // H2O, two hydrogens and one oxygen
    fscanf(file, "%lf", &mass[j]);
    // mass[j]*=G_PER_MOL; // conversion of g/mol to local units. (see su.h)
  }

  fgets(line, sizeof(line), file); // discard rest of the line
  int atomicNr;
  int fixed;
  double x, y, z;
  for (j = 0; j < Ncomponent; j++) {
    char symbol[3];
    fgets(line, sizeof(line), file);
    sscanf(line, "%2s\n", symbol);
    atomicNr = symbol2atomicNumber(symbol);
    fgets(line, sizeof(line), file); // skip one line
    for (i = first[j]; i < first[j + 1]; i++) {
      setMass(i, mass[j]);
      setAtomicNr(i, atomicNr);
      fgets(line, sizeof(line), file);
      sscanf(line, "%lf %lf %lf %d\n", &x, &y, &z, &fixed);
      setPosition(i, 0, x);
      setPosition(i, 1, y);
      setPosition(i, 2, z);
      setFixed(i, static_cast<bool>(fixed));
    }
  }

  fgets(line, sizeof(line), file);
  for (j = 0; j < Ncomponent; j++) {
    fgets(line, sizeof(line), file);
    fgets(line, sizeof(line), file); // skip one line
    for (i = first[j]; i < first[j + 1]; i++) {
      fgets(line, sizeof(line), file);
      sscanf(line, "%lf %lf %lf %d\n", &x, &y, &z, &fixed);
      setVelocity(i, 0, x);
      setVelocity(i, 1, y);
      setVelocity(i, 2, z);
    }
  }

  if (usePeriodicBoundaries) {
    applyPeriodicBoundary(); // Transform the coordinate to use the minimum
                             // image convention.
  }
  //    potential_ = new Potential(parameters_);
  return (true);
}

void Matter::setPotential(std::shared_ptr<PotBase> pot) {
  this->potential = pot;
  recomputePotential = true;
}

size_t Matter::getPotentialCalls() const { return this->npotcalls; }

double Matter::getEnergyVariance() { return this->energyVariance; }

// VectorType Matter::getForceVariance() {
//   return this->variance.segment(1, numberOfFreeAtoms() * 3);
// }

// double Matter::getMaxVariance() { return this->variance.maxCoeff(); }

std::shared_ptr<PotBase> Matter::getPotential() { return this->potential; }

} // namespace eonc
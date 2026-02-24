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
#include "HelperFunctions.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "SurrogatePotential.h"

// To write the R style data frame
#include <fmt/os.h>
#include <memory>
#include <spdlog/spdlog.h>
#include <stdexcept>

using namespace std;

namespace {

const char *elementArray[] = {
    "Unknown", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg",      "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr",      "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr",      "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd",      "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd",      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf",      "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po",      "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  NULL};

// guess the atom type from the atomic mass,
std::string mass2atom(double atomicmass) {
  return elementArray[int(atomicmass + .5)];
}

const int MAXC =
    100; // maximum number of components for functions matter2con and con2matter

int symbol2atomicNumber(char const *symbol) {
  int i = 0;

  while (elementArray[i] != NULL) {
    if (strcmp(symbol, elementArray[i]) == 0) {
      return i;
    }
    i++;
  }
  // invalid symbol
  return -1;
}

char const *atomicNumber2symbol(int n) { return elementArray[n]; }
} // namespace

class MatterObjectiveFunction : public ObjectiveFunction {
public:
  MatterObjectiveFunction(std::shared_ptr<Matter> matterPassed,
                          const Parameters &parametersPassed)
      : ObjectiveFunction(matterPassed, parametersPassed) {}
  ~MatterObjectiveFunction() = default;
  double getEnergy() { return matter->getPotentialEnergy(); }
  VectorXd getGradient(bool fdstep = false) {
    return -matter->getForcesFreeV();
  }
  void setPositions(VectorXd x) { matter->setPositionsFreeV(x); }
  VectorXd getPositions() { return matter->getPositionsFreeV(); }
  int degreesOfFreedom() { return 3 * matter->numberOfFreeAtoms(); }
  bool isConverged() {
    return getConvergence() < params.optimizer_options.converged_force;
  }
  double getConvergence() {
    if (params.optimizer_options.convergence_metric == "norm") {
      return matter->getForcesFreeV().norm();
    } else if (params.optimizer_options.convergence_metric == "max_atom") {
      return matter->maxForce();
    } else if (params.optimizer_options.convergence_metric == "max_component") {
      return matter->getForces().maxCoeff();
    } else {
      SPDLOG_CRITICAL("{} Unknown opt_convergence_metric: {}", "[Matter]"s,
                      params.optimizer_options.convergence_metric);
      std::exit(1);
    }
  }
  VectorXd difference(VectorXd a, VectorXd b) { return matter->pbcV(a - b); }
};

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

  parameters = matter.parameters;

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
//     if(parameters->structure_comparison_options.check_rotation) {
//         return helper_functions::rotationMatch(this, &matter,
//         parameters->structure_comparison_options.distance_difference);
//     }else{
//         return (parameters->structure_comparison_options.distance_difference)
//         > perAtomNorm(matter);
//     }
// }

bool Matter::compare(const Matter &matter, bool indistinguishable) {
  if (nAtoms != matter.numberOfAtoms())
    return false;
  if (parameters->structure_comparison_options.check_rotation &&
      indistinguishable) {
    return helper_functions::sortedR(
        *this, matter,
        parameters->structure_comparison_options.distance_difference);
  } else if (indistinguishable) {
    if (this->numberOfFixedAtoms() == 0 and
        parameters->structure_comparison_options.remove_translation)
      helper_functions::translationRemove(*this, matter);
    return helper_functions::identical(
        *this, matter,
        parameters->structure_comparison_options.distance_difference);
  } else if (parameters->structure_comparison_options.check_rotation) {
    return helper_functions::rotationMatch(
        *this, matter,
        parameters->structure_comparison_options.distance_difference);
  } else {
    if (this->numberOfFixedAtoms() == 0 and
        parameters->structure_comparison_options.remove_translation)
      helper_functions::translationRemove(*this, matter);
    return (parameters->structure_comparison_options.distance_difference) >
           perAtomNorm(matter);
  }
}

// bool Matter::operator!=(const Matter& matter) {
//     return !operator==(matter);
// }

// Returns the distance to the given matter object.
double Matter::distanceTo(const Matter &matter) {
  return pbc(positions - matter.positions).norm();
}

AtomMatrix Matter::pbc(const AtomMatrix &diff) const {
  AtomMatrix ddiff = diff * cellInverse;

  int i, j;
  for (i = 0; i < diff.rows(); i++) {
    for (j = 0; j < 3; j++) {
      ddiff(i, j) = fmod(fmod(ddiff(i, j), 1.0) + 1.5, 1.0) - .5;
    }
  }

  return ddiff * cell;
}

VectorXd Matter::pbcV(const VectorXd &diffVector) const {
  AtomMatrix pbcMatrix =
      pbc(AtomMatrix::Map(diffVector.data(), diffVector.size() / 3, 3));
  return VectorXd(VectorXd::Map(pbcMatrix.data(), diffVector.size()));
}

// Returns the maximum distance between two atoms in the Matter objects.
double Matter::perAtomNorm(const Matter &matter) {
  long i = 0;
  double max_distance = 0.0;

  if (matter.numberOfAtoms() == nAtoms) {
    AtomMatrix diff = pbc(positions - matter.positions);
    for (i = 0; i < nAtoms; i++) {
      max_distance = max(diff.row(i).norm(), max_distance);
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

long int Matter::numberOfAtoms() const { return (nAtoms); }

Matrix3d Matter::getCell() const { return cell; }

void Matter::setCell(const Matrix3d &newCell) { cell = newCell; }

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

// return coordinates of atoms by const reference (zero-copy)
const AtomMatrix &Matter::getPositions() const { return positions; }
// return a modifiable copy of positions
AtomMatrix Matter::getPositionsCopy() const { return positions; }

VectorXd Matter::getPositionsV() const {
  return VectorXd::Map(positions.data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getPositionsFree() const {
  AtomMatrix ret(numberOfFreeAtoms(), 3);
  int i, j = 0;
  for (i = 0; i < nAtoms; i++) {
    if (!isFixed(i)) {
      ret.row(j) = positions.row(i);
      j++;
    }
  }
  return ret;
}

VectorXi Matter::getAtomicNrsFree() const {
  return this->atomicNrs.array() * getFreeV().cast<int>().array();
}

bool Matter::relax(bool quiet, bool writeMovie, bool checkpoint,
                   string prefixMovie, string prefixCheckpoint) {
  auto objf = std::make_shared<MatterObjectiveFunction>(
      std::make_shared<Matter>(*this), *parameters);
  auto optim = helpers::create::mkOptim(
      objf, parameters->optimizer_options.method, *parameters);

  ostringstream min;
  min << prefixMovie;
  if (writeMovie) {
    matter2con(min.str(), false);
  }

  int iteration = 0;
  if (!quiet) {
    SPDLOG_LOGGER_DEBUG(m_log, "{} {:10s}  {:14s}  {:18s}  {:13s}\n",
                        "[Matter]", "Iter", "Step size",
                        parameters->optimizer_options.convergence_metric_label,
                        "Energy");
    SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}\n",
                        "[Matter]", iteration, 0.0, objf->getConvergence(),
                        getPotentialEnergy());
  }

  while (!objf->isConverged() &&
         iteration < parameters->optimizer_options.max_iterations) {

    AtomMatrix pos = getPositions();

    optim->step(parameters->optimizer_options.max_move);
    iteration++;
    setPositionsFreeV(objf->getPositions());

    double stepSize =
        helper_functions::maxAtomMotion(pbc(getPositions() - pos));

    if (!quiet) {
      SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                          "[Matter]", iteration, stepSize,
                          objf->getConvergence(), getPotentialEnergy());
    }

    if (writeMovie) {
      matter2con(min.str(), true);
    }

    if (checkpoint) {
      ostringstream chk;
      chk << prefixCheckpoint << "_cp";
      matter2con(chk.str(), false);
    }
  }

  if (iteration == 0) {
    if (!quiet) {
      SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                          "[Matter]", iteration, 0.0, objf->getConvergence(),
                          getPotentialEnergy());
    }
  }
  //    bool converged =
  //    optimizer->run(parameters->optimizer_options.max_iterations,
  //    parameters->optimizer_options.max_move);
  return objf->isConverged();
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
}

// Same but takes vector instead of n x 3 matrix
void Matter::setPositionsV(const VectorXd &pos) {
  setPositions(AtomMatrix::Map(pos.data(), numberOfAtoms(), 3));
}

void Matter::setPositionsFree(const AtomMatrix &pos) {
  // FIXME: Ensure pos and existing data are in the same form with atom ids
  int i, j = 0;
  for (i = 0; i < nAtoms; i++) {
    if (!isFixed(i)) {
      positions.row(i) = pos.row(j);
      j += 1;
    }
  }
  recomputePotential = true;
}

void Matter::setPositionsFreeV(const VectorXd &pos) {
  setPositionsFree(AtomMatrix::Map(pos.data(), numberOfFreeAtoms(), 3));
}

AtomMatrix Matter::getBiasForces() {
  if (biasPotential != NULL) {
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
// return forces applied on all atoms in array 'force'
AtomMatrix Matter::getForces() {
  computePotential();
  AtomMatrix ret = forces;
  int i;
  for (i = 0; i < nAtoms; i++) {
    if (isFixed[i]) {
      ret.row(i).setZero();
    }
  }
  return ret;
}

VectorXd Matter::getForcesV() {
  return VectorXd::Map(getForces().data(), 3 * numberOfAtoms());
}

AtomMatrix Matter::getForcesFree() {
  AtomMatrix allForces = getForces();
  AtomMatrix ret(numberOfFreeAtoms(), 3);
  int i, j = 0;
  for (i = 0; i < nAtoms; i++) {
    if (!isFixed[i]) {
      ret.row(j) = allForces.row(i);
      j++;
    }
  }
  return ret;
}

VectorXd Matter::getForcesFreeV() {
  return VectorXd::Map(getForcesFree().data(), 3 * numberOfFreeAtoms());
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
}

int Matter::getFixed(long int indexAtom) const { return (isFixed[indexAtom]); }

void Matter::setFixed(long int indexAtom, int isFixed_passed) {
  isFixed[indexAtom] = isFixed_passed;
}

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
  for (long int i = 0; i < nAtoms; i++) {
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

long Matter::getForceCalls() const { return (forceCalls); }

void Matter::resetForceCalls() {
  forceCalls = 0;
  return;
}

// Print atomic coordinate to a .xyz file
void Matter::matter2xyz(std::string filename,
                        bool append /*Append if file already exist*/) {
  FILE *file;
  long int i;
  filename += ".xyz";
  if (append) {
    file = fopen(filename.c_str(), "ab");
  } else {
    file = fopen(filename.c_str(), "wb");
  }
  if (file == 0) {
    cerr << "Can't create file " << filename << endl;
    exit(1);
  }
  fprintf(file, "%ld\nGenerated by EON\n", numberOfAtoms());

  if (usePeriodicBoundaries) {
    applyPeriodicBoundary(); // Transform the coordinate to use the minimum
                             // image convention.
  }

  for (i = 0; i < numberOfAtoms(); i++) {
    fprintf(file, "%s\t%11.6f\t%11.6f\t%11.6f\n",
            atomicNumber2symbol(getAtomicNr(i)), getPosition(i, 0),
            getPosition(i, 1), getPosition(i, 2));
  }
  fclose(file);
}

std::pair<std::array<double, 3>, std::array<double, 3>>
Matter::cell_to_lengths_angles() const {
  std::array<double, 3> lengths;
  lengths[0] = cell.row(0).norm();
  lengths[1] = cell.row(1).norm();
  lengths[2] = cell.row(2).norm();
  std::array<double, 3> angles;
  angles[0] = acos(cell.row(0).dot(cell.row(1)) / lengths[0] / lengths[1]) *
              180.0 / helper_functions::pi;
  angles[1] = acos(cell.row(0).dot(cell.row(2)) / lengths[0] / lengths[2]) *
              180.0 / helper_functions::pi;
  angles[2] = acos(cell.row(1).dot(cell.row(2)) / lengths[1] / lengths[2]) *
              180.0 / helper_functions::pi;
  return {lengths, angles};
}

// Print atomic coordinates to a .con file
bool Matter::matter2con(std::string filename, bool append) {
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 3, "con")) {
    filename += ".con";
  }

  if (usePeriodicBoundaries) {
    applyPeriodicBoundary();
  }

  auto [lengths, angles_deg] = cell_to_lengths_angles();

  // Strip trailing newlines from header buffers for readcon
  auto strip_nl = [](const char *s) -> std::string {
    std::string str(s);
    while (!str.empty() && (str.back() == '\n' || str.back() == '\r'))
      str.pop_back();
    return str;
  };

  readcon::ConFrameBuilder builder(
      {lengths[0], lengths[1], lengths[2]},
      {angles_deg[0], angles_deg[1], angles_deg[2]},
      {strip_nl(headerCon1), strip_nl(headerCon2)},
      {strip_nl(headerCon5), strip_nl(headerCon6)});

  for (long i = 0; i < numberOfAtoms(); i++) {
    builder.add_atom(atomicNumber2symbol(getAtomicNr(i)), getPosition(i, 0),
                     getPosition(i, 1), getPosition(i, 2), getFixed(i) != 0,
                     static_cast<uint64_t>(i), getMass(i));
  }

  auto frame = builder.build();
  std::vector<readcon::ConFrame> frames;
  if (append) {
    try {
      frames = readcon::read_all_frames(filename);
    } catch (...) {
      // File doesn't exist yet, start fresh
    }
  }
  frames.push_back(std::move(frame));
  readcon::ConFrameWriter writer(filename, 17);
  writer.extend(frames);
  return true;
}

// Load atomic coordinates from a .con file via readcon-core (mmap reader)
bool Matter::con2matter(std::string filename) {
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 3, "con")) {
    filename += ".con";
  }
  try {
    auto frame = readcon::read_first_frame(filename);
    return con2matter(frame);
  } catch (const std::exception &e) {
    SPDLOG_LOGGER_ERROR(m_log, "Failed to read {}: {}", filename, e.what());
    return false;
  }
}

// Populate Matter from a parsed readcon frame
bool Matter::con2matter(const readcon::ConFrame &frame) {
  const auto &atoms = frame.atoms();
  const auto &lengths = frame.cell();
  const auto &angles_deg = frame.angles();
  const auto &prebox = frame.prebox_header();
  const auto &postbox = frame.postbox_header();

  // Store headers for round-tripping via matter2con
  snprintf(headerCon1, sizeof(headerCon1), "%s\n", prebox[0].c_str());
  snprintf(headerCon2, sizeof(headerCon2), "%s\n", prebox[1].c_str());
  snprintf(headerCon5, sizeof(headerCon5), "%s\n", postbox[0].c_str());
  snprintf(headerCon6, sizeof(headerCon6), "%s\n", postbox[1].c_str());

  // Build cell matrix from lengths and angles
  double angles[3] = {angles_deg[0], angles_deg[1], angles_deg[2]};
  if (angles[0] == 90.0 && angles[1] == 90.0 && angles[2] == 90.0) {
    cell.setZero();
    cell(0, 0) = lengths[0];
    cell(1, 1) = lengths[1];
    cell(2, 2) = lengths[2];
  } else {
    angles[0] *= helper_functions::pi / 180.0;
    angles[1] *= helper_functions::pi / 180.0;
    angles[2] *= helper_functions::pi / 180.0;

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

  // Also store headerCon4 for convel round-tripping
  snprintf(headerCon4, sizeof(headerCon4), "%f %f %f\n", angles_deg[0],
           angles_deg[1], angles_deg[2]);

  resize(static_cast<long>(atoms.size()));
  for (size_t i = 0; i < atoms.size(); ++i) {
    positions(i, 0) = atoms[i].x;
    positions(i, 1) = atoms[i].y;
    positions(i, 2) = atoms[i].z;
    setMass(i, atoms[i].mass);
    setAtomicNr(i, static_cast<int>(atoms[i].atomic_number));
    setFixed(i, atoms[i].is_fixed ? 1 : 0);
  }

  if (usePeriodicBoundaries) {
    applyPeriodicBoundary();
  }
  recomputePotential = true;
  return true;
}

void Matter::computePotential() {
  if (recomputePotential) {
    if (!potential) {
      throw(std::runtime_error("Whoops, you need a potential.."));
      potential = helper_functions::makePotential(
          parameters->potential_options.potential, *parameters);
    }
    auto surrogatePotential =
        std::dynamic_pointer_cast<SurrogatePotential>(potential);
    if (surrogatePotential) {
      // Surrogate potential case
      auto [freePE, freeForces, vari] = surrogatePotential->get_ef_var(
          this->getPositionsFree(), this->getAtomicNrsFree(), cell);
      // Now populate full structures
      this->potentialEnergy = freePE;
      this->energyVariance = vari;
      for (long idx{0}, jdx{0}; idx < nAtoms; idx++) {
        if (!isFixed(idx)) {
          forces.row(idx) = freeForces.row(jdx);
          jdx++;
        }
      }
    } else {
      // Non-surrogate potential case
      auto [pE, frcs] = potential->get_ef(positions, atomicNrs, cell);
      potentialEnergy = pE;
      forces = frcs;
    }
    forceCalls = forceCalls + 1;
    recomputePotential = false;

    if (isFixed.sum() == 0 && parameters->main_options.removeNetForce) {
      Vector3d tempForce(3);
      tempForce = forces.colwise().sum() / nAtoms;

      for (long int i = 0; i < nAtoms; i++) {
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

  // I think this can be done in one line with the rowwise method
  double maxForce = 0.0;
  for (int i = 0; i < nAtoms; i++) {
    if (getFixed(i)) {
      continue;
    }
    maxForce = max(forces.row(i).norm(), maxForce);
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
  AtomMatrix ret(nAtoms, 3);
  int i, j;
  for (i = 0; i < nAtoms; i++) {
    for (j = 0; j < 3; j++) {
      ret(i, j) = double(!bool(isFixed(i)));
    }
  }
  return ret;
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
  ret.col(0).array() /= masses.array();
  ret.col(1).array() /= masses.array();
  ret.col(2).array() /= masses.array();
  return ret;
}

Matrix<double, Eigen::Dynamic, 1> Matter::getMasses() const { return masses; }

bool Matter::matter2convel(std::string filename) {
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 6, "convel")) {
    filename += ".convel";
  }

  if (usePeriodicBoundaries) {
    applyPeriodicBoundary();
  }

  auto [lengths, angles_deg] = cell_to_lengths_angles();

  auto strip_nl = [](const char *s) -> std::string {
    std::string str(s);
    while (!str.empty() && (str.back() == '\n' || str.back() == '\r'))
      str.pop_back();
    return str;
  };

  readcon::ConFrameBuilder builder(
      {lengths[0], lengths[1], lengths[2]},
      {angles_deg[0], angles_deg[1], angles_deg[2]},
      {strip_nl(headerCon1), strip_nl(headerCon2)},
      {strip_nl(headerCon5), strip_nl(headerCon6)});

  for (long i = 0; i < numberOfAtoms(); i++) {
    builder.add_atom_with_velocity(
        atomicNumber2symbol(getAtomicNr(i)), getPosition(i, 0),
        getPosition(i, 1), getPosition(i, 2), getFixed(i) != 0,
        static_cast<uint64_t>(i), getMass(i), velocities(i, 0),
        velocities(i, 1), velocities(i, 2));
  }

  auto frame = builder.build();
  readcon::ConFrameWriter writer(filename, 6);
  std::vector<readcon::ConFrame> frames;
  frames.push_back(std::move(frame));
  writer.extend(frames);
  return true;
}

bool Matter::convel2matter(std::string filename) {
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 3, "con")) {
    filename += ".con";
  }
  try {
    auto frame = readcon::read_first_frame(filename);
    bool ok = con2matter(frame);
    if (!ok)
      return false;
    if (frame.has_velocities()) {
      const auto &atoms = frame.atoms();
      for (size_t i = 0; i < atoms.size(); ++i) {
        setVelocity(i, 0, atoms[i].vx);
        setVelocity(i, 1, atoms[i].vy);
        setVelocity(i, 2, atoms[i].vz);
      }
    }
    return true;
  } catch (const std::exception &e) {
    SPDLOG_LOGGER_ERROR(m_log, "Failed to read convel {}: {}", filename,
                        e.what());
    return false;
  }
}

void Matter::writeTibble(std::string fname) {
  using namespace fmt::literals;
  auto out = fmt::output_file(fname);
  AtomMatrix fSys = this->getForces();
  double eSys = this->getPotentialEnergy();
  AtomMatrix pos = this->getPositions();
  out.print("x y z fx fy fz energy mass symbol atmID fixed\n");
  for (auto idx{0}; idx < this->numberOfAtoms(); idx++) {
    out.print(
        "{x} {y} {z} {fx} {fy} {fz} {energy} {mass} {symbol} {idx} {fixed}\n",
        "x"_a = pos.row(idx)[0], "y"_a = pos.row(idx)[1],
        "z"_a = pos.row(idx)[2], "fx"_a = fSys.row(idx)[0],
        "fy"_a = fSys.row(idx)[1], "fz"_a = fSys.row(idx)[2], "energy"_a = eSys,
        "mass"_a = this->getMass(idx),
        "symbol"_a = atomicNumber2symbol(this->getAtomicNr(idx)),
        "idx"_a = (idx + 1),
        "fixed"_a =
            this->getFixed(idx)); // NOTE: idx MAY not be the same id as before
  }
  return;
}

void Matter::setPotential(std::shared_ptr<Potential> pot) {
  this->potential = pot;
  recomputePotential = true;
}

size_t Matter::getPotentialCalls() const {
  return this->potential->forceCallCounter;
}

double Matter::getEnergyVariance() { return this->energyVariance; }

// Eigen::VectorXd Matter::getForceVariance() {
//   return this->variance.segment(1, numberOfFreeAtoms() * 3);
// }

// double Matter::getMaxVariance() { return this->variance.maxCoeff(); }

std::shared_ptr<Potential> Matter::getPotential() { return this->potential; }

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
                          std::shared_ptr<Parameters> parametersPassed)
      : ObjectiveFunction(matterPassed, parametersPassed) {}
  ~MatterObjectiveFunction() = default;
  double getEnergy() { return matter->getPotentialEnergy(); }
  VectorXd getGradient(bool fdstep = false) {
    return -matter->getForcesFreeV();
  }
  void setPositions(VectorXd x) { matter->setPositionsFreeV(x); }
  VectorXd getPositions() { return matter->getPositionsFreeV(); }
  int degreesOfFreedom() { return 3 * matter->numberOfFreeAtoms(); }
  bool isConverged() { return getConvergence() < params->optConvergedForce; }
  double getConvergence() {
    if (params->optConvergenceMetric == "norm") {
      return matter->getForcesFreeV().norm();
    } else if (params->optConvergenceMetric == "max_atom") {
      return matter->maxForce();
    } else if (params->optConvergenceMetric == "max_component") {
      return matter->getForces().maxCoeff();
    } else {
      SPDLOG_CRITICAL("{} Unknown opt_convergence_metric: {}", "[Matter]"s,
                      params->optConvergenceMetric);
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
//     if(parameters->checkRotation) {
//         return helper_functions::rotationMatch(this, &matter,
//         parameters->distanceDifference);
//     }else{
//         return (parameters->distanceDifference) > perAtomNorm(matter);
//     }
// }

bool Matter::compare(const Matter &matter, bool indistinguishable) {
  if (nAtoms != matter.numberOfAtoms())
    return false;
  if (parameters->checkRotation && indistinguishable) {
    return helper_functions::sortedR(*this, matter,
                                     parameters->distanceDifference);
  } else if (indistinguishable) {
    if (this->numberOfFixedAtoms() == 0 and parameters->removeTranslation)
      helper_functions::translationRemove(*this, matter);
    return helper_functions::identical(*this, matter,
                                       parameters->distanceDifference);
  } else if (parameters->checkRotation) {
    return helper_functions::rotationMatch(*this, matter,
                                           parameters->distanceDifference);
  } else {
    if (this->numberOfFixedAtoms() == 0 and parameters->removeTranslation)
      helper_functions::translationRemove(*this, matter);
    return (parameters->distanceDifference) > perAtomNorm(matter);
  }
}

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

VectorXd Matter::pbcV(VectorXd diffVector) const {
  AtomMatrix pbcMatrix =
      pbc(AtomMatrix::Map(diffVector.data(), diffVector.size() / 3, 3));
  return VectorXd::Map(pbcMatrix.data(), diffVector.size());
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

void Matter::setCell(Matrix3d newCell) { cell = newCell; }

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
      std::make_shared<Matter>(*this), parameters);
  auto optim =
      helpers::create::mkOptim(objf, parameters->optMethod, parameters);

  ostringstream min;
  min << prefixMovie;
  if (writeMovie) {
    matter2con(min.str(), false);
  }

  int iteration = 0;
  if (!quiet) {
    SPDLOG_LOGGER_DEBUG(m_log, "{} {:10s}  {:14s}  {:18s}  {:13s}\n",
                        "[Matter]", "Iter", "Step size",
                        parameters->optConvergenceMetricLabel, "Energy");
    SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}\n",
                        "[Matter]", iteration, 0.0, objf->getConvergence(),
                        getPotentialEnergy());
  }

  while (!objf->isConverged() && iteration < parameters->optMaxIterations) {

    AtomMatrix pos = getPositions();

    optim->step(parameters->optMaxMove);
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
  //    bool converged = optimizer->run(parameters->optMaxIterations,
  //    parameters->optMaxMove);
  return objf->isConverged();
}

VectorXd Matter::getPositionsFreeV() const {
  return VectorXd::Map(getPositionsFree().data(), 3 * numberOfFreeAtoms());
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
void Matter::setPositionsV(const VectorXd pos) {
  setPositions(AtomMatrix::Map(pos.data(), numberOfAtoms(), 3));
}

void Matter::setPositionsFree(const AtomMatrix pos) {
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

void Matter::setPositionsFreeV(const VectorXd pos) {
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

void Matter::setBiasForces(const AtomMatrix bf) {
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

void Matter::setMasses(VectorXd massesIn) {
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

// Print atomic coordinates to a .con file
bool Matter::matter2con(std::string filename, bool append) {
  bool state;
  FILE *file;
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 3, "con")) {
    filename += ".con";
  };
  if (append) {
    file = fopen(filename.c_str(), "ab");
  } else {
    file = fopen(filename.c_str(), "wb");
  }
  state = matter2con(file);
  fclose(file);
  return (state);
}

bool Matter::matter2con(FILE *file) {
  long int i;
  int j;
  long int Nfix = 0;  // Nfix to store the number of fixed atoms
  int Ncomponent = 0; // used to store the number of components (eg water: two
                      // components H and O)
  int first[MAXC]; // to store the position of the first atom of each component
                   // plus at the end the total number of atoms
  double mass[MAXC];
  long atomicNrs[MAXC];
  first[0] = 0;

  if (usePeriodicBoundaries) {
    applyPeriodicBoundary(); // Transform the coordinate to use the minimum
                             // image convention.
  }

  if (numberOfAtoms() > 0) {
    if (getFixed(0))
      Nfix = 1; // count the number of fixed atoms
    mass[0] = getMass(0);
    atomicNrs[0] = getAtomicNr(0);
  };
  j = 0;
  for (i = 1; i < numberOfAtoms(); i++) {
    if (getFixed(i))
      Nfix++; // count the number of fixed atoms
    if (getAtomicNr(i) !=
        atomicNrs[j]) { // check if there is a second component
      j++;
      if (j >= MAXC) {
        SPDLOG_LOGGER_ERROR(m_log,
                            "Does not support more than {} components and the "
                            "atoms must be ordered by component.",
                            MAXC);
        return false;
      };
      mass[j] = getMass(i);
      atomicNrs[j] = getAtomicNr(i);
      first[j] = i;
    }
  }
  first[j + 1] = numberOfAtoms();
  Ncomponent = j + 1;

  fputs(headerCon1, file);
  fputs(headerCon2, file);
  double lengths[3];
  lengths[0] = cell.row(0).norm();
  lengths[1] = cell.row(1).norm();
  lengths[2] = cell.row(2).norm();
  fprintf(file, "%f\t%f\t%f\n", lengths[0], lengths[1], lengths[2]);
  double angles[3];
  angles[0] =
      acos(cell.row(0).dot(cell.row(1)) / lengths[0] / lengths[1]) * 180 / M_PI;
  angles[1] =
      acos(cell.row(0).dot(cell.row(2)) / lengths[0] / lengths[2]) * 180 / M_PI;
  angles[2] =
      acos(cell.row(1).dot(cell.row(2)) / lengths[1] / lengths[2]) * 180 / M_PI;
  fprintf(file, "%f\t%f\t%f\n", angles[0], angles[1], angles[2]);
  fputs(headerCon5, file);
  fputs(headerCon6, file);

  fprintf(file, "%d\n", Ncomponent);
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%d ", first[j + 1] - first[j]);
  }
  fprintf(file, "\n");
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%f ", mass[j]);
  }
  fprintf(file, "\n");
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%s\n", atomicNumber2symbol(atomicNrs[j]));
    fprintf(file, "Coordinates of Component %d\n", j + 1);
    for (i = first[j]; i < first[j + 1]; i++) {
      fprintf(file, "%22.17f %22.17f %22.17f %d %4ld\n", getPosition(i, 0),
              getPosition(i, 1), getPosition(i, 2), getFixed(i), i);
    }
  }
  return true;
}

// Load atomic coordinates from a .con file
bool Matter::con2matter(std::string filename) {
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
  state = con2matter(file);
  fclose(file);
  return (state);
}

bool Matter::con2matter(FILE *file) {
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

  long int i;
  int j;

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
  int Ncomponent; // Number of components or different types of atoms  (eg
                  // water: two components H and O)
  if (sscanf(line, "%d", &Ncomponent) == 0) {
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

  // stores the position of the first atom of each element
  // 'MAXC+1': the last element is used to store the total number of atom
  long int first[MAXC + 1];
  long int Natoms = 0;
  first[0] = 0;

  // Now we want to know the number of atom of each type. Ex with H2O, two
  // hydrogens and one oxygen
  fgets(line, sizeof(line), file);
  // split at either space or tab
  char *split = strtok(line, " \t");
  for (j = 0; j < Ncomponent; j++) {
    if (split == NULL) {
      SPDLOG_LOGGER_ERROR(
          m_log, "input con file does not list the number of each component");
      return false;
    }
    if (sscanf(split, "%ld", &Natoms) != 1) {
      SPDLOG_LOGGER_ERROR(
          m_log, "input con file does not list the number of each component");
      return false;
    }
    first[j + 1] = Natoms + first[j];
    // split at either space or tab
    split = strtok(NULL, " \t");
  }

  // Set the total number of atoms, and allocates memory
  resize(first[Ncomponent]);

  double mass[MAXC];
  fgets(line, sizeof(line), file);
  // split at either space or tab
  split = strtok(line, " \t");

  for (j = 0; j < Ncomponent; j++) {
    // Now we want to know the number of atom of each type. Ex with
    // H2O, two hydrogens and one oxygen
    if (split == NULL) {
      SPDLOG_LOGGER_ERROR(m_log, "input con file does not list enough masses");
      return false;
    }
    // *1* seems like a bug as a result of copying and pasting from above
    // *1*       if(sscanf(split, "%ld", &Natoms)!=1)
    if (sscanf(split, "%lf", &mass[j]) != 1) {
      SPDLOG_LOGGER_ERROR(m_log, "input con file does not list enough masses");
      return false;
    }
    // *1*       sscanf(line, "%lf", &mass[j]);
    // split at either space or tab
    split = strtok(NULL, " \t");
  }

  int atomicNr;
  int fixed;
  double x, y, z;
  for (j = 0; j < Ncomponent; j++) {
    char symbol[3];
    fgets(line, sizeof(line), file);
    sscanf(line, "%2s", symbol);
    atomicNr = symbol2atomicNumber(symbol);
    fgets(line, sizeof(line), file); // skip one line
    for (i = first[j]; i < first[j + 1]; i++) {
      setMass(i, mass[j]);

      setAtomicNr(i, atomicNr);
      fgets(line, sizeof(line), file);
      if (strlen(line) < 6) {
        SPDLOG_LOGGER_ERROR(m_log, "error parsing position in con file");
        return false;
      }

      sscanf(line, "%lf %lf %lf %d\n", &x, &y, &z, &fixed);
      positions(i, 0) = x;
      positions(i, 1) = y;
      positions(i, 2) = z;
      setFixed(i, static_cast<bool>(fixed));
    }
  }
  if (usePeriodicBoundaries) {
    // Transform the coordinate to use the minimum
    // image convention.
    applyPeriodicBoundary();
  }
  // potential_ = new Potential(parameters_);
  recomputePotential = true;
  return (true);
}

void Matter::computePotential() {
  if (recomputePotential) {
    if (!potential) {
      throw(std::runtime_error("Whoops, you need a potential.."));
      potential =
          helper_functions::makePotential(parameters->potential, parameters);
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

    if (isFixed.sum() == 0 && parameters->removeNetForce) {
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

void Matter::setAtomicNrs(const VectorXi atmnrs) {
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

void Matter::setVelocities(const AtomMatrix v) {
  velocities = v.array() * getFree().array();
}

void Matter::setForces(const AtomMatrix f) {
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
  bool state;
  FILE *file;
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 6, "convel")) {
    filename += ".convel";
  }
  file = fopen(filename.c_str(), "w");
  state = matter2convel(file);
  fclose(file);
  return (state);
}

bool Matter::matter2convel(FILE *file) {
  long int i;
  int j;
  long int Nfix = 0;  // Nfix to store the number of fixed atoms
  int Ncomponent = 0; // used to store the number of components (eg water: two
                      // components H and O)
  int first[MAXC]; // to store the position of the first atom of each component
                   // plus at the end the total number of atoms
  double mass[MAXC];
  long atomicNrs[MAXC];
  first[0] = 0;

  if (usePeriodicBoundaries) {
    applyPeriodicBoundary(); // Transform the coordinate to use the minimum
                             // image convention.
  }

  if (numberOfAtoms() > 0) {
    if (getFixed(0))
      Nfix = 1; // count the number of fixed atoms
    mass[0] = getMass(0);
    atomicNrs[0] = getAtomicNr(0);
  }
  j = 0;
  for (i = 1; i < numberOfAtoms(); i++) {
    if (getFixed(i))
      Nfix++; // count the number of fixed atoms
    if (getAtomicNr(i) !=
        atomicNrs[j]) { // check if there is a second component
      j++;
      if (j >= MAXC) {
        SPDLOG_LOGGER_ERROR(m_log,
                            "Does not support more than {} components and the "
                            "atoms must be ordered by component.",
                            MAXC);
        return false;
      }
      mass[j] = getMass(i);
      atomicNrs[j] = getAtomicNr(i);
      first[j] = i;
    }
  }
  first[j + 1] = numberOfAtoms();
  Ncomponent = j + 1;

  fputs(headerCon1, file);
  fputs(headerCon2, file);
  double lengths[3];
  lengths[0] = cell.row(0).norm();
  lengths[1] = cell.row(1).norm();
  lengths[2] = cell.row(2).norm();
  fprintf(file, "%f\t%f\t%f\n", lengths[0], lengths[1], lengths[2]);
  double angles[3];
  angles[0] =
      acos(cell.row(0).dot(cell.row(1)) / lengths[0] / lengths[1]) * 180 / M_PI;
  angles[1] =
      acos(cell.row(0).dot(cell.row(2)) / lengths[0] / lengths[2]) * 180 / M_PI;
  angles[2] =
      acos(cell.row(1).dot(cell.row(2)) / lengths[1] / lengths[2]) * 180 / M_PI;
  fprintf(file, "%f\t%f\t%f\n", angles[0], angles[1], angles[2]);
  fputs(headerCon5, file);
  fputs(headerCon6, file);

  fprintf(file, "%d\n", Ncomponent);
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%d ", first[j + 1] - first[j]);
  }
  fprintf(file, "\n");
  for (j = 0; j < Ncomponent; j++) {
    // mass[j]/=G_PER_MOL; // GH: I don't understand why we need to convert the
    // mass units
    fprintf(file, "%f ", mass[j]);
  }
  fprintf(file, "\n");
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%s\n", atomicNumber2symbol(atomicNrs[j]));
    fprintf(file, "Coordinates of Component %d\n", j + 1);
    for (i = first[j]; i < first[j + 1]; i++) {
      fprintf(file, "%11.6f\t%11.6f\t%11.6f\t%d\t%ld\n", getPosition(i, 0),
              getPosition(i, 1), getPosition(i, 2), getFixed(i), i);
    }
  }
  fprintf(file, "\n");
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%s\n", atomicNumber2symbol(atomicNrs[j]));
    fprintf(file, "Velocities of Component %d\n", j + 1);
    for (i = first[j]; i < first[j + 1]; i++) {
      fprintf(file, "%11.6f\t%11.6f\t%11.6f\t%d\t%ld\n", velocities(i, 0),
              velocities(i, 1), velocities(i, 2), getFixed(i), i);
    }
  }
  return true;
}

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

  long int i;
  int j;

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
  int Ncomponent; // Number of components or different types of atoms. For
                  // instance H2O has two components (H and O).
  if (sscanf(line, "%d", &Ncomponent) == 0) {
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
  long int first[MAXC + 1];
  long int Natoms = 0;
  first[0] = 0;

  // Now we want to know the number of atom of each type. Ex with H2O, two
  // hydrogens and one oxygen
  for (j = 0; j < Ncomponent; j++) {
    fscanf(file, "%ld", &Natoms);
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

#include "Matter.h"

#include "BondBoost.h"
#include "HelperFunctions.h"
#include "Log.h"
#include "Optimizer.h"
#include "StringHelpers.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <utility>

#include <fmt/core.h>
#include <fmt/printf.h>

using namespace std;
namespace fs = std::filesystem;

static const std::string LOG_PREFIX = "[Matter]"s;

// This is anonymous to only allow access from within this file
namespace {

// guess the atom type from the atomic mass,
std::string_view mass2atom(double atomicmass) { return elements::symbols[int(atomicmass + .5)]; }

const int MAXC = 100; // maximum number of components for functions matter2con and con2matter

size_t symbol2atomicNumber(std::string_view const symbol) {
    for (size_t idx = 0; auto elem : elements::symbols) {
        if (symbol.compare(elem) == 0) {
            return idx;
        }
        ++idx;
    }
    return -1;
}

const std::string_view atomicNumber2symbol(int n) { return elements::symbols[n]; }
} // namespace

Matter::Matter(Parameters *parameters) { initializeDataMembers(parameters); }

Matter::Matter(Parameters *parameters, const long int nAtoms) {
    resize(nAtoms); // prepare memory for nAtoms
    initializeDataMembers(parameters);
}

void Matter::initializeDataMembers(Parameters *params) {
    nAtoms = 0;
    biasPotential = NULL;
    cell.resize(3, 3);
    cell.setZero();
    cellInverse.resize(3, 3);
    cellInverse.setZero();
    usePeriodicBoundaries = true;
    recomputePotential = true;
    forceCalls = 0;
    parameters = params;
    potential = NULL;
}

Matter::Matter(const Matter &matter) { operator=(matter); }

Matter::~Matter() {
    //    if (potential!=NULL){
    //        delete potential;
    //    }
}

const Matter &Matter::operator=(const Matter &matter) {
    resize(matter.numberOfAtoms());
    nAtoms = matter.nAtoms;

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

    potentialEnergy = matter.potentialEnergy;
    recomputePotential = matter.recomputePotential;

    headerCon1 = matter.headerCon1;
    headerCon2 = matter.headerCon2;
    headerCon4 = matter.headerCon4;
    headerCon5 = matter.headerCon5;
    headerCon6 = matter.headerCon6;

    return *this;
}

// The == comparison considers identity. This is crucial for process search.
// bool Matter::operator==(const Matter& matter) {
//    if(parameters->checkRotation) {
//        return helper_functions::rotationMatch(this, &matter, parameters->distanceDifference);
//    }else{
//        return (parameters->distanceDifference) > perAtomNorm(matter);
//    }
//}

bool Matter::compare(const Matter *matter, bool indistinguishable) {
    if (nAtoms != matter->numberOfAtoms())
        return false;
    if (parameters->checkRotation && indistinguishable) {
        return helper_functions::sortedR(this, matter, parameters->distanceDifference);
    } else if (indistinguishable) {
        if (this->numberOfFixedAtoms() == 0 and parameters->removeTranslation)
            helper_functions::translationRemove(this, matter);
        return helper_functions::identical(this, matter, parameters->distanceDifference);
    } else if (parameters->checkRotation) {
        return helper_functions::rotationMatch(this, matter, parameters->distanceDifference);
    } else {
        if (this->numberOfFixedAtoms() == 0 and parameters->removeTranslation)
            helper_functions::translationRemove(this, matter);
        return (parameters->distanceDifference) > perAtomNorm(*matter);
    }
}

// bool Matter::operator!=(const Matter& matter) {
//    return !operator==(matter);
//}

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
    AtomMatrix pbcMatrix = pbc(AtomMatrix::Map(diffVector.data(), diffVector.size() / 3, 3));
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

        //        if (potential)
        //        {
        //            delete potential;
        //            potential = NULL;
        //        }

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

bool Matter::relax(
    bool quiet, bool writeMovie, bool checkpoint, string prefixMovie, string prefixCheckpoint) {
    MatterObjectiveFunction objf(this, parameters);
    Optimizer *optimizer = Optimizer::getOptimizer(&objf, parameters);

    ostringstream min;
    min << prefixMovie;
    if (writeMovie) {
        matter2con(min.str(), false);
    }

    int iteration = 0;
    const char *forceLabel = parameters->optConvergenceMetricLabel.c_str();
    if (!quiet) {
        log("%s %10s  %14s  %18s  %13s\n",
            LOG_PREFIX.c_str(),
            "Iter",
            "Step size",
            forceLabel,
            "Energy");
        log("%s %10i  %14.5e  %18.5e  %13.5f\n",
            LOG_PREFIX.c_str(),
            iteration,
            0.0,
            objf.getConvergence(),
            getPotentialEnergy());
    }

    while (!objf.isConverged() && iteration < parameters->optMaxIterations) {

        AtomMatrix pos = getPositions();

        optimizer->step(parameters->optMaxMove);
        iteration++;

        double stepSize = helper_functions::maxAtomMotion(pbc(getPositions() - pos));

        if (!quiet) {
            log("%s %10i  %14.5e  %18.5e  %13.5f\n",
                LOG_PREFIX.c_str(),
                iteration,
                stepSize,
                objf.getConvergence(),
                getPotentialEnergy());
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
            log("%s %10i  %14.5e  %18.5e  %13.5f\n",
                LOG_PREFIX.c_str(),
                iteration,
                0.0,
                objf.getConvergence(),
                getPotentialEnergy());
        }
    }
    //    bool converged = optimizer->run(parameters->optMaxIterations, parameters->optMaxMove);
    delete optimizer;
    return objf.isConverged();
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

void Matter::setBiasPotential(BondBoost *bondBoost) { biasPotential = bondBoost; }

void Matter::setBiasForces(const AtomMatrix bf) { biasForces = bf.array() * getFree().array(); }
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

VectorXd Matter::getForcesV() { return VectorXd::Map(getForces().data(), 3 * numberOfAtoms()); }

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
double Matter::distance(long int index1, long int index2) const {
    return pbc(positions.row(index1) - positions.row(index2)).norm();
}

// return projected distance between the atoms with index1 and index2 on asix (0-x,1-y,2-z)
double Matter::pdistance(long int index1, long int index2, int axis) const {
    Matrix<double, 1, 3> ret;
    ret.setZero();
    ret(0, axis) = positions(index1, axis) - positions(index2, axis);
    ret = pbc(ret);
    return ret(0, axis);
}

// return the distance atom with index has moved between the current Matter object and the Matter
// object passed as argument
double Matter::distance(const Matter &matter, long int index) const {
    return pbc(positions.row(index) - matter.getPositions().row(index)).norm();
}

double Matter::getMass(long int indexAtom) const { return (masses[indexAtom]); }

void Matter::setMass(long int indexAtom, double mass) { masses[indexAtom] = mass; }

void Matter::setMasses(VectorXd massesIn) {
    for (int i = 0; i < nAtoms; i++) {
        masses[i] = massesIn[i];
    }
}

long int Matter::getAtomicNr(long int indexAtom) const { return (atomicNrs[indexAtom]); }

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
//}

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

double Matter::getMechanicalEnergy() { return getPotentialEnergy() + getKineticEnergy(); }

long int Matter::numberOfFreeAtoms() const { return nAtoms - isFixed.sum(); }

long int Matter::numberOfFixedAtoms() const { return isFixed.sum(); }

long int Matter::getForceCalls() const { return (forceCalls); }

void Matter::resetForceCalls() {
    forceCalls = 0;
    return;
}

// Print atomic coordinate to a .xyz file
void Matter::matter2xyz(std::string filename, bool append /*Append if file already exist*/) {
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
        applyPeriodicBoundary(); // Transform the coordinate to use the minimum image convention.
    }

    for (i = 0; i < numberOfAtoms(); i++) {
        fprintf(file,
                "%s\t%11.6f\t%11.6f\t%11.6f\n",
                atomicNumber2symbol(getAtomicNr(i)).data(),
                getPosition(i, 0),
                getPosition(i, 1),
                getPosition(i, 2));
    }
    fclose(file);
}

// Print atomic coordinates to a .con file
bool Matter::matter2con(std::string filename, bool append) {
    bool state;
    fs::path filePath = filename;
    std::ofstream file;
    if (not(filePath.extension() == ".con")) {
        filePath.replace_extension(".con");
    }
    if (append) {
        file.open(filePath, std::ios::app);
    } else {
        file.open(filePath);
    }
    state = matter2con(file);
    file.close();
    return (state);
}

bool Matter::matter2con(std::ofstream &file) {
    long int i;
    int j;
    long int Nfix = 0; // Nfix to store the number of fixed atoms
    int Ncomponent
        = 0;         // used to store the number of components (eg water: two components H and O)
    int first[MAXC]; // to store the position of the first atom of each component plus at the end
                     // the total number of atoms
    double mass[MAXC];
    long atomicNrs[MAXC];
    first[0] = 0;

    if (usePeriodicBoundaries) {
        applyPeriodicBoundary(); // Transform the coordinate to use the minimum image convention.
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
            Nfix++;                           // count the number of fixed atoms
        if (getAtomicNr(i) != atomicNrs[j]) { // check if there is a second component
            j++;
            if (j >= MAXC) {
                std::cerr << "Does not support more than " << MAXC
                          << " components and the atoms must be ordered by component.\n";
                return false;
            };
            mass[j] = getMass(i);
            atomicNrs[j] = getAtomicNr(i);
            first[j] = i;
        }
    }
    first[j + 1] = numberOfAtoms();
    Ncomponent = j + 1;

    file << headerCon1 << std::endl;
    file << headerCon2 << std::endl;
    double lengths[3];
    lengths[0] = cell.row(0).norm();
    lengths[1] = cell.row(1).norm();
    lengths[2] = cell.row(2).norm();
    file << fmt::sprintf("%f\t%f\t%f\n", lengths[0], lengths[1], lengths[2]);
    double angles[3];
    angles[0] = acos(cell.row(0).dot(cell.row(1)) / lengths[0] / lengths[1]) * 180 / M_PI;
    angles[1] = acos(cell.row(0).dot(cell.row(2)) / lengths[0] / lengths[2]) * 180 / M_PI;
    angles[2] = acos(cell.row(1).dot(cell.row(2)) / lengths[1] / lengths[2]) * 180 / M_PI;
    file << fmt::sprintf("%f\t%f\t%f\n", angles[0], angles[1], angles[2]);
    file << headerCon5 << std::endl;
    file << headerCon6 << std::endl;
    file << fmt::sprintf("%d\n", Ncomponent);

    for (j = 0; j < Ncomponent; j++) {
        file << fmt::sprintf("%d ", first[j + 1] - first[j]);
    }
    file << std::endl;
    for (j = 0; j < Ncomponent; j++) {
        file << fmt::format("{} ", mass[j]);
    }
    file << std::endl;
    for (j = 0; j < Ncomponent; j++) {
        file << fmt::sprintf("%s\n", atomicNumber2symbol(atomicNrs[j]));
        file << fmt::sprintf("Coordinates of Component %d\n", j + 1);
        for (i = first[j]; i < first[j + 1]; i++) {
            file << fmt::sprintf("%22.17f %22.17f %22.17f %d %4ld\n",
                                 getPosition(i, 0),
                                 getPosition(i, 1),
                                 getPosition(i, 2),
                                 getFixed(i),
                                 i + 1);
        }
    }
    return true;
}

// Load atomic coordinates from a .con file
bool Matter::con2matter(std::string filename) {
    bool state{false};
    fs::path filePath = filename;
    std::ifstream file{filePath};
    if (not(filePath.extension() == ".con")) {
        filePath.replace_extension(".con");
    }
    if (not file.is_open()) {
        std::cerr << fmt::format("File {} was not found.\n", filePath.string());
    } else {
        state = con2matter(file);
        file.close();
    }
    return (state);
}

bool Matter::con2matter(std::ifstream &file) {
    std::string line; // Temporary string
    size_t Ncomponent{0};
    std::getline(file, headerCon1);

    //    if (strchr(headerCon1,'\r')) {
    //        /* Files created on Windows or on Mac with Excell have carriage returns (\r) instead
    //        of or along with the new line charater (\n). C recognises only the \n as the end of
    //        line. */ cerr << "A carriage return ('\\r') has been detected. To work correctly, new
    //        lines should be indicated by the new line character (\\n)."; return false; // return
    //        false for error
    //    }

    std::getline(file, headerCon2);

    // The third line contains the length of the periodic cell
    std::getline(file, line);
    std::vector<double> lengths = helper_functions::get_val_from_string<double>(line, 3);

    // The fourth line contains the angles of the cell vectors
    std::getline(file, line);
    std::vector<double> angles = helper_functions::get_val_from_string<double>(line, 3);
    // TODO: Raise if nelements is provided but does not match the return.

    // Matter::cell assignment
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
    /*Matter::*/ cellInverse = cell.inverse();

    std::getline(file, headerCon5);
    std::getline(file, headerCon6);

    // Number of components or different types of atoms  (e.g. water: two components H and O)
    std::getline(file, line);
    if (line.empty()) {
        std::cout
            << "The number of components seems to be missing. One component is assumed instead\n"s;
        Ncomponent = 1;
    } else {
        // Guaranteed to be 1 element, so this is fine
        auto vtest = helper_functions::get_val_from_string<double>(line, 1)[0];
        if (vtest < 0) {
            std::cerr << "con2atoms does not support negative counts for the atoms.\n";
            return false;
        } else if (vtest == 0) {
            std::cout
                << "The number of components cannot be read. One component is assumed instead\n"s;
            Ncomponent = 1;
        } else {
            Ncomponent = static_cast<size_t>(vtest);
        }
    }

    // Use a vector of pairs to keep track of the components and their count
    std::vector<std::pair<size_t, size_t>> ncomp_count;

    // Now we want to know the number of atom of each type.
    // e.g with H2O, two Hydrogen atoms and one Oxygen atom
    std::getline(file, line);
    auto ncomps = helper_functions::get_val_from_string<size_t>(line);
    if (not(ncomps.size() == Ncomponent)) {
        std::cerr << "input con file does not list the number of each component";
        std::cerr << fmt::format("found {} components, {} types\n", Ncomponent, ncomps.size());
        return false;
    }
    for (size_t idx{1}; auto compnum : ncomps) { // Additional initializations C++20
        ncomp_count.push_back(std::make_pair(idx, compnum));
        ++idx;
    }
    // Asign nAtoms
    /*Matter::*/ nAtoms = std::accumulate(ncomps.begin(), ncomps.end(), 0);

    // Now we can initialize the fields in Matter
    positions = AtomMatrix::Constant(nAtoms, 3, 0);  // TEMP
    velocities = AtomMatrix::Constant(nAtoms, 3, 0); // TEMP
    forces = AtomMatrix::Constant(nAtoms, 3, 0);     // TEMP
    masses = VectorXd::Constant(nAtoms, 0);          // TEMP
    atomicNrs = VectorXi::Constant(nAtoms, 0);       // TEMP
    isFixed = VectorXi::Constant(nAtoms, false);     // TEMP

    // Get unique masses
    // These will be used to eventually populate the total masses
    std::getline(file, line);
    auto masses_unique = helper_functions::get_val_from_string<double>(line);
    if (not(masses_unique.size() == Ncomponent)) {
        std::cerr << "input con file does not list the masses of each component";
        std::cerr << fmt::format(
            "found {} components, {} types\n", Ncomponent, masses_unique.size());
        return false;
    }

    // Get atomic number and continue
    // The idea is we have both the number of components and the component index
    // in ncomp_count so we can use this information to validate
    std::vector<double> atomic_nrs;
    AtomMatrix _pos = AtomMatrix::Constant(nAtoms, 3, 0); // TEMP
    std::vector<double> masses_std;                       // Will be mapped to masses
    for (size_t idx{0}, compid{0}; auto comp : ncomp_count) {
        compid = comp.first - 1;
        // Element {BLAH}
        std::getline(file, line);
        auto sym = symbol2atomicNumber(line);
        // Coordinates of component BLAH
        std::getline(file, line); // Skip line with component number
        for (size_t a{0}; a < comp.second; a++) {
            // Now parse:
            // x y z is_fixed atom_index
            std::getline(file, line);
            auto tmp = helper_functions::get_val_from_string<double>(line, 5);
            for (size_t b{0}; b < 3; b++) {
                _pos(idx, b) = tmp[b]; // x y z
            }
            /*Matter::*/ setFixed(idx, static_cast<bool>(tmp[3])); // is_fixed
            if (( tmp[4] != idx + 1 ) and (tmp[4] != idx)) {
                std::cerr << fmt::format(
                    "Atoms in con file are not numbered from 1...N:\n expected:\t{} or {} got:\t{}",
                    idx + 1, idx,
                    tmp[4]);
                return false;
            } else {
                /*Matter::*/ setMass(idx, masses_unique[compid]);
                /*Matter::*/ setAtomicNr(idx, sym);
            }
            ++idx;
        }
    }
    /*Matter::*/ setPositions(_pos); // Finalize

    if (usePeriodicBoundaries) {
        applyPeriodicBoundary(); // Transform the coordinate to use the minimum image convention.
    }

    // potential_ = new Potential(parameters_);
    /*Matter::*/ recomputePotential = true;
    return true;
}

void Matter::computePotential() {
    if (recomputePotential) {
        if (!potential) {
            potential = Potential::getPotential(parameters);
        }

        // TODO: Handle the number of system images better
        forces = potential->force(nAtoms, positions, atomicNrs, &potentialEnergy, cell, 1);
        forceCalls = forceCalls + 1;
        recomputePotential = false;

        if (isFixed.sum() == 0 and parameters->removeNetForce) {
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

VectorXd Matter::getFreeV() const { return VectorXd::Map(getFree().data(), 3 * numberOfAtoms()); }

AtomMatrix Matter::getVelocities() const { return velocities.array() * getFree().array(); }

void Matter::setVelocities(const AtomMatrix v) { velocities = v.array() * getFree().array(); }

void Matter::setForces(const AtomMatrix f) { forces = f.array() * getFree().array(); }

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
    fs::path filePath = filename;
    std::ofstream file;
    if (not(filePath.extension() == ".convel")) {
        filePath.replace_extension(".convel");
    }
    file.open(filePath);
    state = matter2convel(file);
    file.close();
    return (state);
}

bool Matter::matter2convel(std::ofstream &file) {
    long int i;
    int j;
    long int Nfix = 0; // Nfix to store the number of fixed atoms
    int Ncomponent
        = 0;         // used to store the number of components (eg water: two components H and O)
    int first[MAXC]; // to store the position of the first atom of each component plus at the end
                     // the total number of atoms
    double mass[MAXC];
    long atomicNrs[MAXC];
    first[0] = 0;

    if (usePeriodicBoundaries) {
        applyPeriodicBoundary(); // Transform the coordinate to use the minimum image convention.
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
            Nfix++;                           // count the number of fixed atoms
        if (getAtomicNr(i) != atomicNrs[j]) { // check if there is a second component
            j++;
            if (j >= MAXC) {
                std::cerr << "Does not support more than " << MAXC
                          << " components and the atoms must be ordered by component.\n";
                return false;
            }
            mass[j] = getMass(i);
            atomicNrs[j] = getAtomicNr(i);
            first[j] = i;
        }
    }
    first[j + 1] = numberOfAtoms();
    Ncomponent = j + 1;

    file << headerCon1 << std::endl;
    file << headerCon2 << std::endl;
    double lengths[3];
    lengths[0] = cell.row(0).norm();
    lengths[1] = cell.row(1).norm();
    lengths[2] = cell.row(2).norm();
    file << fmt::sprintf("%f\t%f\t%f\n", lengths[0], lengths[1], lengths[2]);
    double angles[3];
    angles[0] = acos(cell.row(0).dot(cell.row(1)) / lengths[0] / lengths[1]) * 180 / M_PI;
    angles[1] = acos(cell.row(0).dot(cell.row(2)) / lengths[0] / lengths[2]) * 180 / M_PI;
    angles[2] = acos(cell.row(1).dot(cell.row(2)) / lengths[1] / lengths[2]) * 180 / M_PI;
    file << fmt::sprintf("%f\t%f\t%f\n", angles[0], angles[1], angles[2]);
    file << headerCon5 << std::endl;
    file << headerCon6 << std::endl;

    file << fmt::sprintf("%d\n", Ncomponent);
    for (j = 0; j < Ncomponent; j++) {
        file << fmt::sprintf("%d ", first[j + 1] - first[j]);
    }
    file << fmt::sprintf("\n");
    for (j = 0; j < Ncomponent; j++) {
        // mass[j]/=G_PER_MOL; // GH: I don't understand why we need to convert the mass units
        file << fmt::sprintf("%f ", mass[j]);
    }
    file << fmt::sprintf("\n");
    for (j = 0; j < Ncomponent; j++) {
        file << fmt::sprintf("%s\n", atomicNumber2symbol(atomicNrs[j]).data());
        file << fmt::sprintf("Coordinates of Component %d\n", j + 1);
        for (i = first[j]; i < first[j + 1]; i++) {
            file << fmt::sprintf("%11.6f\t%11.6f\t%11.6f\t%d\t%ld\n",
                                 getPosition(i, 0),
                                 getPosition(i, 1),
                                 getPosition(i, 2),
                                 getFixed(i),
                                 i);
        }
    }
    file << fmt::sprintf("\n");
    for (j = 0; j < Ncomponent; j++) {
        file << fmt::sprintf("%s\n", atomicNumber2symbol(atomicNrs[j]).data());
        file << fmt::sprintf("Velocities of Component %d\n", j + 1);
        for (i = first[j]; i < first[j + 1]; i++) {
            file << fmt::sprintf("%11.6f\t%11.6f\t%11.6f\t%d\t%ld\n",
                                 velocities(i, 0),
                                 velocities(i, 1),
                                 velocities(i, 2),
                                 getFixed(i),
                                 i);
        }
    }
    return true;
}

// TODO: Test
bool Matter::convel2matter(std::string filename) {
    bool state{false};
    fs::path filePath = filename;
    std::ifstream file{filePath};
    if (not(filePath.extension() == ".convel")) {
        filePath.replace_extension(".convel");
    }
    if (not file.is_open()) {
        std::cerr << fmt::format("File {} was not found.\n", filePath.string());
    } else {
        state = convel2matter(file);
        file.close();
    }
    return (state);
}
bool Matter::convel2matter(std::ifstream &file) {
    std::string line; // Temporary string
    size_t Ncomponent{0};
    std::getline(file, headerCon1);

    std::getline(file, headerCon2);

    // The third line contains the length of the periodic cell
    std::getline(file, line);
    std::vector<double> lengths = helper_functions::get_val_from_string<double>(line, 3);

    // The fourth line contains the angles of the cell vectors
    std::getline(file, headerCon4);
    std::vector<double> angles = helper_functions::get_val_from_string<double>(headerCon4, 3);

    // Matter::cell assignment
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
    /*Matter::*/ cellInverse = cell.inverse();

    std::getline(file, headerCon5);
    std::getline(file, headerCon6);

    // Number of components or different types of atoms  (e.g. water: two components H and O)
    std::getline(file, line);
    if (line.empty()) {
        std::cout
            << "The number of components seems to be missing. One component is assumed instead\n"s;
        Ncomponent = 1;
    } else {
        // Guaranteed to be 1 element, so this is fine
        auto vtest = helper_functions::get_val_from_string<double>(line, 1)[0];
        if (vtest < 0) {
            std::cerr << "con2atoms does not support negative counts for the atoms.\n";
            return false;
        } else if (vtest == 0) {
            std::cout
                << "The number of components cannot be read. One component is assumed instead\n"s;
            Ncomponent = 1;
        } else {
            Ncomponent = static_cast<size_t>(vtest);
        }
    }

    // Use a vector of pairs to keep track of the components and their count
    std::vector<std::pair<size_t, size_t>> ncomp_count;

    // Now we want to know the number of atom of each type.
    // e.g with H2O, two Hydrogen atoms and one Oxygen atom
    std::getline(file, line);
    auto ncomps = helper_functions::get_val_from_string<size_t>(line);
    if (not(ncomps.size() == Ncomponent)) {
        std::cerr << "input con file does not list the number of each component";
        std::cerr << fmt::format("found {} components, {} types\n", Ncomponent, ncomps.size());
        return false;
    }
    for (size_t idx{1}; auto compnum : ncomps) { // Additional initializations C++20
        ncomp_count.push_back(std::make_pair(idx, compnum));
        ++idx;
    }
    // Asign nAtoms
    /*Matter::*/ nAtoms = std::accumulate(ncomps.begin(), ncomps.end(), 0);

    // Now we can initialize the fields in Matter
    positions = AtomMatrix::Constant(nAtoms, 3, 0);  // TEMP
    velocities = AtomMatrix::Constant(nAtoms, 3, 0); // TEMP
    forces = AtomMatrix::Constant(nAtoms, 3, 0);     // TEMP
    masses = VectorXd::Constant(nAtoms, 0);          // TEMP
    atomicNrs = VectorXi::Constant(nAtoms, 0);       // TEMP
    isFixed = VectorXi::Constant(nAtoms, false);     // TEMP

    // Get unique masses
    // These will be used to eventually populate the total masses
    std::getline(file, line);
    auto masses_unique = helper_functions::get_val_from_string<double>(line);
    if (not(masses_unique.size() == Ncomponent)) {
        std::cerr << "input con file does not list the masses of each component";
        std::cerr << fmt::format(
            "found {} components, {} types\n", Ncomponent, masses_unique.size());
        return false;
    }

    // Get atomic number and continue
    // The idea is we have both the number of components and the component index
    // in ncomp_count so we can use this information to validate
    std::vector<double> atomic_nrs;
    AtomMatrix _pos = AtomMatrix::Constant(nAtoms, 3, 0); // TEMP
    std::vector<double> masses_std;                       // Will be mapped to masses
    for (size_t idx{0}, compid{0}; auto comp : ncomp_count) {
        compid = comp.first - 1;
        // Element {BLAH}
        std::getline(file, line);
        auto sym = symbol2atomicNumber(line);
        // Coordinates of component BLAH
        std::getline(file, line); // Skip line with component number
        for (size_t a{0}; a < comp.second; a++) {
            // Now parse:
            // x y z is_fixed atom_index
            std::getline(file, line);
            auto tmp = helper_functions::get_val_from_string<double>(line, 5);
            for (size_t b{0}; b < 3; b++) {
                _pos(idx, b) = tmp[b]; // x y z
            }
            /*Matter::*/ setFixed(idx, static_cast<bool>(tmp[3])); // is_fixed
            if (( tmp[4] != idx + 1 ) and (tmp[4] != idx)) {
                std::cerr << fmt::format(
                    "Atoms in con file are not numbered from 1...N:\n expected:\t{} or {} got:\t{}",
                    idx + 1, idx,
                    tmp[4]);
                return false;
            } else {
                /*Matter::*/ setMass(idx, masses_unique[compid]);
                /*Matter::*/ setAtomicNr(idx, sym);
            }
            ++idx;
        }
    }
    /*Matter::*/ setPositions(_pos); // Finalize
    // Now reading velocities
    for (size_t idx{0}, compid{0}; auto comp : ncomp_count) {
        // Element {BLAH}
        std::getline(file, line);
        std::getline(file, line); // Skip line with component number
        for (size_t a{0}; a < comp.second; a++) {
            // Now parse:
            // x y z is_fixed atom_index
            std::getline(file, line);
            auto tmp = helper_functions::get_val_from_string<double>(line, 3);
            setVelocity(idx, 0, tmp[0]); // x
            setVelocity(idx, 1, tmp[1]); // y
            setVelocity(idx, 2, tmp[2]); // z
            ++idx;
        }
    }

    if (usePeriodicBoundaries) {
        /*Matter::*/ applyPeriodicBoundary(); // Transform the coordinate to use the minimum image
                                              // convention.
    }
    //    potential_ = new Potential(parameters_);
    return true;
}

// MatterObjective
double MatterObjectiveFunction::getEnergy() { return matter->getPotentialEnergy(); }
double MatterObjectiveFunction::getConvergence() {
    if (parameters->optConvergenceMetric == "norm") {
        return matter->getForcesFreeV().norm();
    } else if (parameters->optConvergenceMetric == "max_atom") {
        return matter->maxForce();
    } else if (parameters->optConvergenceMetric == "max_component") {
        return matter->getForces().maxCoeff();
    } else {
        log("%s Unknown opt_convergence_metric: %s\n",
            LOG_PREFIX.c_str(),
            parameters->optConvergenceMetric.c_str());
        exit(1);
    }
}
bool MatterObjectiveFunction::isConverged() {
    return getConvergence() < parameters->optConvergedForce;
}
int MatterObjectiveFunction::degreesOfFreedom() { return 3 * matter->numberOfFreeAtoms(); }
VectorXd MatterObjectiveFunction::getPositions() { return matter->getPositionsFreeV(); }
VectorXd MatterObjectiveFunction::getGradient(bool fdstep) {
    return -matter->getForcesFreeV();
}
VectorXd MatterObjectiveFunction::difference(VectorXd a, VectorXd b) {
    return matter->pbcV(a - b);
}
void MatterObjectiveFunction::setPositions(VectorXd x) { matter->setPositionsFreeV(x); }

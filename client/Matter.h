#ifndef MATTER_H
#define MATTER_H

#include "Eigen.h"
#include "Parameters.h"
#include "Potential.h"
#include "SurrogatePotential.h"
#include <memory>
#include <spdlog/sinks/basic_file_sink.h>

// This is a forward declaration of BondBoost to avoid a circular dependency.
class BondBoost;

#include <string>

/* Data describing an atomic structure. This class has been devised to handle
 * information about an atomic structure such as positions, velocities, masses,
 * etc. It also allow to associate a forcefield for the structure through a
 * pointer to function (potential()). The class can read and save data to a .con
 * file (atom2con() and con2atom()). It can also save to a .xyz file
 * (atom2xyz()).*/

class Matter {
public:
  ~Matter() = default;
  Matter(std::shared_ptr<Potential> pot, std::shared_ptr<Parameters> params)
      : potential{pot},
        usePeriodicBoundaries{true},
        recomputePotential{true},
        forceCalls{0},
        parameters{params},
        nAtoms{0},
        positions{Eigen::MatrixXd::Zero(0, 3)},
        velocities{Eigen::MatrixXd::Zero(0, 3)},
        forces{Eigen::MatrixXd::Zero(0, 3)},
        biasForces{Eigen::MatrixXd::Zero(0, 3)},
        biasPotential{nullptr},
        masses{Eigen::VectorXd::Zero(0)},
        atomicNrs{Eigen::VectorXi::Zero(0)},
        isFixed{Eigen::VectorXi::Zero(0)},
        cell{Eigen::Matrix3d::Zero()},
        cellInverse{Eigen::Matrix3d::Zero()},
        energyVariance{0.0},
        potentialEnergy{0.0} {

    m_log = spdlog::get("combi");
  } // the number of atoms shall be set later
    // using resize()
  // TODO: This is a placeholder, it delegates to the standard constructor
  // Matter(std::shared_ptr<Parameters> parameters)
  //     : Matter(parameters.get()) {
  // } // the number of atoms shall be set later using resize()
  // Matter(Parameters *parameters,
  //        long int nAtoms);      // prepare the object for use with nAtoms
  //        atoms
  Matter(const Matter &matter);                  // create a copy of matter
  const Matter &operator=(const Matter &matter); // copy the matter object
  // bool operator==(const Matter& matter); // true if differences in positions
  // are below differenceDistance bool operator!=(const Matter& matter); //
  // inverse of ==
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
  void setCell(Matrix3d newCell);
  double getPosition(long int atom, int axis)
      const; // return the position of an atom along one of the axis
  void setPosition(
      long int atom, int axis,
      double position); // set the position of atom along axis to position
  void setVelocity(
      long int atom, int axis,
      double velocity); // set the velocity of atom along axis to velocity
  bool relax(bool quiet = false, bool writeMovie = false,
             bool checkpoint = false, string prefixMovie = string(),
             string prefixCheckpoint = string());

  AtomMatrix pbc(AtomMatrix diff) const;
  VectorXd pbcV(VectorXd diff) const;

  size_t getPotentialCalls() const;
  AtomMatrix getPositions() const; // return coordinates of atoms in array pos
  VectorXd getPositionsV() const;
  AtomMatrix
  getPositionsFree() const; // return coordinates of free atoms in array pos
  VectorXd getPositionsFreeV() const;
  void
  setPositions(const AtomMatrix pos); // update Matter with the new positions of
                                      // the free atoms given in array pos
  void setPositionsV(const VectorXd pos);
  void setPositionsFree(
      const AtomMatrix pos); // update Matter with the new positions of the free
                             // atoms given in array pos
  void setPositionsFreeV(const VectorXd pos);

  AtomMatrix getVelocities() const;
  void setVelocities(const AtomMatrix v);
  void setBiasForces(const AtomMatrix bf);
  void setBiasPotential(BondBoost *bondBoost);
  void setForces(const AtomMatrix f);
  AtomMatrix getAccelerations();

  AtomMatrix getForces(); // return forces applied on all atoms in array force
  AtomMatrix getBiasForces();
  VectorXd getForcesV();
  AtomMatrix getForcesFree();
  VectorXd getForcesFreeV();

  double getMass(long int atom) const; // return the mass of the atom specified
  void setMass(long int atom, double mass); // set the mass of an atom
  void setMasses(VectorXd massesIn);        // set the mass of an atom
  long getAtomicNr(
      long int atom) const; // return the atomic number of the atom specified
  void setAtomicNr(long int atom,
                   long atomicNr);          // set the atomic number of an atom
  VectorXi getAtomicNrs() const;            // Get the vector of atomic numbers
  VectorXi getAtomicNrsFree() const;        // Get the vector of atomic numbers
  void setAtomicNrs(const VectorXi atmnrs); // Get the vector of atomic numbers

  int getFixed(long int atom)
      const; // return true if the atom is fixed, false if it is movable
  void setFixed(long int atom,
                int isFixed); // set the atom to fixed (true) or movable (false)
  // void setPotentialEnergy(double);
  double getEnergyVariance();
  // Eigen::VectorXd getForceVariance();
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

  long
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
  VectorXd getFreeV() const;
  Eigen::Matrix<double, Eigen::Dynamic, 1> getMasses() const;

private:
  shared_ptr<spdlog::logger> m_log;
  std::shared_ptr<Potential>
      potential; // pointer to function calculating the energy and forces
  bool usePeriodicBoundaries; // boolean telling periodic boundaries are used
  mutable bool recomputePotential; // boolean indicating if the potential energy
                                   // and forces need to be recalculated
  mutable long
      forceCalls; // keep track of how many force calls have been performed

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
  std::shared_ptr<Parameters> parameters;
  long nAtoms;
  AtomMatrix positions;
  AtomMatrix velocities;
  AtomMatrix forces;
  AtomMatrix biasForces;
  BondBoost *biasPotential;
  VectorXd masses;
  VectorXi atomicNrs;
  VectorXi isFixed; // array of bool, false for movable atom, true for fixed
  Matrix3d cell;
  Matrix3d cellInverse;
  mutable double energyVariance;
  mutable double potentialEnergy;
};

#endif

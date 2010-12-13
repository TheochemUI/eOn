//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef MATTER_H
#define MATTER_H

#include "Parameters.h"
#include "Potentials.h"

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

/* This structure is used as a inherited by Matter. It contains data about an atomic structure. It is used when a fast and direct access to the private data members is required in Matter. It should not be used alone.  Use carefully. Avoid it and use Matter instead if possible. */

struct MatterPrivateData {
    Parameters *parameters;
 
    long int nAtoms; // number of atoms
    Matrix<double, Eigen::Dynamic, 3> positions; // positions
    Matrix<double, Eigen::Dynamic, 3> velocities; // velocities
    Matrix<double, Eigen::Dynamic, 3> forces; // forces
    Matrix<double, Eigen::Dynamic, 1> masses; // masses
    Matrix<int, Eigen::Dynamic, 1> atomicNrs; // atomic numbers
    Matrix<int, Eigen::Dynamic, 1> isFixed; // array of bool, false for movable atom, true for fixed
    Matrix<double, 3, 3> cellBoundaries; // boundaries of the cell
    mutable double potentialEnergy; // potential energy
};

/* Data describing an atomic structure. This class has been devised to handle information about an atomic structure such as positions, velocities, masses, etc. It also allow to associate a forcefield for the structure through a pointer to function (potential()). The class can read and save data to a .con file (atom2con() and con2atom()). It can also save to a .xyz file (atom2xyz()).*/

class Matter : private MatterPrivateData {
public:
    Matter(Parameters *parameters); // the number of atoms shall be set later using resize()
    Matter(Parameters *parameters, long int nAtoms); // prepare the object for use with nAtoms atoms
    Matter(const Matter& matter); // create a copy of matter
    ~Matter(); // Destructor
    const Matter& operator=(const Matter& matter); // copy the matter object
    bool operator==(const Matter& matter); // true if differences in positions are below differenceDistance
    double distanceTo(const Matter& matter); // the distance to the given matter object
    double perAtomNorm(const Matter& matter); // the maximum distance between two atoms in the Matter objects
    void setPotential(); // set potential function to use
    void resize(long int nAtoms); // set or reset the number of atoms
    long int numberOfAtoms() const; // return the number of atoms
    Vector3d getBoundary(int axis) const; // return the length of the periodic cell for the axis specified
    double getBoundary(int axis1, int axis2) const; // return the length of the periodic cell for the axis specified
    void setBoundary(int axis, Vector3d); // set the length of the periodic cell for the axis specified
    void setBoundary(int axis1, int axis2, double val); // set the length of the periodic cell for the axis specified
    void activatePeriodicBoundaries(); // activate the periodic boundary conditions
    void deactivatePeriodicBoundaries(); // deactivate periodic boundary conditions
    double getPosition(long int atom, int axis) const; // return the position of an atom along one of the axis
    void setPosition(long int atom, int axis, double position); // set the position of atom along axis to position

    Matrix<double, Eigen::Dynamic, 3> pbc(Matrix<double, Eigen::Dynamic, 3> diff) const;
    
    Matrix<double, Eigen::Dynamic, 3> getPositions() const; // return coordinates of free atoms in array pos
    void setPositions(const Matrix<double, Eigen::Dynamic, 3> pos); // update Matter with the new positions of the free atoms given in array pos
    
    Matrix<double, Eigen::Dynamic, 3> getVelocities() const; 
    void setVelocities(const Matrix<double, Eigen::Dynamic, 3> v); 
    void setForces(const Matrix<double, Eigen::Dynamic, 3> f);
    Matrix<double, Eigen::Dynamic, 3> getAccelerations(); 

    Matrix<double, Eigen::Dynamic, 3> getForces(); // return forces applied on all atoms in array force 

    double getMass(long int atom) const; // return the mass of the atom specified
    void setMass(long int atom, double mass); // set the mass of an atom
    long getAtomicNr(long int atom) const; // return the atomic number of the atom specified
    void setAtomicNr(long int atom, long atomicNr); // set the atomic number of an atom

    int getFixed(long int atom) const; // return true if the atom is fixed, false if it is movable
    void setFixed(long int atom, int isFixed); // set the atom to fixed (true) or movable (false)
    double getPotentialEnergy(); // return the potential energy
    double getKineticEnergy() const; // return the Kinetic energy
    double getMechanicalEnergy(); // return the mechanical energy (i.e. kinetic plus potential energy)

    double distance(long index1, long index2) const; // return the distance between two atoms in same configuration
    double pdistance(long index1, long index2, int axis) const;
    double distance(const Matter& matter, long index) const; // the distance between the same atom in two cofigurations

    long int numberOfFreeAtoms() const; // return the number of free (or movable) atoms

    long getForceCalls() const; // return how many force calls that have been performed
    void resetForceCalls(); // zeroing the value of force calls

    bool isItConverged(double convergeCriterion);
    double maxForce(void);

    bool con2matter(std::string filename); // read con file into Matter, return true if successful
    bool con2matter(FILE *file); // read con file and load data into Matter, return true if successful
    bool matter2con(std::string filename) const; // print con file from data in Class Matter
    bool matter2con(FILE *file) const; // print con file from data in Class Matter
    void matter2xyz(std::string filename, bool append=false) const; // print xyz file from data in Matter
    
    Matrix<double, Eigen::Dynamic, 3> getFree() const;
    Matrix<double, Eigen::Dynamic, 1> getMasses() const;

private:
    Potential *potential; // pointer to function calculating the energy and forces
    bool usePeriodicBoundaries; // boolean telling periodic boundaries are used
    mutable bool recomputePotential; // boolean indicating if the potential energy and forces need to be recalculated
    mutable long forceCalls; // keep track of how many force calls have been performed

    // CON file header information, which is not used in the eon code
    char headerCon1[STRING_SIZE];
    char headerCon2[STRING_SIZE];
    char headerCon4[STRING_SIZE];
    char headerCon5[STRING_SIZE];
    char headerCon6[STRING_SIZE];
    
    void computePotential();
    void initialiseDataMembers(Parameters *parameters);
    void clearMemory(); // clear dynamically allocated memory
    void applyPeriodicBoundary();
    void applyPeriodicBoundary(double & component, int axis);
    void applyPeriodicBoundary(Matrix<double, Eigen::Dynamic, 3> & diff);
};

#endif

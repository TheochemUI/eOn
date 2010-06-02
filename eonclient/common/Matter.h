/*
 *===============================================
 *  Matter.h
 *  eon2
 *-----------------------------------------------
 *  Created by Jean-Claude C. Berthet 2006.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
*/
#ifndef Matter_H
#define Matter_H

#include "system_unit.h" // unit converters
#include "Parameters.h"
#include "Potentials.h"
#include "Constants.h"
/**This structure is used as a inherited by Matter. It contains data about an atomic structure. It is used when a fast and direct access to the private data members is required in Matter. It should not be used alone.
@warning Use carefully. Avoid it and use Matter instead if possible.
@see Matter */
struct MatterPrivateData {
    Parameters *parameters_;
    
	long int nAtoms_; ///< Number of atoms
	double *positions_; ///< Array of 3*n elements containing the positions.
	/** Array of 3*n elements containing positions. This array is used when alternate coordinates need to be stored. This is the case when constrains are applied.
	@see	constraint() and UpdateCons()*/
	double *positionsBefore_;
	double *velocities_; ///< Array of 3*n elements containing the velocities.
	double *forces_; ///< Array of 3*n elements containing the forces.
	double *masses_; ///< Array of masses.
	long *atomicNrs_; ///< Array of atomic numbers.
	int *isFixed_; ///< Array of bool, false for movable atom, true for fixed.
	double cellBoundaries_[3]; ///< Boundaries of the cell.
	mutable double potentialEnergy_; ///< Potential energy.
};
/** Data describing an atomic structure. This class has been devised to handle information about an atomic structure such as positions, velocities, masses, etc. It also allow to associate a forcefield for the structure through a pointer to function (potential()). The class can read and save data to a .con file (atom2con() and con2atom()). It can also save to a .xyz file (atom2xyz()).*/
class Matter : private MatterPrivateData {
	/** Pointer to function to constraint the molecules. #constraint_ must point to a function which corrects the position and velocities of the atoms in order to get the system complying with constraints. Constraint algorithms usually requires 'reference' coordinates which means the coordinates before the modifications leading to a unconstrained coordinates were applied. Matter will provide these reference coordinates to constraint through the argument @em ref.
	@see	The function is called by UpdateCons(), UpdateAcc() and UpdateForce().
	@param[in]		nAtoms				Number of atoms(including fixed atoms)
	@param[in]		positionsBefore		Positions before the unconstrained operations
	@param[in,out]	positions				Coordinates to correct
	@param[in,out]	velocities				Velocities to correct
	@param[in,out]	cellBoundaries		Array of length 3 containing the parameters of the periodic cell.*/
	typedef void (*Constraints)(const long int nAtoms, const double positionsBefore[], double positions[], double velocities[], const double cellBoundaries[]);
public:
	Matter(Parameters *parameters);///< The number of atoms shall be set later using resize().
	Matter(Parameters *parameters, long int nAtoms);///< Prepare the object for use with nAtoms atoms.
	Matter(const Matter& matter);///< Create a copy of @em matter.
	~Matter();///< Destructor.
	const Matter& operator=(const Matter& matter);///< The object 'matter' is copied into this object.
    bool operator==(const Matter& matter);///< Returns true if all differences in positions are bellow getMaxDifferencePos, they are then considered equal.
    void setPotential();///< Set function to use as for potential. 
    void resize(long int nAtoms);///< Set or reset the number of atoms.
	long int numberOfAtoms() const;///< Return the number of atoms. 
	double getBoundary(int axis) const;///< Return the length of the periodic cell for the axis specified (axis=0, 1 or 2).
	void setBoundary(int axis, double length);///< Set the length of the periodic cell for the axis specified (axis=0, 1 or 2).
	void activatePeriodicBoundaries();/**< Activate the periodic boundary conditions and keep coordinates in the minimum image convention. 
		When the function is called, coordinates are recalculated to fit the minimum image convention (i.e. coordinates within[-cellBoundaries_[X]/2;+cellBoundaries_[X]/2], etc...). Subsequently, coordinates are also recalculated each time these are modified.
		@see desactivatePeriodicBoundaries()*/
	void desactivatePeriodicBoundaries();///< Desactivate periodic boundary conditions. Coordinates will no longer be recalculated to comply with the minimum image convention.
	void setConstraints(Constraints constraints);/**< Set the function to use constraints are needed.
		@param	constraints	set to zero to remove any constraints.*/
	double getPosition(long int atom, int axis) const;///< Return the position of an atom along one of the axis.
	void setPosition(long int atom, int axis, double position);///< Set or the position of @em atom along @em axis to @em position.

    void getPositions(double pos[]) const; //return coordinates of free atoms in array 'pos'
    void setPositions(const double pos[]); //Update Matter with the new positions of the free atoms given in array 'pos'
    void getForces(double forces[]) const; //return forces applied on all atoms in array 'force' 

    double getMass(long int atom) const;///< Return the mass of the atom specified.
	void setMass(long int atom, double mass);///< Set the mass of an atom.
    long getAtomicNr(long int atom) const;///< Return the atomic number of the atom specified.
    void setAtomicNr(long int atom, long atomicNr);///< Set the atomic number of an atom.

    int getFixed(long int atom) const;///< Return true if the atom is fixed, false if it is movable.
	void setFixed(long int atom, int isFixed);///< Set the atom to fixed (true) or movable (false)
	double potentialEnergy() const;///< Return the potential energy.
	double kineticEnergy() const;///< Return the Kinetic energy.
	double mechanicalEnergy() const;///< Return the mechanical energy (i.e. kinetic plus potential energy). 

    double distance(long index1, long index2) const; ///< Return distance between two atoms in same configuration.
    double distance(const Matter& matter, long index) const; ///< Return distance difference for the same atom in two cofigurations.

    long int numberOfFreeAtoms() const;///< Return the number of free (or movable) atoms.
	void getFreePositions(double positions[]) const;///< Store positions of free atoms in an array.
	void setFreePositions(const double positions[]);///< Change the positions of the free atoms to values provided in an array.
	bool getFreeVelocities(double velocities[]) const;/**< Stored the velocities of the free atoms in an array. 
	@return	Returns true if successfull.*/
	void setFreeVelocities(const double velocities[]);///< Change the velocities of the free atoms to values provided in an array.
	void getFreeAccelerations(double accelerations[]) const;///< Return the accelarations of the free atoms in an array.
	void getFreeForces(double forces[]) const;///< Return forces applied on free atoms in an array.
	void getFreeMasses(double masses[]) const;///< Return the masses of free atoms in an array.
    
    long getForceCalls() const;///< Return how many force calls that have been performed.
    void resetForceCalls();///< Zeroing the value of force calls.

	void updateForces(double positions[], double velocities[], double forces[]);
	void updateAccelerations(double positions[], double velocities[], double accelerations[]);

	bool con2matter(std::string filename);/**< Read '.con' file and load data into Matter. Support up to ten components (ten types of atoms). 
	@return Returns true if successfull.*/
    bool con2matter(FILE *file);/**< Read '.con' file and load data into Matter. Support up to ten components (ten types of atoms). 
    @return Returns true if successfull.*/    
	bool matter2con(std::string filename) const;///< Print @em .con file from data in Class Matter. 
    bool matter2con(FILE *file) const;///< Print @em .con file from data in Class Matter. 
	void matter2xyz(std::string filename, bool append=false /*Append if file already exists*/) const;/**< Print @em .xyz file based on data stored in Class Matter. 
    @param append	Matters if file @filename already exists. When true, append new data to the existing file.*/

private:
    Potentials *potential_;/// Pointer to function calculating the energy and forces.
	Constraints constraints_;/// Pointer to function constraining some of the coordinates of the structure.
	bool usePeriodicBoundaries_;///< Boolean telling periodic boundaries are used.
	mutable bool computePotential_; ///< Boolean telling if the potential energy and forces need to be recalculated.
    mutable long forceCalls_; ///< Keeping track of how many force calls that have been performed.

    char headerCon1_[100];///< To contain headerline 1 from CON file that is not used in this code.
    char headerCon2_[100];///< To contain headerline 2 from CON file that is not used in this code.
    char headerCon4_[100];///< To contain headerline 4 from CON file that is not used in this code.
    char headerCon5_[100];///< To contain headerline 5 from CON file that is not used in this code.
    char headerCon6_[100];///< To contain headerline 6 from CON file that is not used in this code.
    
	void computePotential() const;
	void initialiseDataMembers(Parameters *parameters);
	void clearMemory(); ///< Clear all dynamically allocated memory.
	void applyPeriodicBoundary();
	void applyPeriodicBoundary(int axis);
	void applyPeriodicBoundary(long atom, int axis);
	void applyConstraints();
};

#endif

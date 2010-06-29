/*
 *===============================================
 *  Potential.h
 *  eon2
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 10/4/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef POTENTIALS
#define POTENTIALS

#include "Constants.h"
#include "Parameters.h"
#include "PotentialsInterface.h"
#include "potentials/NewPotential/NewPotential.h"
//#include "VASP.h"
//#include "EDIP.h"
//#include "potentials/EMT/EffectiveMediumTheory.h"
//#include "LJBinary.h"
#include "potentials/Morse/Morse.h"
//#include "LJ.h"
//#include "Lenosky.h"
//#include "SW.h"
//#include "Tersoff.h"
//#include "LJBinary.h"
#include "potentials/Aluminum/Aluminum.h"


/** Class serving as a wrapper between the force calculator and the Atoms object. It is here it is decided which potential that is going to be used! Might seem redundant but it is nessesary in order to have a clean code in the Matter class*/
class Potentials{
public:
    Potentials(Parameters *parameters);///< Constructor
    
    ~Potentials();///< Destructor
    
    /** A similar function should be provided by the force calculator.
    @param[in]    nAtoms      The number of atoms.
    @param[in]    *positions  Pointer to the array containing the atoms positions, should have size 3*N.
    @param[in]    *atomicNrs  Pointer to the array containing the atomic numbers, should have size N.
    @param[out]   *forces     Pointer to the array where the forces should be stored, should have size 3*N.
    @param[out]   *energy     Pointer where the total energy should be stored.
    @param[in]    *box        Pointer to the array containing the size of the supercell, Box[0] is length in the x directions, Box[1] is length in the y directions, Box[2] is length in the z directions.*/
    void force(long nAtoms, const double *positions, const long *atomicNrs, double *forces, double *energy, const double *box);
    
private:
    PotentialsInterface *interface_;
    Parameters *parameters_;
};

#endif

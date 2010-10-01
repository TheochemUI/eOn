/*
 *===============================================
 *  PotentialsInterface.h
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 1/3/07.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef POTENTIALS_INTERFACE_H
#define POTENTIALS_INTERFACE_H

/** Defining the interface for the potential.*/
class PotentialsInterface{
public:
    virtual ~PotentialsInterface();
    void virtual initialize() = 0;
    void virtual cleanMemory() = 0;
    //One should implement the method cleanMemory instead a destructor 
    //as the desctructor call cause troubles when linking towards the potential libraries

    
    /* A similar function should be provided by the force calculator.
    @param[in]    nAtoms      The number of atoms.
    @param[in]    *positions  Pointer to the array containing the atoms positions,
                              should have size 3*N.
    @param[in]    *atomicNrs  Pointer to the array containing the atomic numbers,
                              should have size N.
    @param[out]   *forces     Pointer to the array where the forces should be stored,
                              should have size 3*N.
    @param[out]   *energy     Pointer where the total energy should be stored.
    @param[in]    *box        Pointer to the array containing the size of the
                              supercell, Box[0] is length in the x directions,
                              Box[1] is length in the y directions,
                              Box[2] is length in the z directions.
    */    
    void virtual force(long nAtoms, const double *positions, const long *atomicNrs, 
                       double *forces, double *energy, const double *box) = 0;
};
#endif

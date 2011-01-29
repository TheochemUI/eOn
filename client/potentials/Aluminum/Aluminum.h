//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
//-----------------------------------------------------------------------------------
#ifndef ALUMINUM_POTENTIAL
#define ALUMINUM_POTENTIAL

#include "../../Potential.h"

    /** External function implemented in Fortran. Calculate interactions between molecules of water using forcefield Aluminum.
    @param[in]	N           Number of atoms.
    @param[in]	R           Array to positions of the atoms in Angstrom.
    @param[out]	F           Array used to return the forces resulting from interactions between molecules. Forces are in eV/Angstrom.
    @param[out]	U           Pointer to energy in eV.
    @param[in]  bx, by, bz  Pointer to box dimensions in Angstrom.
    */
    
extern "C" 
{
    void force_(const long int *N, const double *R, double *F, double *U, const double* bx, const double* by, const double* bz);
}    
extern "C" 
{
    void potinit_();
}    

/** Aluminum potential.*/
class Aluminum : public Potential
{     
public:
    Aluminum(void);
	
    // To satify interface
    void initialize(void);    
    void cleanMemory(void);    
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
};
#endif


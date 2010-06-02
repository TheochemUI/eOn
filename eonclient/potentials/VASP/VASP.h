
#ifndef VASP_POTENTIAL
#define VASP_POTENTIAL

#include <iostream>
#include <cstdio>
#include "stdlib.h"

#include "system_unit.h" // unit converters
#include "PotentialsInterface.h"

/** VASP potential.*/
class VASP : public PotentialsInterface{
    
public:
// Functions
	// constructor
    VASP(void);
	
    // Just to satify interface
    void initialize() {};
    void cleanMemory(void);    
    void force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box);

private:
    void writePositionsToFile(long N, const double *R, long const *atomicNrs, const double *box);
    void readForcesFromFile(long N, double *F, double *U);
};
#endif


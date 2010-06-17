
#ifndef LENNARD_JONES_BINARY
#define LENNARD_JONES_BINARY

//#include <stdlib.h>

#include <math.h> 
#include <iostream>
#include "system_unit.h" // unit converters
#include "PotentialsInterface.h"

/** Lennard Jones Binary potential.*/
class LJBinary : public PotentialsInterface{

private:
//	Variables
    double u0;
    double cuttOffR;
    double psi;
    
    double cuttOffU;
    
public:
// Functions
	// constructor
    LJBinary(void);
	LJBinary(double r0Recieved, double u0Recieved, double psiRecieved);
	
    // Just to satify interface
    void initialize() {};
    void cleanMemory(void);    
    
    void force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box);
    void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);
    double pbc_mdr(const double *r1, const double *r2, double bx, double by, double bz); 
    void pbc_vdr(const double *r1, const double *r2, double dr[3], double bx, double by, double bz); 

};
#endif


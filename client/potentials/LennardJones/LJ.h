
#ifndef LENNARD_JONES
#define LENNARD_JONES

#include <math.h> 
#include <iostream>
#include "../../system_unit.h" // unit converters
#include "../../PotentialsInterface.h"

/** Lennard Jones potential.*/
class LJ : public PotentialsInterface{

private:
//	Variables
    double u0;
    double cuttOffR;
    double psi;
    
    double cuttOffU;
    
public:
// Functions
	// constructor
    LJ(void);
	LJ(double r0Recieved, double u0Recieved, double psiRecieved);
	
    // Just to satify interface
    void initialize() {};
    void cleanMemory(void);    
    
    void force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box);
    void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);
};
#endif


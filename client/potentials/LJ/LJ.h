//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef LENNARD_JONES
#define LENNARD_JONES

#include <math.h> 
#include <iostream>
//#include "../../system_unit.h" // unit converters
#include "../../Potential.h"

/** Lennard Jones potential.*/
class LJ : public Potential
{

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

        ~LJ();
	
        // Just to satisfy interface
        void initialize() {};
        void cleanMemory();    
        
        void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
        void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);
};
#endif


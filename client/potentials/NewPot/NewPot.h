//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef NEWPOT_INTERFACE
#define NEWPOT_INTERFACE

#include "../../Potential.h"

/** Template to use if user want to provide potential. */
class NewPot : public Potential{

private:
//	Variables
    double fake1;
    double fake2;
    
public:
// Functions
	// constructor and destructor
    NewPot(Parameters *p);
	
    // To satisfy interface
    void initialize(void);    
    void cleanMemory(void);    

    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
};
#endif


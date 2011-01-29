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

#ifndef bopfox_POTENTIAL
#define bopfox_POTENTIAL

#include "../../Potential.h"

class bopfox : public Potential
{

    public:
        bopfox(void);
		~bopfox();
        void initialize() {};
        void cleanMemory(void);    
        void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);


    private:
        void writeFOX(long N, const double *R, int const *atomicNrs, const double *box);
        void readFU(long N, double *F, double *U);
        

};

#endif


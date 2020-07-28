//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef PYAMFF_POTENTIAL
#define PYAMFF_POTENTIAL

#include "../../Potential.h"
#include <string>

class PyAMFF : public Potential
{

    public:
        PyAMFF(void);
		~PyAMFF();
		void initialize() {};
		void cleanMemory(void);
        void force(long N, const double *R, const int *atomicNrs, 
                   double *F, double *U, const double *box);
        std::string num2sym(int Z);

    private:
        void writePOSCAR(long N, const double *R, const int *atomicNrs,
                         const double *box);
        void readFU(long N, double *F, double *U);
};

#endif


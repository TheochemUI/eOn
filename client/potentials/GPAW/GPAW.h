//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef GPAW_POTENTIAL
#define GPAW_POTENTIAL

#include "../../Potential.h"

//Forward declaration of MPI::Intercom to avoid including the MPI header
//in this header file that is included by Potential.h. This avoids needing
//to use mpic++ for any other file in the program.
namespace MPI
{
    class Intercomm;
}

class GPAW : public Potential
{

    public:
        GPAW(void);
		~GPAW();
        void initialize(){};
        void cleanMemory(void);    
        void force(long N, const double *R, const int *atomicNrs, 
                   double *F, double *U, const double *box);

    private:
        static MPI::Intercomm gpawComm;
        static int instances;
        static bool firstRun;
};

#endif

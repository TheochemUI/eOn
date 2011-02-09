//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef VASP_POTENTIAL
#define VASP_POTENTIAL

#include "../../Potential.h"

class VASP : public Potential
{

    public:
        VASP(void);
		~VASP();
        void initialize() {};
        void cleanMemory(void);    
        void force(long N, const double *R, const int *atomicNrs, 
                   double *F, double *U, const double *box);


    private:
        void writePOSCAR(long N, const double *R, const int *atomicNrs,
                         const double *box);
        void readFU(long N, double *F, double *U);
        void spawnVASP();
        bool vaspRunning();
        static bool firstRun;
        static long vaspRunCount;
        static pid_t vaspPID;
};

#endif


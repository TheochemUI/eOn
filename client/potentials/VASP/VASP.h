
#ifndef VASP_POTENTIAL
#define VASP_POTENTIAL

#include <iostream>
#include <cstdio>
#include "stdlib.h"

#include "../../PotentialsInterface.h"

class VASP : public PotentialsInterface
{

    public:
        VASP(void);
        void initialize() {};
        void cleanMemory(void);    
        void force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box);


    private:
        void writeNEWCAR(long N, const double *R, long const *atomicNrs, const double *box);
        void readFU(long N, double *F, double *U);
        long vaspRunCount;

};

#endif


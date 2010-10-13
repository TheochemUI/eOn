#ifndef VASP_POTENTIAL
#define VASP_POTENTIAL

#include "../../PotentialsInterface.h"

class VASP : public PotentialsInterface
{

    public:
        VASP(void);
		~VASP();
        void initialize() {};
        void cleanMemory(void);    
        void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);


    private:
        void writeNEWCAR(long N, const double *R, int const *atomicNrs, const double *box);
        void readFU(long N, double *F, double *U);
        static long vaspRunCount;
		static bool vaspRunning;

};

#endif


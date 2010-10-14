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
        void force(long N, const double *R, const int *atomicNrs, 
                   double *F, double *U, const double *box);


    private:
        void writeNEWCAR(long N, const double *R, const int *atomicNrs,
                         const double *box);
        void readFU(long N, double *F, double *U);
        void spawnVASP();
        bool vaspRunning();
        static bool firstRun;
        static long vaspRunCount;
        static pid_t vaspPID;
};

#endif


#ifndef bopfox_POTENTIAL
#define bopfox_POTENTIAL

#include "../../PotentialsInterface.h"

class bopfox : public PotentialsInterface
{

    public:
        bopfox(void);
		~bopfox();
        void initialize() {};
        void cleanMemory(void);    
        void force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box);


    private:
        void writeFOX(long N, const double *R, long const *atomicNrs, const double *box);
        void readFU(long N, double *F, double *U);

};

#endif


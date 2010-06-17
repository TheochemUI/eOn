#ifndef NEWPOTENTIAL_INTERFACE
#define NEWPOTENTIAL_INTERFACE

#include "../../PotentialsInterface.h"

/** Template to use if user want to provide potential.*/
class NewPotential : public PotentialsInterface{

private:
//	Variables
    double fake1;
    double fake2;
    
public:
// Functions
	// constructor and destructor
    NewPotential(void);
	
    // To satify interface
    void initialize(void);    
    void cleanMemory(void);    

    void force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box);
};
#endif


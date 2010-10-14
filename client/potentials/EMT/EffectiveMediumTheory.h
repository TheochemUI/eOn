// serves as an interface between emt potentials provided by CamposASE and dynamics provided by EON

#ifndef EFFECTIVE_MEDIUM_THEORY
#define EFFECTIVE_MEDIUM_THEORY

#include"Asap/Atoms.h"
#include"Asap/EMT.h"
#include"Asap/SuperCell.h"
#include"Asap/Vec.h"

//#include "../../system_unit.h" // unit converters
#include "../../PotentialsInterface.h"

/** EMT potential. Inspect the EMT_parms.h to see what the EMT potential is hardcoded to describe.*/
class EffectiveMediumTheory : public PotentialsInterface {

private:
//	Variables
	long numberOfAtoms;
	bool periodicity[3];
	Atoms *AtomsObj;
	EMT *EMTObj;
	SuperCell *SuperCellObj;

public:
// Functions
	// constructor and destructor
    EffectiveMediumTheory(void);
    ~EffectiveMediumTheory(void) {};
    void cleanMemory(void);
    
    // To satify interface
    void initialize() {};
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
};
#endif

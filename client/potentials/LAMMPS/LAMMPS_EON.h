// serves as an interface between LAMMPS potentials maintained by SANDIA

#ifndef LAMMPS_EON
#define LAMMPS_EON

#include "string.h"

#include "filesFromLAMMPS/mpi.h"
#include "filesFromLAMMPS/lammps.h"         // these are LAMMPS include files
#include "filesFromLAMMPS/input.h"
#include "filesFromLAMMPS/atom.h"
#include "filesFromLAMMPS/library.h"

//#include "LAMMPSsrc/mpi.h"
//#include "LAMMPSsrc/lammps.h"         // these are LAMMPS include files
//#include "LAMMPSsrc/input.h"
//#include "LAMMPSsrc/atom.h"
//#include "LAMMPSsrc/library.h"

#include "../../Potential.h"

static LAMMPS_NS::LAMMPS* LAMMPSObj = NULL;

/** LAMMPS potentials.*/
class lammps_eon : public Potential {

private:
//	Variables
	long numberOfAtoms;
// Functions
	void makeNewLAMMPS(long N, const double *box);

public:
// Functions
	// constructor and destructor
    lammps_eon(void);
    ~lammps_eon(void) {};
    void cleanMemory(void);
    
    // To satify interface
    void initialize() {};
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
};
#endif

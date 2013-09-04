// serves as an interface between LAMMPS potentials maintained by SANDIA

#ifndef LAMMPS_EON
#define LAMMPS_EON

#include "../../Potential.h"
#include "../../Parameters.h"

class lammps_eon : public Potential {

private:
	long numberOfAtoms;
	double oldBox[9];
	void *LAMMPSObj;
	void makeNewLAMMPS(long N, const double *R,  const int *atomicNrs, const double *box);
	Parameters *parameters;

public:
    lammps_eon(Parameters *p);
    ~lammps_eon(void);
    void cleanMemory(void);
    
    void initialize() {};
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
};
#endif

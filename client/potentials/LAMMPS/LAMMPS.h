// serves as an interface between LAMMPS potentials maintained by SANDIA

#ifndef LAMMPS
#define LAMMPS

#include "../../Potential.h"
#include "../../Parameters.h"

class lammps : public Potential {

private:
	long numberOfAtoms;
	double oldBox[9];
	void *LAMMPSObj;
	void makeNewLAMMPS(long N, const double *R,  const int *atomicNrs, const double *box);
    Parameters *parameters;
    bool realunits;

public:
    lammps(Parameters *p);
    ~lammps(void);
    void cleanMemory(void);
    
    void initialize() {};
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
};
#endif

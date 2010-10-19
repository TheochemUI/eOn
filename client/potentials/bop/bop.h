#ifndef BOP_POTENTIAL
#define BOP_POTENTIAL

#include "../../PotentialsInterface.h"

    /** External function implemented in Fortran. Calculate interactions between molecules of water using forcefield bop.
    @param[in]	N           Number of atoms.
    @param[in]	R           Array to positions of the atoms in Angstrom.
    @param[out]	F           Array used to return the forces resulting from interactions between molecules. Forces are in eV/Angstrom.
    @param[out]	U           Pointer to energy in eV.
    @param[in]  bx, by, bz  Pointer to box dimensions in Angstrom.
    */
    
extern "C" 
{
    //void boplib_example_eam_(void);
    void bopini_(void);
    void boplib_calc_ef_(long *N, const double *R, double *box, double *U, double *F);    
    //void force_(const long int *N, const double *R, double *F, double *U, const double* bx, const double* by, const double* bz);
}    
//extern "C" 
//{
//    void potinit_();
//}    

/** bop potential.*/
class bop : public PotentialsInterface
{     
public:
    bop(void);
	
    // To satify interface
    void initialize(void);    
    void cleanMemory(void);    
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
    void writeFOX(long N, const double *R, int const *atomicNrs, const double *box);
    static bool initialized;
};
#endif


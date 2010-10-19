#include "bop.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

bool  bop::initialized = false;

const char *elements[] = {"Unknown", "H","He","Li","Be","B","C","N","O",
           "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",
           "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
           "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag",
           "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
           "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta",
           "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
           "Fr","Ra","Ac","Th","Pa","U", NULL};

bop::bop(void)
{
    return;
}

void bop::initialize(void)
{
    return;
}

void bop::cleanMemory(void)
{
    return;
}

void bop::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box)
{
    assert(N > 1);
    if(!initialized)
    {
        writeFOX(N, R, atomicNrs, box);
        bopini_();
        initialized = true;
    }

    double *atomEnergies = new double[N];
    // Initialize positions.

    double bopbox[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    bopbox[0] = box[0];
    bopbox[4] = box[1];
    bopbox[8] = box[2];

    // Call the FU function.
    boplib_calc_ef_(&N, R, bopbox, atomEnergies, F);      
    
    // bopfox returns energy in an array of per-atom energies.  Sum to get total energy.
    *U = 0.0;
    for(int i = 0; i < N; i++)
    {
        *U += atomEnergies[i];
    }
    
    // Cleanup
    delete atomEnergies;
    return;
}

void bop::writeFOX(long N, const double *R, int const *atomicNrs, const double *box)
{
    FILE *struc = fopen("struc.bx", "w");
    fprintf(struc, "StrucName = struc\n");
    fprintf(struc, "aLat = 1.0\n");
    fprintf(struc, "a1 =  %.8lf   %.8lf   %.8lf\n", box[0], 0.0, 0.0);
    fprintf(struc, "a2 =  %.8lf   %.8lf   %.8lf\n", 0.0, box[1], 0.0);
    fprintf(struc, "a3 =  %.8lf   %.8lf   %.8lf\n", 0.0, 0.0, box[2]);
    fprintf(struc, "coord = cartesian\n");
    for(int i = 0; i < N; i++)
    {
        fprintf(struc, "%s   %.8lf   %.8lf   %.8lf   %d\n", elements[atomicNrs[i]], R[i * 3 + 0], R[i * 3 + 1], R[i * 3 + 2], i);
    }
    fclose(struc);
    return;
}


//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "bop.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

bool  bop::initialized = false;
bool  bop::firstforce = false;

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
    // Redirect stdout to llout.
    fpos_t pos;
    int fd;
    fflush(stdout);
    fgetpos(stdout, &pos);
    fd = dup(fileno(stdout));
//    freopen("/dev/null", "w", stdout);
    freopen("./ll_out", "w", stdout);

    assert(N > 1);
    if(!initialized)
    {
        writeFOX(N, R, atomicNrs, box);
        bopini_();
        initialized = true;
    }

    double *atomEnergies = new double[N];
    // Initialize positions.

    //JR start: this does not compile, use old definition    
   // double *bopbox=box;
    double bopbox[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    bopbox[0] = box[0];
    bopbox[4] = box[4];
    bopbox[8] = box[8];
   //JR end
    
    //JR start : include a flag to quit eon if SC not converged in bopfox
    int scfconv;
    scfconv = 0;

    // Call the FU function.
    //boplib_calc_ef_(&N, R, bopbox, atomEnergies, F);      
    boplib_calc_ef_(&N, R, bopbox, atomEnergies, F, &scfconv); 
    if(scfconv == -1){
	    //print into bopfox output
	    printf("\n\n  !!ERROR!! SCF not converged in bopfox, exiting program..\n\n");
	    fflush(stdout);
	    //redirect back to standart out and print again
	    dup2(fd, fileno(stdout));
	    close(fd);
	    printf("\n\n  !!ERROR!! SCF not converged in bopfox, exiting program..\n\n");
	    exit(1);
    }
    //JR end    
    
   
    //JR: print out current structure and energies and forces
    //writeFOX(N, R, atomicNrs, box);
    //tsse_print_();
    //if(!firstforce)
    //{
      //  system("mv struc.bx struc-unrelaxed.bx");
      //  system("mv struc.EnFo.bx struc-unrelaxed.EnFo.bx");
      //  firstforce = true;
    //}
    //JR end
    
    // Redirect stdout to... stdout?? Something like that.
    fflush(stdout);
    dup2(fd, fileno(stdout));
    close(fd);
    clearerr(stdout);
    fsetpos(stdout, &pos);

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
    //JR start: modify printing of structure file in case of non-orthogonal box!
    fprintf(struc, "a1 =  %.8lf   %.8lf   %.8lf\n", box[0], box[1], box[2]);
    fprintf(struc, "a2 =  %.8lf   %.8lf   %.8lf\n", box[3], box[4], box[5]);
    fprintf(struc, "a3 =  %.8lf   %.8lf   %.8lf\n", box[6], box[7], box[8]);
    //JR end
    fprintf(struc, "coord = cartesian\n");
    for(int i = 0; i < N; i++)
    {
        fprintf(struc, "%s   %.8lf   %.8lf   %.8lf   %d\n", elements[atomicNrs[i]], R[i * 3 + 0], R[i * 3 + 1], R[i * 3 + 2], i);
    }
    //JR start: quick fix for magnetic Fe
    fprintf(struc, "magnetisation = defined\n");
    for(int i=0; i < N; i++){
        fprintf(struc, "0.0 0.0 2.5\n");
    }
    fprintf(struc, "\n");
    //JR end
    fclose(struc);
    return;
}


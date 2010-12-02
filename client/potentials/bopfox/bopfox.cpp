//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <unistd.h>
#ifdef WIN32
#include <windows.h>
#define sleep(n) Sleep(1000 * n)
#endif

#include "bopfox.h"

const char *elementNames[] = {"Unknown", "H","He","Li","Be","B","C","N","O",
           "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",
           "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
           "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag",
           "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
           "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta",
           "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
           "Fr","Ra","Ac","Th","Pa","U", NULL};

bopfox::bopfox(void)
{
    return;
}

void bopfox::cleanMemory(void)
{
    return;
}

bopfox::~bopfox()
{
	cleanMemory();
}


void bopfox::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box)
{
    writeFOX(N, R, atomicNrs, box);
	printf("bopfox force call\n");
    system("bopfox >> llout 2>&1");
    readFU(N, F, U);
    return;
}


void bopfox::writeFOX(long N, const double *R, int const *atomicNrs, const double *box)
{
    FILE *struc = fopen("struc.bx", "w");
    fprintf(struc, "StrucName = struc\n");
    fprintf(struc, "aLat = 1.0\n");
    fprintf(struc, "a1 =  %.8lf   %.8lf   %.8lf\n", box[0], box[1], box[2]);
    fprintf(struc, "a2 =  %.8lf   %.8lf   %.8lf\n", box[3], box[4], box[5]);
    fprintf(struc, "a3 =  %.8lf   %.8lf   %.8lf\n", box[6], box[7], box[8]);
    fprintf(struc, "coord = cartesian\n");
    for(int i = 0; i < N; i++)
    {
        fprintf(struc, "%s   %.8lf   %.8lf   %.8lf   %d\n", elementNames[atomicNrs[i]], R[i * 3 + 0], R[i * 3 + 1], R[i * 3 + 2], i);
    }
    fclose(struc);
    return;
}


void bopfox::readFU(long N, double *F, double *U)
{
    FILE *FU;
    FU = fopen("struc.EnFo.bx", "r");
    fscanf(FU, "%lf", U);
    for(int i = 0; i < N; i++)
    {
        fscanf(FU, "%lf %lf %lf", &F[i * 3 + 0], &F[i * 3 + 1], &F[i * 3 + 2]);
    }
    fclose(FU);
    return;
}






















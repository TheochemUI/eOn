#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <unistd.h>
#ifdef WIN32
#include <windows.h>
#define sleep(n) Sleep(1000 * n)
#endif

#include "VASP.h"

VASP::VASP(void)
{
    printf("INITED\n");
	vaspRunCount = 0;
    return;
}

void VASP::cleanMemory(void)
{
    printf("CLEANED\n");
	system("echo LABORT = .TRUE. > STOPCAR");
    return;
}


void VASP::force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box)
{
    writeNEWCAR(N, R, atomicNrs, box);
    if(vaspRunCount == 0)
    {
        system("vasp >> vaspOutput");
    }
    while(access("FU", F_OK) == -1)
    {
        sleep(1);
    }
    readFU(N, F, U);
    system("rm FU");
    vaspRunCount++;
    return;
}


void VASP::writeNEWCAR(long N, const double *R, long const *atomicNrs, const double *box)
{
    // Positions are scaled 
    long i = 0;
    long i_old = 0;
    FILE *NEWCAR;
    
    if(vaspRunCount == 0)
    {
        NEWCAR = fopen("POSCAR","w");
    }
    else
    {
        NEWCAR = fopen("NEWCAR","w");
    }

    // header line (treated as a comment)
    i_old = 0;
    fprintf(NEWCAR, "%li ", atomicNrs[0]);
    for(i = 0; i < N; i++)
    {
        if(atomicNrs[i] != atomicNrs[i_old])
        {
            fprintf(NEWCAR, "%li ", atomicNrs[i]);
            i_old = i;
        }
    }
    fprintf(NEWCAR, ": Atomic numbers\n");
    
    // boundary box
    fprintf(NEWCAR, "1.0\n");
    fprintf(NEWCAR, " %.8lf\t%.8lf\t%.8lf\n", box[0], 0.0, 0.0);
    fprintf(NEWCAR, " %.8lf\t%.8lf\t%.8lf\n", 0.0, box[1], 0.0);
    fprintf(NEWCAR, " %.8lf\t%.8lf\t%.8lf\n", 0.0, 0.0, box[2]);

    // the number of atoms of each of the the different atomic types
    i_old = 0;
    for(i = 0; i < N; i++)
    {
        if(atomicNrs[i] != atomicNrs[i_old])
        {
            fprintf(NEWCAR, "%li ", i - i_old);
            i_old = i;
        }
    }
    fprintf(NEWCAR, "%li\n", N - i_old);

    // coordinates for all atoms
    fprintf(NEWCAR, "Cartesian\n");
    for(i = 0; i < N; i++)
    {
        fprintf(NEWCAR, "%.19lf\t%.19lf\t%.19lf\t F F F\n", R[i * 3 + 0], R[i * 3 + 1],  R[i * 3 + 2]);
    }
    fclose(NEWCAR);
    return;
}


void VASP::readFU(long N, double *F, double *U)
{
    FILE *FU;
    FU = fopen("FU", "r");
    
    fscanf(FU, "%lf", U);
    
    for(int i = 0; i < N; i++)
    {
        fscanf(FU, "%lf %lf %lf", &F[i * 3 + 0], &F[i * 3 + 1], &F[i * 3 + 2]);
    }
   fclose(FU);
    return;
}






















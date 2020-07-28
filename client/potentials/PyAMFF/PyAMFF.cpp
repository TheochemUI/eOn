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
#include <errno.h>
#include <Python.h>
#include <string>

#ifdef WIN32
#include <windows.h>
#define sleep(n) Sleep(1000 * n)
//#define popen _popen
#else
#include <sys/wait.h>
#include <fcntl.h>
#endif

#include "PyAMFF.h"

PyAMFF::PyAMFF(void)
{
    // deleting leftovers from previous run
    system("rm -f POSCAR");  
    system("rm -f train.traj");  
    Py_Initialize();
    return;
}

void PyAMFF::cleanMemory(void)
{
    return;
}


std::string PyAMFF::num2sym(int Z){
    const char *elementArray[] = {"Unknown", "H","He","Li","Be","B","C","N","O",
           "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",
           "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
           "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag",
           "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
           "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta",
           "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
           "Fr","Ra","Ac","Th","Pa","U", NULL};
    std::string name = elementArray[Z];
    return name;
}

PyAMFF::~PyAMFF()
{
    Py_Finalize();
    cleanMemory();
}


void PyAMFF::force(long N, const double *R, const int *atomicNrs, double *F, 
                 double *U, const double *box)
{
    writePOSCAR(N, R, atomicNrs, box);

    FILE *fd = fopen("calc.py", "r");
    PyRun_SimpleFileEx(fd, "calc.py",1);

    readFU(N, F, U);
    remove("FU");
    remove("POSCAR");
    remove("train.traj");
    return;
}


void PyAMFF::writePOSCAR(long N, const double *R, const int *atomicNrs, 
                       const double *box)
{
    // Positions are scaled 
    long i = 0;
    long i_old = 0;
    //string sym;
    FILE *POSCAR;
    
    POSCAR = fopen("POSCAR","w");

    // header line (treated as a comment)
    i_old = 0;
    string sym = num2sym(atomicNrs[i]);

    fprintf(POSCAR, "%s ", sym.c_str());

    for(i = 0; i < N; i++)
    {
        if(atomicNrs[i] != atomicNrs[i_old])
        {
            string sym = num2sym(atomicNrs[i]);
            fprintf(POSCAR, "%s ",sym.c_str());
           
            i_old = i;
        }
    }
    fprintf(POSCAR, "\n");
    
    // boundary box
    fprintf(POSCAR, "1.0\n");
    fprintf(POSCAR, " %.8f\t%.8f\t%.8f\n", box[0], box[1], box[2]);
    fprintf(POSCAR, " %.8f\t%.8f\t%.8f\n", box[3], box[4], box[5]);
    fprintf(POSCAR, " %.8f\t%.8f\t%.8f\n", box[6], box[7], box[8]);

    // the number of atoms of each different atomic type
    i_old = 0;
    for(i = 0; i < N; i++)
    {
        if(atomicNrs[i] != atomicNrs[i_old])
        {
            fprintf(POSCAR, "%li ", i - i_old);
            i_old = i;
        }
    }
    fprintf(POSCAR, "%li\n", N - i_old);

    // coordinates for all atoms
    fprintf(POSCAR, "Cartesian\n");
    for(i = 0; i < N; i++)
    {
        fprintf(POSCAR, "%.19f\t%.19f\t%.19f\t T T T\n", R[i * 3 + 0], R[i * 3 + 1],  R[i * 3 + 2]);
    }
    fclose(POSCAR);

    FILE *fd = fopen("convert.py", "r");
    PyRun_SimpleFileEx(fd, "convert.py",1);

    return;
}


void PyAMFF::readFU(long N, double *F, double *U)
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

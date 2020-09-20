//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------
#include "PyAMFF.h"
#include <set>
#include <iostream>
#include <vector>

extern "C" void calc_eon(long *nAtoms, const double [], const double [], const int [], double [], double *U, int *num_elements, int []);


PyAMFF::PyAMFF(void)
{
    return;
}


void PyAMFF::cleanMemory(void)
{
    return;
}


void PyAMFF::force(long N, const double *R, const int *atomicNrs, double *F, 
                 double *U, const double *box)
{
//    int i;
//    const char *atomicSymbols[N];
//    int numUnique;

//    for (i=0; i < N; i++)
//    {
//        atomicSymbols[i] = atomicNumber2symbol(atomicNrs[i]);
//        cout << atomicNumber2symbol(atomicNrs[i]) << endl;
//        cout << "symbols" << endl/;
//        cout << symbols[i] << endl;
//        cout << i << endl;
//    }
    vector<int> unique_atomicNrs;
    for (int i=0; i<N; i++)
    {
        int j;
        for (j=0; j<i; j++)
            if (atomicNrs[i] == atomicNrs[j])
                break;
        if (i == j)
            unique_atomicNrs.push_back(atomicNrs[i]);
    }

    int num_elements;
    num_elements = unique_atomicNrs.size();
    int unique[num_elements];
    copy(unique_atomicNrs.begin(),unique_atomicNrs.end(),unique);
    
    calc_eon(&N, R, box, atomicNrs, F, U, &num_elements, unique);
    

    return;
}




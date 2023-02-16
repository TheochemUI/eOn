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

extern "C" void read_mlffParas(long *nAtoms, int *num_elements, int *max_fps, const int [], int []);

extern "C" void cleanup();

extern "C" void nncleanup();

extern "C" void prepfNN( long *nAtoms, int *num_elements, int *max_fps, const int [], int []);

PyAMFF::PyAMFF(void)
{
    new_pyamff = true; 
    return;
}

PyAMFF::~PyAMFF(){
    cleanMemory();
}

void PyAMFF::cleanMemory(void)
{
    if(new_pyamff != true ){
        new_pyamff = true;
        //also call NN cleanup
        cleanup();
        cout << "fp cleaned up" << endl;
        nncleanup();
        cout << "nn cleaned up" << endl;
    }
    return;
}

std::pair<double, AtomMatrix> ExtPot::get_ef(const AtomMatrix pos,
                                             const VectorXi atmnrs,
                                             const Matrix3d m_box) {
    double energy{0};
    long N{pos.rows()};
    AtomMatrix forces{Eigen::MatrixXd::Zero(N, 3)};
    const double *R = pos.data();
    const double *box = m_box.data();
    const int *atomicNrs = atmnrs.data();
    double *F = forces.data();
    double *U = &energy;

    int max_fps;
    max_fps = 100;
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
    //cout << "num_unique_elements" << num_elements << endl;
    int unique[num_elements];
    copy(unique_atomicNrs.begin(),unique_atomicNrs.end(),unique);
    

    if (new_pyamff == true){
        cout << "reading mlff in c" << endl;
        read_mlffParas(&N, &num_elements, &max_fps, atomicNrs, unique);
        cout << "prepping fnn in c" << endl;
        prepfNN(&N, &num_elements, &max_fps, atomicNrs, unique);
        cout << "fnn prepped!" << endl;
    }
    new_pyamff = false;
    calc_eon(&N, R, box, atomicNrs, F, U, &num_elements, unique);
    

    return std::make_pair(energy, forces);
}

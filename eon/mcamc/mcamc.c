#include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "qd/dd_real.h"

using namespace Eigen;

typedef Matrix<dd_real,Dynamic,Dynamic> MatrixXdd;
typedef Matrix<dd_real,Dynamic,1> VectorXdd;

extern "C" {
    void solve(int, double *, int, double *, double *, double *, double *);
}

void solve(int Qsize, double *Qflat, int Rcols, double *Rflat, double *c_in, double *B, double *t) {

    int i, row, col;
    unsigned int oldcw;
    fpu_fix_start(&oldcw);

    // Build an Eigen vector out of c_in.
    VectorXdd c = VectorXdd(Qsize);
    for (i = 0; i < Qsize; i++) {
        c(i) = 1.0 / c_in[i];
    };

    // Convert the Qflat array into an Eigen matrix.
    MatrixXdd Q = MatrixXdd(Qsize, Qsize);
    for (row = 0; row < Qsize; row++) {
        for (col = 0; col < Qsize; col++) {
            Q(row, col) =  Qflat[row*Qsize + col];    
        }
    }

    // Convert the Rflat array into an Eigen matrix.
    MatrixXdd R = MatrixXdd(Qsize, Rcols);
    for (row = 0; row < Qsize; row++) {
        for (col = 0; col < Rcols; col++) {
            R(row,col) =  Rflat[row*Rcols + col];    
        }
    }

    MatrixXdd A = MatrixXdd(Qsize, Qsize);
    A.setIdentity();
    A = A - Q;

    PartialPivLU<MatrixXdd> lu = A.partialPivLu();

    VectorXdd t_calc = lu.solve(c);
    MatrixXdd B_calc = lu.solve(R);

    // Store the solution t_calc into the array t.
    for (i = 0; i < Qsize; i++) {
        t[i] = t_calc[i].x[0];
    }

    // Store the solution B_calc into the array B.
    for (row = 0; row < Qsize; row++) {
        for (col = 0; col < Rcols; col++) {
            B[row*Rcols+col] = B_calc(row, col).x[0];
        }
    }

    fpu_fix_end(&oldcw);

}


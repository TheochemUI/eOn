#include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <qd/dd_real.h>
#include <qd/qd_real.h>

using namespace Eigen;


extern "C" {
    void solve_float(int, double *, int, double *, double *, double *, double *);
    void solve_double(int, double *, int, double *, double *, double *, double *);
    void solve_double_double(int, double *, int, double *, double *, double *, double *);
    void solve_quad_double(int, double *, int, double *, double *, double *, double *);
}

template<class T> void solve_general(int Qsize, double *Qflat, int Rcols, double *Rflat, double *c_in, double *B, double *t) {
    typedef Matrix<T,Dynamic,Dynamic> MatrixXdd;
    typedef Matrix<T,Dynamic,1> VectorXdd;

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
    MatrixXdd R = MatrixXdd(Qsize, Rcols);

    for (row = 0; row < Qsize; row++) {
        T total = 0;
        for (col = 0; col < Qsize; col++) {
            Q(row, col) =  Qflat[row*Qsize + col];    
            total += Q(row, col);
        }
        for (col = 0; col < Rcols; col++) {
            R(row,col) =  Rflat[row*Rcols + col];    
            total += R(row, col);
        }
        Q.row(row) /= total;
        R.row(row) /= total;
    }

    MatrixXdd A = MatrixXdd(Qsize, Qsize);
    A.setIdentity();
    A = A - Q;

    PartialPivLU<MatrixXdd> lu = A.partialPivLu();

    VectorXdd t_calc = lu.solve(c);
    MatrixXdd B_calc = lu.solve(R);

    // Store the solution t_calc into the array t.
    for (i = 0; i < Qsize; i++) {
        t[i] = t_calc[i];
    }

    // Store the solution B_calc into the array B.
    for (row = 0; row < Qsize; row++) {
        for (col = 0; col < Rcols; col++) {
            B[row*Rcols+col] = B_calc(row, col);
        }
    }

    fpu_fix_end(&oldcw);

}

void solve_float(int Qsize, double *Qflat, int Rcols, double *Rflat, 
                  double *c_in, double *B, double *t) {
    solve_general<float>(Qsize, Qflat, Rcols, Rflat, c_in, B, t);
}

void solve_double(int Qsize, double *Qflat, int Rcols, double *Rflat, 
                  double *c_in, double *B, double *t) {
    solve_general<double>(Qsize, Qflat, Rcols, Rflat, c_in, B, t);
}

void solve_double_double(int Qsize, double *Qflat, int Rcols, double *Rflat, 
                  double *c_in, double *B, double *t) {
    solve_general<dd_real>(Qsize, Qflat, Rcols, Rflat, c_in, B, t);
}

void solve_quad_double(int Qsize, double *Qflat, int Rcols, double *Rflat, 
                  double *c_in, double *B, double *t) {
    solve_general<qd_real>(Qsize, Qflat, Rcols, Rflat, c_in, B, t);
}


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>

#include "QSC.h"
#include "Parameters.h"

/* 
 * This is an implementation of the Quantum Sutton-Chen Potential,
 * which is an EAM type potential with the following functional form:
 * EAM functional: F_i(rho_i) = c*sqrt(rho_i)
 * EAM density:    rho_i(r_ij) = sum_(i!=j) (a/r_ij)^m
 * Pair potential: V(r_ij) = (a/r_ij)^n
 */


QSC::QSC()
{
    return;
}

void QSC::initialize()
{
}

void QSC::cleanMemory()
{
}


void QSC::force(long N, const double *R, const long *atomicNrs, double *F,
                double *U, const double *box)
{
    *U = 0.0;
    double *rho = new double[N];

    /* Potential and Density Calculation */
    for (int i=0; i<N; i++) {
        double repulsive=0.0;
        double attractive=0.0;
        rho[i]=0.0;
        for (int j=0; j<N; j++) {
            if (i==j) continue;
            double r_ij = distance(box, R, i, j);
            qsc_parameters p;
            p = get_qsc_parameters(atomicNrs[i], atomicNrs[j]); 
            if (r_ij > 2*p.a) continue;
            repulsive += p.epsilon*pair_potential(r_ij, p.a, p.n); 
            rho[i] += pair_potential(r_ij, p.a, p.m);
        }
        repulsive *= 0.5;
        qsc_parameters p=get_qsc_parameters(atomicNrs[i], atomicNrs[i]); 
        attractive = p.c * p.epsilon * sqrt(rho[i]);
        *U += repulsive - attractive;
        //printf("c=%f,e=%f,a=%f,m=%f,n=%f\n",p.c,p.epsilon,p.a,p.m,p.n);
        //printf("attractive=%f repulsive=%f\n", attractive, repulsive);
    }

    //printf("%10.4f %10.4f\n", distance(box, R, 0, 1), *U);

    /* Forces */
    for(int i=0;i<N;i++){
        F[ 3*i ] = 0;
        F[3*i+1] = 0;
        F[3*i+2] = 0;
    }
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if (i==j) continue;
            qsc_parameters p;
            p = get_qsc_parameters(atomicNrs[i], atomicNrs[j]); 
            double r_ij = distance(box, R, i, j);
            if (r_ij > 2*p.a) continue;
            double mag_force = -p.epsilon/r_ij*
                               (p.n*pair_potential(r_ij, p.a, p.n)-
                                - p.c*p.m/2.0*(sqrt(rho[i])*sqrt(rho[j]))*
                                pair_potential(r_ij, p.a, p.n));
            F[3*i]   += mag_force * (R[3*i]   - R[3*j]  )/r_ij;
            F[3*i+1] += mag_force * (R[3*i+1] - R[3*j+1])/r_ij;
            F[3*i+2] += mag_force * (R[3*i+2] - R[3*j+2])/r_ij;
        }
        //printf("F[%2i] = %10.4g\n", 3*i, F[3*i]);
        //printf("F[%2i] = %10.4g\n", 3*i+1,F[3*i+1]);
        //printf("F[%2i] = %10.4g\n", 3*i+2,F[3*i+2]);
    }
    //printf("F1x = %10.4f F2x = %10.4f\n", F[0], F[3]);

    delete rho;
}

double QSC::pair_potential(double r, double a, double n)
{
    return pow(a/r, n);
}

double QSC::distance(const double *box, const double *R, int i, int j)
{
    double diffR, diffRX, diffRY, diffRZ;

    diffRX = R[3*i]   - R[3*j];
    diffRY = R[3*i+1] - R[3*j+1];
    diffRZ = R[3*i+2] - R[3*j+2];

    diffRX = diffRX-box[0]*floor(diffRX/box[0]+0.5); 
    diffRY = diffRY-box[1]*floor(diffRY/box[1]+0.5);
    diffRZ = diffRZ-box[2]*floor(diffRZ/box[2]+0.5);
    
    diffR = sqrt(diffRX*diffRX+diffRY*diffRY+diffRZ*diffRZ);

    return diffR;
}


QSC::qsc_parameters QSC::get_qsc_parameters(int a, int b)
{
    qsc_parameters p;
    int ia, ib;
    bool founda=false, foundb=false;

    for (int i=0;i<NPARAMS;i++) {
        if (qsc_element_params[i].Z == a) {
            founda = true;
            ia = i;
        }
        if (qsc_element_params[i].Z == b) {
            foundb = true;
            ib = i;
        }
    }

    if (founda == false and foundb == false) {
        /* This sucks. We need to have a way to alert user that the 
         * parameters are missing */
        throw 14324;
    }

    p.epsilon = sqrt(qsc_element_params[ia].epsilon * 
                     qsc_element_params[ib].epsilon);
    p.a = (qsc_element_params[ia].a + qsc_element_params[ib].a)/2.0;
    p.m = (qsc_element_params[ia].m + qsc_element_params[ib].m)/2.0;
    p.n = (qsc_element_params[ia].n + qsc_element_params[ib].n)/2.0;
    p.c = qsc_element_params[ia].c;

    return p;
}

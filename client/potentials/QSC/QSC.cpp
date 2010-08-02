#include <math.h>
#include <stdio.h>
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

    /* Zero out Forces and densities */
    for(int i=0;i<N;i++){
        F[3*i  ] = 0;
        F[3*i+1] = 0;
        F[3*i+2] = 0;
        rho[i] = 0.0;
    }

    /* Calculate the local density (rho[i]) for each atom 
     * and the potential energy (U). */
    for (int i=0; i<N; i++) {
        double pair_term=0.0;
        qsc_parameters p_ij, p_ii;
        for (int j=i+1; j<N; j++) {
            double r_ij = distance(box, R, i, j);
            /* Get the mixed parameters */
            p_ij = get_qsc_parameters(atomicNrs[i], atomicNrs[j]); 
            if (r_ij > 2*p_ij.a) continue;
            double delta_rho = pair_potential(r_ij, p_ij.a, p_ij.m);
            rho[i] += delta_rho;
            rho[j] += delta_rho;
            pair_term += p_ij.epsilon*pair_potential(r_ij, p_ij.a, p_ij.n);
        }
        /* Get the parameters for element i */
        p_ii = get_qsc_parameters(atomicNrs[i], atomicNrs[i]); 
        double embedding_term = p_ii.c * p_ii.epsilon * sqrt(rho[i]);
        *U += pair_term - embedding_term;
    }

    /* Forces Calculation */
    for (int i=0; i<N; i++) {
        for (int j=i+1; j<N; j++) {
            qsc_parameters p;
            p = get_qsc_parameters(atomicNrs[i], atomicNrs[j]); 

            double r_ij = distance(box, R, i, j);
            if (r_ij > 2*p.a) continue;

            double Fij;
            Fij =       p.epsilon *
                        (p.n*pair_potential(r_ij, p.a, p.n) -
                         (p.c*p.m*0.5*(pow(rho[i],-0.5)+pow(rho[j],-0.5)) *
                         pair_potential(r_ij, p.a, p.m)))/(r_ij*r_ij);

            double diffx,diffy,diffz, Fijx, Fijy, Fijz;
            diffx = R[3*i  ] - R[3*j  ];
            diffy = R[3*i+1] - R[3*j+1];
            diffz = R[3*i+2] - R[3*j+2];
            /* Orthogonal PBC */
            diffx = diffx-box[0]*floor(diffx/box[0]+0.5); 
            diffy = diffy-box[1]*floor(diffy/box[1]+0.5);
            diffz = diffz-box[2]*floor(diffz/box[2]+0.5);

            Fijx = Fij * diffx;
            Fijy = Fij * diffy;
            Fijz = Fij * diffz;

            F[3*i]   += Fijx;
            F[3*i+1] += Fijy;
            F[3*i+2] += Fijz;
            F[3*j]   -= Fijx;
            F[3*j+1] -= Fijy;
            F[3*j+2] -= Fijz;
        }
    }

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

    /* Orthogonal PBC */
    diffRX = diffRX-box[0]*floor(diffRX/box[0]+0.5); 
    diffRY = diffRY-box[1]*floor(diffRY/box[1]+0.5);
    diffRZ = diffRZ-box[2]*floor(diffRZ/box[2]+0.5);
    
    diffR = sqrt(diffRX*diffRX+diffRY*diffRY+diffRZ*diffRZ);

    return diffR;
}


QSC::qsc_parameters QSC::get_qsc_parameters(int element_a, int element_b)
{
    qsc_parameters p;
    int i=0, ia=-1, ib=-1, Z;

    while (true) {
        Z = qsc_element_params[i].Z;

        /* -1 is the element number of the parameter structure at
         * the end of the array. */
        if (Z == -1) {
            /* Need better way to tell user there are missing
             * parameters. */
            throw 14324;
        }

        if (Z == element_a) {
            ia = i;
        } 

        if (Z == element_b) {
            ib = i;
        }

        /* If we have found both parameters we are done. */
        if (ia != -1 and ib != -1) {
            break;
        }
        i++;
    }

    /* Mixing rules */
    p.epsilon = sqrt(qsc_element_params[ia].epsilon * 
                     qsc_element_params[ib].epsilon);
    p.a = (qsc_element_params[ia].a + qsc_element_params[ib].a)/2.0;
    p.m = (qsc_element_params[ia].m + qsc_element_params[ib].m)/2.0;
    p.n = (qsc_element_params[ia].n + qsc_element_params[ib].n)/2.0;
    p.c = qsc_element_params[ia].c;

    return p;
}

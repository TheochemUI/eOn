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
    cutoff = 6.0;
    verlet_skin = 0.01;
    init=false;
    return;
}

void QSC::initialize()
{
}

void QSC::cleanMemory()
{
    if (init==true) {
        delete oldR;
    }
}

void QSC::new_vlist(long N, const double *R, const double *box)
{
    vlist = new int*[N];
    nlist = new int[N];
    distances = new struct distance*[N];
    double rv = cutoff+verlet_skin;

    for (int i=0; i<N; i++) {
        vlist[i] = new int[N];
        distances[i] = new struct distance[N];
        nlist[i] = 0;
        for (int j=i+1;j<N; j++) {
            calc_distance(box, R, i, j, &distances[i][j]);
            if (distances[i][j].r <= rv) {
                vlist[i][nlist[i]] = j;
                nlist[i] += 1;
            }
        }
    }
}

void QSC::calc_distance(const double *box, const double *R, int i, int j, struct distance *d)
{
        double diffRX = R[3*i]   - R[3*j];
        double diffRY = R[3*i+1] - R[3*j+1];
        double diffRZ = R[3*i+2] - R[3*j+2];

        /* Orthogonal PBC */
        diffRX = diffRX-box[0]*floor(diffRX/box[0]+0.5); 
        diffRY = diffRY-box[1]*floor(diffRY/box[1]+0.5);
        diffRZ = diffRZ-box[2]*floor(diffRZ/box[2]+0.5);
        
        d->r = sqrt(diffRX*diffRX+diffRY*diffRY+diffRZ*diffRZ);
        d->d[0] = diffRX;
        d->d[1] = diffRY;
        d->d[2] = diffRZ;
}

void QSC::update_vlist(long N, const double *R, const double *box) 
{
    bool update=false;
    for (int i=0;i<3*N;i++) {
        double diff = oldR[i]-R[i];
        diff = diff-box[0]*floor(diff/box[0]+0.5); 
        if (fabs(diff) > verlet_skin) {
            update=true;
            break;
        }
    }

    if (update==true) {
        new_vlist(N, R, box);
    }else{
        for (int i=0; i<N; i++) {
            for (int k=0; k<nlist[i]; k++) {
                int j = vlist[i][k];
                calc_distance(box, R, i, j, &distances[i][j]);
            }
        }
    }
}

void QSC::force(long N, const double *R, const long *atomicNrs, double *F,
                double *U, const double *box)
{
    if (init==false) {
        init = true;
        new_vlist(N, R, box);
        oldR = new double[3*N];
    }else{
        update_vlist(N, R, box);
    }

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
        qsc_parameters p_ij, p_ii, p_jj;
        /* Get the parameters for element i */
        p_ii = get_qsc_parameters(atomicNrs[i], atomicNrs[i]); 
        for (int k=0; k<nlist[i]; k++) {
            int j = vlist[i][k];
            /* Get the parameters */
            p_ij = get_qsc_parameters(atomicNrs[i], atomicNrs[j]);
            if (distances[i][j].r > cutoff) continue;
            p_jj = get_qsc_parameters(atomicNrs[j], atomicNrs[j]);

            /* Take care of density */
            double delta_rho;
            delta_rho = pair_potential(distances[i][j].r, p_jj.a, p_jj.m);
            rho[i] += delta_rho;

            delta_rho = pair_potential(distances[i][j].r, p_ii.a, p_ii.m);
            rho[j] += delta_rho;

            /* Repulsive pair term */
            pair_term += p_ij.epsilon*pair_potential(distances[i][j].r, p_ij.a, p_ij.n);
        }
        double embedding_term = p_ii.c * p_ii.epsilon * sqrt(rho[i]);
        *U += pair_term - embedding_term;
    }

    /* Forces Calculation */
    for (int i=0; i<N; i++) {
        for (int k=0; k<nlist[i]; k++) {
            int j = vlist[i][k];
            qsc_parameters p_ii, p_ij, p_jj;
            p_ii = get_qsc_parameters(atomicNrs[i], atomicNrs[i]); 
            p_ij = get_qsc_parameters(atomicNrs[i], atomicNrs[j]); 
            p_jj = get_qsc_parameters(atomicNrs[j], atomicNrs[j]); 

            double r_ij = distances[i][j].r;
            if (distances[i][j].r > cutoff) continue;

            double Fij;
            Fij  = p_ij.epsilon*p_ij.n*pair_potential(r_ij, p_ij.a, p_ij.n);
            Fij -= p_ii.epsilon*p_ii.c*p_jj.m*0.5*pow(rho[i],-0.5)*
                   pair_potential(r_ij, p_jj.a, p_jj.m);
            Fij -= p_jj.epsilon*p_jj.c*p_ii.m*0.5*pow(rho[j],-0.5)*
                   pair_potential(r_ij, p_ii.a, p_ii.m);
            Fij /= r_ij;


            double Fijx = Fij * distances[i][j].d[0]/r_ij;
            double Fijy = Fij * distances[i][j].d[1]/r_ij;
            double Fijz = Fij * distances[i][j].d[2]/r_ij;

            F[3*i]   += Fijx;
            F[3*i+1] += Fijy;
            F[3*i+2] += Fijz;
            F[3*j]   -= Fijx;
            F[3*j+1] -= Fijy;
            F[3*j+2] -= Fijz;
        }
    }

    delete rho;

    for (int i=0;i<3*N;i++) {
        oldR[i] = R[i];
    }
}

inline double QSC::pair_potential(double r, double a, double n)
{
    return pow(a/r, n);
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
    p.a = 0.5*(qsc_element_params[ia].a + qsc_element_params[ib].a);
    p.m = 0.5*(qsc_element_params[ia].m + qsc_element_params[ib].m);
    p.n = 0.5*(qsc_element_params[ia].n + qsc_element_params[ib].n);
    p.c = qsc_element_params[ia].c;

    return p;
}

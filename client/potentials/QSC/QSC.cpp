#include <math.h>
#include <stdio.h>
#include <stdlib.h>
int compare_ints(const void *a,const void *b);

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
    verlet_skin = 0.5;
    init=false;
    return;
}

QSC::~QSC()
{
    //cleanMemory();
}

void QSC::initialize() {}
void QSC::initialize(long N, const double *R, const int *atomicNrs,
                     const double *box) {
        /* Determine what elements are in the system to build a table
         * of the parameters. */
        int *sorted_atomic_numbers = new int[N];
        for (int i=0;i<N;i++) {
            sorted_atomic_numbers[i] = atomicNrs[i];
        }
        qsort(sorted_atomic_numbers, N, sizeof(int), compare_ints);
        int nunique=0;
        unique_elements[0] = sorted_atomic_numbers[0];
        for (int i=0;i<N;i++) {
            if (sorted_atomic_numbers[i] != unique_elements[nunique]) {
                nunique++;
                unique_elements[nunique] = sorted_atomic_numbers[i];
            }
        }
        nunique++;
        delete sorted_atomic_numbers;

        largest_element_num = unique_elements[nunique-1];

        qsc_param_cache = new qsc_parameters*[largest_element_num+1];
        /* Set all of them to NULL so we know which to delete later */
        for (int i=0; i<largest_element_num+1; i++) {
            qsc_param_cache[i] = NULL;
        }

        for (int i=0;i<nunique;i++) {
            int ui = unique_elements[i];
            qsc_param_cache[ui] = new qsc_parameters[largest_element_num+1];
            for (int j=0;j<nunique;j++) {
                int uj = unique_elements[j];
                qsc_param_cache[ui][uj] = get_qsc_parameters(ui, uj);
            }
        }

        /* Allocate memory */
        natomstoclear = N; //this is ugly but i need this in cleanMemory
        rho = new double[N];
        sqrtrho = new double[N];
        vlist = new int*[N];
        nlist = new int[N];
        distances = new struct distance*[N];
        oldR = new double[3*N];
        V = new double*[N];
        phi = new double*[N];
        for (int i=0;i<N;i++) {
            vlist[i] = new int[N];
            distances[i] = new struct distance[N];
            V[i] = new double[N];
            phi[i] = new double[N];
        }

        /* Make a new verlet list and populate distances[][] */
        new_vlist(N, R, box);

        init = true;
}


void QSC::cleanMemory()
{
    if (init==true) {
        delete rho;
        delete sqrtrho;
        delete oldR;
        delete nlist;
        for (int i=0; i<natomstoclear; i++) {
            delete vlist[i];
            delete distances[i];
            delete V[i];    
            delete phi[i];
        }
        for (int i=0; i<largest_element_num+1; i++) {
            if (qsc_param_cache[i] != NULL) {
                delete qsc_param_cache[i];
            }
        }
        delete qsc_param_cache;
        delete vlist;
        delete distances;
        delete V;    
        delete phi;

        init = false;
    }
}

void QSC::new_vlist(long N, const double *R, const double *box)
{
    double rv = cutoff+verlet_skin;

    for (int i=0; i<N; i++) {
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

    for (int i=0;i<3*N;i++) {
        oldR[i] = R[i];
    }
}

int compare_ints(const void *a,const void *b) 
{
    return ( *(int*)a - *(int*)b );
}

void QSC::energy(long N, const double *R, const int *atomicNrs, double *U,
                 const double *box) {
    if (init==false) {
        initialize(N, R, atomicNrs, box);
    }else{
        update_vlist(N, R, box);
    }

    *U = 0.0;
    for(int i=0;i<N;i++){
        rho[i] = 0.0;
    }

    /* Calculate the local density (rho[i]) for each atom 
     * and the potential energy (U). */
    for (int i=0; i<N; i++) {
        double pair_term=0.0;
        qsc_parameters p_ij, p_ii, p_jj;
        /* Get the parameters for element i */
        p_ii = qsc_param_cache[atomicNrs[i]][atomicNrs[i]];
        for (int k=0; k<nlist[i]; k++) {
            int j = vlist[i][k];
            if (distances[i][j].r > cutoff) continue;
            /* Get the parameters */
            p_ij = qsc_param_cache[atomicNrs[i]][atomicNrs[j]];
            p_jj = qsc_param_cache[atomicNrs[j]][atomicNrs[j]];

            /* Take care of density */
            phi[i][j] = pair_potential(distances[i][j].r, p_jj.a, p_jj.m);
            rho[i] += phi[i][j];
            phi[j][i] = pair_potential(distances[i][j].r, p_ii.a, p_ii.m);
            rho[j] += phi[j][i];

            /* Repulsive pair term */
            V[i][j] = p_ij.epsilon*
                      pair_potential(distances[i][j].r, p_ij.a, p_ij.n);

            pair_term += V[i][j];
        }
        sqrtrho[i] = sqrt(rho[i]);
        double embedding_term = p_ii.c * p_ii.epsilon * sqrtrho[i];
        *U += pair_term - embedding_term;
    }

}

void QSC::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, const double *box)
{
    if (init==false) {
        initialize(N, R, atomicNrs, box);
    }else{
        update_vlist(N, R, box);
    }

    energy(N, R, atomicNrs, U, box);


    /* Zero out Forces */
    for(int i=0;i<3*N;i++){
        F[i] = 0.0;
    }

    /* Forces Calculation */
    for (int i=0; i<N; i++) {
        for (int k=0; k<nlist[i]; k++) {
            int j = vlist[i][k];
            qsc_parameters p_ii, p_ij, p_jj;
            p_ii = qsc_param_cache[atomicNrs[i]][atomicNrs[i]];
            p_ij = qsc_param_cache[atomicNrs[i]][atomicNrs[j]];
            p_jj = qsc_param_cache[atomicNrs[j]][atomicNrs[j]];

            double r_ij = distances[i][j].r;
            if (distances[i][j].r > cutoff) continue;

            double Fij;
            Fij  = p_ij.n*V[i][j];
            Fij -= p_ii.epsilon*p_ii.c*p_jj.m*0.5*(1.0/sqrtrho[i])*phi[i][j];
            Fij -= p_jj.epsilon*p_jj.c*p_ii.m*0.5*(1.0/sqrtrho[j])*phi[j][i];
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
}

inline double QSC::pair_potential(double r, double a, double n)
{
    return pow(a/r, n);
}

void QSC::calc_distance(const double *box, const double *R, int i, int j, 
                        struct distance *d)
{
        double diffRX = R[3*i]   - R[3*j];
        double diffRY = R[3*i+1] - R[3*j+1];
        double diffRZ = R[3*i+2] - R[3*j+2];

        /* Orthogonal PBC */
        diffRX = diffRX-box[0]*floor(diffRX/box[0]+0.5); 
        diffRY = diffRY-box[4]*floor(diffRY/box[4]+0.5);
        diffRZ = diffRZ-box[8]*floor(diffRZ/box[8]+0.5);
        
        d->r = sqrt(diffRX*diffRX+diffRY*diffRY+diffRZ*diffRZ);
        d->d[0] = diffRX;
        d->d[1] = diffRY;
        d->d[2] = diffRZ;
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
    if (ia==ib) {
        return qsc_element_params[ia];
    }
    p.epsilon = sqrt(qsc_element_params[ia].epsilon * 
                     qsc_element_params[ib].epsilon);
    p.a = 0.5*(qsc_element_params[ia].a + qsc_element_params[ib].a);
    p.m = 0.5*(qsc_element_params[ia].m + qsc_element_params[ib].m);
    p.n = 0.5*(qsc_element_params[ia].n + qsc_element_params[ib].n);
    p.c = qsc_element_params[ia].c;
    p.Z = 0;

    return p;
}

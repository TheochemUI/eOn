#include "math.h"
#include "stdlib.h"
#include "stdio.h"

// Fn. declarations
void force(int, double *, double *, double *, double, double, double, int *);
void pbc_vdr(double*, double*, double*, double, double, double);
double pbc_mdr(double*, double*, double, double, double);

// The same routine as energy above, but also stores
// the forces.  It is cheaper to only call
// energy() if that is all that is required
void force(int n, double *x, double *f, double *u, double bx, double by, double bz, int *type) 
{

    double r_cut = 2.5;           // Cutoff distance
    *u = 0.0;

    double dr[3], dr6, dr2, mdr, rc2 = r_cut * r_cut, ff[3], force;
    double Etot, ELJ;
    int a, b, tpa, tpb;

    // Define half the box length //
    double C12[2][2], C6[2][2];   // 4.0 * eps * sig^12 and 4.0 * eps * sig^6 respectively

    // Set up the LJ parameters //
    double sigA, sigB, epsA, epsB, epsAB, sigAB, sig3, sig6;
    sigA = epsA = 1.0;
    epsAB = 1.5;
    sigAB = 0.8;
    epsB = 0.5;
    sigB = 0.88;
    sig3 = sigA * sigA * sigA;
    sig6 = sig3 * sig3;
    C6[0][0] = 4.0 * epsA * sig6;
    C12[0][0] = C6[0][0] * sig6;
    sig3 = sigB * sigB * sigB;
    sig6 = sig3 * sig3;
    C6[1][1] = 4.0 * epsB * sig6;
    C12[1][1] = C6[1][1] * sig6;
    sig3 = sigAB * sigAB * sigAB;
    sig6 = sig3 * sig3;
    C6[0][1] = 4.0 * epsAB * sig6;
    C12[0][1] = C6[0][1] * sig6;
    C6[1][0] = C6[0][1];
    C12[1][0] = C12[0][1];
    // Set up the shift forces and energies //
    double Esh[2][2], Fsh[2][2];  // Energy and force shifting values
    double dr3;
    dr3 = r_cut * r_cut * r_cut;
    dr6 = dr3 * dr3;
    Esh[0][0] = (C12[0][0] / dr6 - C6[0][0]) / dr6;
    Esh[1][1] = (C12[1][1] / dr6 - C6[1][1]) / dr6;
    Esh[0][1] = (C12[0][1] / dr6 - C6[0][1]) / dr6;
    Esh[1][0] = Esh[0][1];
    Fsh[0][0] = 6.0 * (2.0 * C12[0][0] / dr6 - C6[0][0]) / dr6 / r_cut;
    Fsh[1][1] = 6.0 * (2.0 * C12[1][1] / dr6 - C6[1][1]) / dr6 / r_cut;
    Fsh[0][1] = 6.0 * (2.0 * C12[0][1] / dr6 - C6[0][1]) / dr6 / r_cut;
    Fsh[1][0] = Fsh[0][1];

    int i;
    for(i = 0; i < n * 3; i++)
    {
        f[i] = 0.0;
    }

    for (a = 0; a < n - 1; a++) 
    {
        tpa = type[a];
        for (b = a + 1; b < n; b++) 
        {
            tpb = type[b];
            pbc_vdr(&x[a * 3], &x[b * 3], dr, bx, by, bz);
            dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if (dr2 > rc2)
                continue;
            mdr = sqrt(dr2);
            dr6 = dr2 * dr2 * dr2;
            ELJ = ((C12[tpa][tpb] / dr6 - C6[tpa][tpb]) / dr6) - Esh[tpa][tpb] - Fsh[tpa][tpb] * (r_cut - mdr);
            force = (6.0 * (2.0 * C12[tpa][tpb] / dr6 - C6[tpa][tpb]) / dr6 / dr2) - Fsh[tpa][tpb] / mdr;
            ff[0] = force * dr[0];
            ff[1] = force * dr[1];
            ff[2] = force * dr[2];
            f[a * 3 + 0] += ff[0];
            f[a * 3 + 1] += ff[1];
            f[a * 3 + 2] += ff[2];
            f[b * 3 + 0] -= ff[0];
            f[b * 3 + 1] -= ff[1];
            f[b * 3 + 2] -= ff[2];
            *u += ELJ;
        }
    }
}



// Returns the magnitude of the vector returned by pbc_vdr
double pbc_mdr(double *r1, double *r2, double bx, double by, double bz) 
{
    double dr[3], mdr;
    pbc_vdr(r1, r2, dr, bx, by, bz);
    mdr = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    return sqrt(mdr);
}



// Sets dr = r1 - r2 taking into account the PBCs
void pbc_vdr(double *r1, double *r2, double dr[3], double bx, double by, double bz) 
{
    double bxh = 0.5 * bx;
    double byh = 0.5 * by;
    double bzh = 0.5 * bz;
    dr[0] = r1[0] - r2[0];
    dr[1] = r1[1] - r2[1];
    dr[2] = r1[2] - r2[2];
    while (dr[0] >= bxh) dr[0] -= bx;
    while (dr[0] < -bxh) dr[0] += bx;
    while (dr[1] >= byh) dr[1] -= by;
    while (dr[1] < -byh) dr[1] += by;
    while (dr[2] >= bzh) dr[2] -= bz;
    while (dr[2] < -bzh) dr[2] += bz;
}














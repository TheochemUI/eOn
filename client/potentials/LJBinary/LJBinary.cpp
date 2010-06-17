#include "LJBinary.h"

LJBinary::LJBinary(){
    // Values from Andri
    this->setParameters(7830.0*8.61738573e-5, 2.2*2.47, 2.47);
}

LJBinary::LJBinary(double u0Recieved, double cuttOffRRecieved, double psiRecieved){
    this->setParameters(u0Recieved, cuttOffRRecieved, psiRecieved);
    return;
}

void LJBinary::cleanMemory(void){
    return;
}

// General Functions
void LJBinary::setParameters(double u0Recieved, double cuttOffRRecieved, double psiRecieved){
    u0 = u0Recieved;
    psi = psiRecieved;
    
    cuttOffR = cuttOffRRecieved;
    cuttOffU = 4*u0*(pow(psi/cuttOffR,12)-pow(psi/cuttOffR,6));
    return;
}

void LJBinary::force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box)
{
//    fprintf(stderr, "client_eon: entered LJBinary.force()\n");

//    fprintf(stderr, "client_eon: setting up LJ parameters\n");
    double r_cut = 2.5;           // Cutoff distance
    *U = 0.0;
	
    double bx = box[0];
    double by = box[1];
    double bz = box[2];

    double dr[3], dr6, dr2, mdr, rc2 = r_cut * r_cut, ff[3], force;
    double Etot, ELJ;
    long a, b;
    int tpa, tpb;

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

//    fprintf(stderr, "client_eon: zeroing the force vector\n");
    int i;
    for(i = 0; i < N * 3; i++)
    {
        F[i] = 0.0;
    }


//    fprintf(stderr, "client_eon: entering LJBinary a for loop\n");
    for (a = 0; a < N - 1; a++) 
    {
        tpa = atomicNrs[a] - 1;
//        fprintf(stderr, "client_eon: entering LJBinary b for loop\n");
        for (b = a + 1; b < N; b++) 
        {
            tpb = atomicNrs[b] - 1;
            pbc_vdr(&R[a * 3], &R[b * 3], dr, bx, by, bz);
            dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if (dr2 > rc2)
                continue;
            mdr = sqrt(dr2);
            dr6 = dr2 * dr2 * dr2;
//            fprintf(stderr, "client_eon: calculating ELJ @ (%d, %d)\n", tpa, tpb);
//            fprintf(stderr, "client_eon: C12[0][0]: %d\n", C12[0][0]);
//            fprintf(stderr, "client_eon: C12[tpa][tpb]: %d\n", C12[tpa][tpb]);
//            fprintf(stderr, "client_eon: dr6: %d\n", dr6);
//            fprintf(stderr, "client_eon: C12[tpa][tpb] / dr6: %d\n", C12[tpa][tpb] / dr6);
//            fprintf(stderr, "client_eon: C6[tpa][tpb]: %d\n", C6[tpa][tpb]);
//            fprintf(stderr, "client_eon: ((C12[tpa][tpb] / dr6 - C6[tpa][tpb]) / dr6): %d\n", ((C12[tpa][tpb] / dr6 - C6[tpa][tpb]) / dr6));
//            fprintf(stderr, "client_eon: -Fsh[tpa][tpb] * (r_cut - mdr): %d\n", -Fsh[tpa][tpb] * (r_cut - mdr));
            ELJ = ((C12[tpa][tpb] / dr6 - C6[tpa][tpb]) / dr6) - Esh[tpa][tpb] - Fsh[tpa][tpb] * (r_cut - mdr);
//            fprintf(stderr, "client_eon: Calculated ELJ\n");
            force = (6.0 * (2.0 * C12[tpa][tpb] / dr6 - C6[tpa][tpb]) / dr6 / dr2) - Fsh[tpa][tpb] / mdr;
            ff[0] = force * dr[0];
            ff[1] = force * dr[1];
            ff[2] = force * dr[2];
            F[a * 3 + 0] += ff[0];
            F[a * 3 + 1] += ff[1];
            F[a * 3 + 2] += ff[2];
            F[b * 3 + 0] -= ff[0];
            F[b * 3 + 1] -= ff[1];
            F[b * 3 + 2] -= ff[2];
            *U += ELJ;
        }
//        fprintf(stderr, "client_eon: exiting LJBinary b for loop\n");
    }
//    fprintf(stderr, "client_eon: exiting LJBinary a for loop\n");

//    fprintf(stderr, "client_eon: exiting LJBinary.force()\n");

}

// Returns the magnitude of the vector returned by pbc_vdr
double LJBinary::pbc_mdr(const double *r1, const double *r2, double bx, double by, double bz) 
{
    double dr[3], mdr;
    pbc_vdr(r1, r2, dr, bx, by, bz);
    mdr = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    return sqrt(mdr);
}

// Sets dr = r1 - r2 taking into account the PBCs
void LJBinary::pbc_vdr(const double *r1, const double *r2, double dr[3], double bx, double by, double bz) 
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


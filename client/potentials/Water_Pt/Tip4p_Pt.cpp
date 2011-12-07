#include "Tip4p_Pt.hpp"


Tip4p_Pt::Tip4p_Pt() :
Potential(),
forcefields::ZhuPhilpott<>(8.5, 1.0)
{}

void Tip4p_Pt::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box)
{
    double diagbox[3];
    diagbox[0]=box[0];
    diagbox[1]=box[4];
    diagbox[2]=box[8];
    int i=0;
    while (atomicNrs[i] == 1)
        i+=2;
    computeHH_O_Pt_(i/2, N - i*3/2, R, F, *U, diagbox, 0);
}

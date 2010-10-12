#include "zhu_philpott_for_eon.hpp"

ZpIce::ZpIce() :
    PotentialsInterface(),
    forcefields::ZhuPhilpott<>(8.5, 1.0)
{}

void ZpIce::force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box)
{
    int i=0;
    while (atomicNrs[i] == 1)
        i+=2;
    computeHH_O_Pt_(i/2, N - i*3/2, R, F, *U, box, 0);
}

Tip4p::Tip4p() :
    PotentialsInterface(),
    forcefields::Tip4p(8.5, 1.0)
{}

void Tip4p::force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box)
{
    computeHH_O_(N, R, F, *U, box, 0);
}

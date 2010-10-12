/** @file
Wrapper for Fortran
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
*/
#include <fstream>
#include "zhu_philpott.hpp"
#include "zhu_philpott_parameters.hpp"

namespace {
    forcefields::ZhuPhilpott<forcefields::zhu_philpott_parameters::Iceland> potential(8.5, 1.0);
}

extern "C"
void zhu_philpott_compute_hh_o_pt_(const int & nWater, const int & nPt, const double r[], double f[], double & energy, double const box[])
{
    int const n=nWater*9+nPt*3;
    double R[n];
    for (int i=0; i < n; ++i)
        R[i]=r[i] - r[nWater*9+2];
    potential.computeHH_O_Pt_(nWater, nPt, R, f, energy, box, 0);
}

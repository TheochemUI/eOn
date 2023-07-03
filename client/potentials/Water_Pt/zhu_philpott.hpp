/** @file
Potential for water and Platinum
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
*/
#pragma once
#ifndef FORCEFIELDS_ZHU_PHILPOTT_HPP
#define FORCEFIELDS_ZHU_PHILPOTT_HPP
#include "../Water/spce_ccl.hpp"
#include "zhu_philpott_parameters.hpp"

namespace forcefields {
template <class P = zhu_philpott_parameters::Standard>
class ZhuPhilpott : public SpceCcl, private P {
public:
  ZhuPhilpott();
  ZhuPhilpott(double cutoff, double switchingWidth);
  ZhuPhilpott(ZhuPhilpott const &);
  void operator=(ZhuPhilpott const &);
  ~ZhuPhilpott();
  void computeHH_O_Pt_(const int nWater, const int nPt, const double r[],
                       double f[], double &energy, double const b[],
                       bool const fixed[]);
  void computeHH_O_(const int nWater, const double r[], double f[],
                    double &energy, double const b[], bool const fixed[]);
  int nPlatinum() const;
  void setPlatinum(int nPlatinum, double const positions[]);
  static char const *getName();

private:
  template <int H, int O, int H3, int O3>
  void computeTemplate(const int nWater, const double (*const rh1)[H3],
                       const double (*const rh2)[H3],
                       const double (*const ro)[O3], double (*const fh1)[H3],
                       double (*const fh2)[H3], double (*const fo)[O3],
                       const int nPt, const double rPt[][3], double fPt[][3],
                       double &energy, double const b[],
                       bool const (*const xh1)[H] = 0,
                       bool const (*const xh2)[H] = 0,
                       bool const (*const xo)[O] = 0, bool const *xPt = 0);
  void interactWithCorePt(Water &water, int const nPt, double const rPt[][3],
                          double fPt[][3], double &energy);
  void interactionPtO(double const R1[], double const R2[], double F1[],
                      double F2[], double &energy);
  void interactionPtH(double const R1[], double const R2[], double F1[],
                      double F2[], double &energy);
  void anisotropic(const double distance[], double force[], double &energy,
                   double const epsilon, double const sigma,
                   double const alpha);
  void isotropic10(double const distance, double &force, double &energy,
                   double const epsilon, double const sigma, double const C10);
  void interactWithImage(Water &w1, Water &w2, double &U);
  void coulombWithCutoff(Water &w1, Water &w2, double &u,
                         double const relativePermittivity);
  void coulombFull(Water &w1, Water &w2, double &U,
                   double const relativePermittivity);
  int nPlatinum_;
  double *positions_, *forces_;
};
template class ZhuPhilpott<zhu_philpott_parameters::Standard>;
template class ZhuPhilpott<zhu_philpott_parameters::Iceland>;
} // namespace forcefields

#endif

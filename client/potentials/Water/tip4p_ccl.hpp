#pragma once
#ifndef FORCEFIELDS_TIP4P_CCL_HPP
#define FORCEFIELDS_TIP4P_CCL_HPP
#include "ccl.hpp"

/** @file
TIP4P potential for water
@author Jean-Claude C. Berthet
@date 2006-2007
University of Iceland
*/

namespace forcefields {
class Tip4p : public Ccl {
public:
  Tip4p();
  Tip4p(double cutoff, double switchingWidth);
  ~Tip4p() {}
  void computeHH_O_(const int nAtoms, const double R[], double F[], double &U,
                    const double b[]);
  void computeHH_O_(const int nAtoms, const double R[], double F[], double &U,
                    const double b[], const bool fixed[]);
  static char const *getName();

private:
  struct Water;
  template <int H, int O>
  void computeTemplate(const int nMolecules, const double (*const rh1)[H * 3],
                       const double (*const rh2)[H * 3],
                       const double (*const ro)[O * 3],
                       double (*const fh1)[H * 3], double (*const fh2)[H * 3],
                       double (*const fo)[O * 3], double &energy,
                       double const b[], bool const (*const xh1)[H] = 0,
                       bool const (*const xh2)[H] = 0,
                       bool const (*const xo)[O] = 0);
  void coulombWithCutoff(Water &w1, Water &w2, double &U);
  void coulombFull(Water &w1, Water &w2, double &U);
  void lennardJonesWithCutoff(Water &w1, Water &w2, double &U);
};
} // namespace forcefields
#endif

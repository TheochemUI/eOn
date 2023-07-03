#pragma once
#ifndef FORCEFIELDS_CCL_HPP
#define FORCEFIELDS_CCL_HPP
#include "potential_base.hpp"

/** @file
Potential CCL Table II for intramolecular interactions in water.
@author Jean-Claude C. Berthet
@date 2008
University of Iceland
*/
namespace forcefields {
class Ccl : public PotentialBase {
public:
  Ccl() {}
  void computeHH_O_(const int nAtoms, const double R[], double F[], double &U,
                    const double b[]);
  void computeHH_O_(const int nAtoms, const double R[], double F[], double &U,
                    const double b[], const bool fixed[]);
  char const *getName() const;

protected:
  Ccl(double cutoff, double switchingWidth);
  void intramolecular(double const rh1[], double const rh2[], double const ro[],
                      double fh1[], double fh2[], double fo[], double &energy);
  /// @f$ \rho @f$ parameter in paper. @see forcefields::Ccl
  struct Rho {
    double _1; ///< rho^1
    double _2; ///< rho^2
    double _3; ///< rho^3
    double n[3];
  };
  /// @f$ \Delta \theta @f$ parameter in paper. @see forcefields::Ccl
  struct Dtheta {
    double _1; ///< Delta theta
    double _2; ///< (Delta theta)^2
    double _3; ///< (Delta theta)^3
    double n1[3];
    /// Same as #n1 but for hydrogen 2.
    double n2[3];
  };
  void initialiseRho(Vector3 const &v, Rho &r);
  void initialiseDtheta(Vector3 const &v1, Vector3 const &v2, Dtheta &dth,
                        double const thetaEquilibrium);
  void intramolecular(Rho const &ro1, Rho const &ro2, Dtheta const &dth,
                      double &energy, double fh1[], double fh2[], double fo[]);
  static double const re_;
  static double const thetae_;

private:
  template <int H, int O>
  void computeTemplate(const int nMolecules, const double (*const rh1)[H * 3],
                       const double (*const rh2)[H * 3],
                       const double (*const ro)[O * 3],
                       double (*const fh1)[H * 3], double (*const fh2)[H * 3],
                       double (*const fo)[O * 3], double &energy,
                       double const b[], bool const (*const xh1)[H] = 0,
                       bool const (*const xh2)[H] = 0,
                       bool const (*const xo)[O] = 0);
  void ro_2(Rho const &ro1, Rho const &ro2, double &energy, double &d1,
            double &d2);
  void ro1_ro2(Rho const &s1, Rho const &s2, double &energy, double &d1,
               double &d2);
  void ro_theta(Rho const &s1, Rho const &s2, Dtheta const &s3, double &energy,
                double &d1, double &d2, double &d3);
  void theta_2(Dtheta const &s3, double &energy, double &d3);

  void ro_3(Rho const &s1, Rho const &s2, double &energy, double &d1,
            double &d2);
  void ro_ro1_ro2(Rho const &s1, Rho const &s2, double &energy, double &d1,
                  double &d2);
  void ro_2_theta(Rho const &s1, Rho const &s2, Dtheta const &s3,
                  double &energy, double &d1, double &d2, double &d3);
  void ro1_ro2_theta(Rho const &s1, Rho const &s2, Dtheta const &s3,
                     double &energy, double &d1, double &d2, double &d3);
  void ro_theta_2(Rho const &s1, Rho const &s2, Dtheta const &s3,
                  double &energy, double &d1, double &d2, double &d3);
  void theta_3(Dtheta const &s3, double &energy, double &d3);

  void ro_4(Rho const &s1, Rho const &s2, double &energy, double &d1,
            double &d2);
  void ro1_ro2_ro_2(Rho const &s1, Rho const &s2, double &energy, double &d1,
                    double &d2);
  void ro1_2_ro2_2(Rho const &s1, Rho const &s2, double &energy, double &d1,
                   double &d2);
  void ro_3_theta(Rho const &s1, Rho const &s2, Dtheta const &s3,
                  double &energy, double &d1, double &d2, double &d3);
  void ro_ro1_ro2_theta(Rho const &s1, Rho const &s2, Dtheta const &s3,
                        double &energy, double &d1, double &d2, double &d3);
  void ro_2_theta_2(Rho const &s1, Rho const &s2, Dtheta const &s3,
                    double &energy, double &d1, double &d2, double &d3);
  void ro1_ro2_theta_2(Rho const &s1, Rho const &s2, Dtheta const &s3,
                       double &energy, double &d1, double &d2, double &d3);
  // ro_theta_3 = 0
  void theta_4(Dtheta const &s3, double &energy, double &d3);
};
} // namespace forcefields
#endif

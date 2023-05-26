/** @file
Basic  tools to write potentials.
@author Jean-Claude C. Berthet
@date 2006-2008
University of Iceland
*/
#include "potential_base.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

/** @class forcefields::PotentialBase
@brief Functions and tools commonly used by potentials.
Contains functions for quadratic retrainst, Lennard-Jones, Coulomb, to manage
the periodic boundaries and interface for cutoff. The system of units used is
(eV, Angstrom, fs, e). Rules to combine different Lennard-Jones parameters were
taken from:\n Unlike Lennard-Jones Parameters for Vapor-Liquid Equilibria,
Thorsten Schnabel, Jadran Vrabec , Hans Hasse, Institut fur Technische
Thermodynamik und Thermische Verfahrenstechnik, Universitat Stuttgart, D-70550
Stuttgart, Germany, http://www.itt.uni-stuttgart.de/~schnabel/CR.pdf.
*/
using namespace std;
using namespace forcefields;

/** Non bond interaction cutoff.
When the distance between two molecules is over getCutoff(), van der Waals and
Coulomb interactions between the two molecules are ignored.
@see getSwitchingWidth().*/
PotentialBase::PotentialBase() : cutoff_(6.5), switchingWidth_(2.0) {
  periods_[0] = 0.0;
  periods_[1] = 0.0;
  periods_[2] = 0.0;
}

PotentialBase::PotentialBase(double cutoff, double switchingWidth)
    : cutoff_(cutoff), switchingWidth_(switchingWidth) {
  if (switchingWidth_ > cutoff_) {
    cerr << "Error: getSwitchingWidth() > getCutoff()" << endl;
    exit(EXIT_FAILURE);
  };
  periods_[0] = 0.0;
  periods_[1] = 0.0;
  periods_[2] = 0.0;
}

double PotentialBase::getCutoff() const { return cutoff_; }

void PotentialBase::setCutoff(double cutoff) { cutoff_ = cutoff; }

/** Width of the switching zone.
In the switching zone the potential is replaced by @f$ V=U\times S(x) @f$ where
@em U is the real potential @em V the switching potential and @em S the
switching function. Variable x=0 at the beginning of the switching zone and x=1
at the end the function S is @f$ S=2x^3-3x^2+1 @f$ which has the following
properties: S(0)=1, S(1)=0, S'(0)=0 and S'(1)=0. The switching when the distance
between the two molecules is getCutoff()-getSwitchingWidth() and ends when the
distance is getCutoff();
@see getCutoff().
@{ */
double PotentialBase::getSwitchingWidth() const { return switchingWidth_; }

void PotentialBase::setSwitchingWidth(double width) { switchingWidth_ = width; }
// @}

/** Minimum image representation.
@return Minimum image representation of @a r. The value returned is @f$
\frac{-period}{2} \le r \le \frac{period}{2} @f$
@param[in] r Cyclic coordinate.
@param[in] period Period of the coordinates.
*/
double PotentialBase::applyPeriodicity0(double r, double const period) {
  if (not isinf(period)) {
    double n = r / period + 0.5;
    // This is slightly more efficient than using function floor().
    int m = (int)n;
    if (n < 0)
      --m;
    r -= m * period;
  };
  return r;
}

/** Minimum image representation.
@return Minimum image representation of @a r. The value returned is @f$
\frac{-period}{2} \le r \le \frac{period}{2} @f$.
@param[in] r      Three-dimension vector.
@param[in] periods      Three-dimension vector containing periods along each of
the axes.
*/
void PotentialBase::applyPeriodicity0(double r[], double const periods[]) {
  for (int i = 0; i < 3; ++i)
    r[i] = applyPeriodicity0(r[i], periods[i]);
}

/** Calculate the base of an isosceles triangle.
An isosceles triangle is a triangle with two sides of equal length. If you know
the length of the two equal sides and the angle between the two, the function
can calculate the length of the third side.
@param[in]  length      Length of the two equal sides of the triangle.
@param[in]  angle       Angle between the two equal side in radiant.
*/
double PotentialBase::isoscelesBase(double length, double angle) {
  return 2.0 * length * sin(angle / 2.0);
}

/** Compute quadratic restraints between two atoms.
@f$ U=k(r-r_0)^2 @f$
@param[in]              R1[]    Positions of atom 1.
@param[in]              R2[]    Positions of atom 2.
@param[in,out]       F1[]    Incremented by force on atom 1.
@param[in,out]       F2[]    Incremented by force on atom 2.
@param[in,out]       u      Incremented by energy.
@param[in]              k      Potential curvature
@param[in]              req    Distance of the two atoms at equilibrium.
@warning                 @a F1, @a F2 and @a u are incremented.*/
void PotentialBase::restrainLength(const double R1[], const double R2[],
                                   double F1[], double F2[], double &u,
                                   const double k, const double req) {
  double f, fx;
  double r[3], r1, r2;
  distance(R2, R1, r, r1, r2);
  double d = r1 - req;
  f = -2.0 * k * d;
  u += k * d * d;
  fx = f * r[0] / r1;
  F1[0] -= fx;
  F2[0] += fx;
  fx = f * r[1] / r1;
  F1[1] -= fx;
  F2[1] += fx;
  fx = f * r[2] / r1;
  F1[2] -= fx;
  F2[2] += fx;
}

/** Angular quadratic restraints.
@f$ U=k(a-a_0)^2 @f$
@param[in]              "r1[], r2[], r3[]"    Positions of atoms.
@param[in,out]       "f1[], f2[], f3[]"     Force on atoms. The forces are
incremeneted by the function.
@param[in,out]       u      Incremented by energy.
@param[in]              k      Potential curvature
@param[in]              aeq    Angle made by the three atoms 1,2,3 at
equilibrium.
*/
void PotentialBase::restrainAngle(double const r1[], double const r2[],
                                  double const r3[], double f1[], double f2[],
                                  double f3[], double &u, const double k,
                                  const double aeq) {
  double r12[3], r12_1, r12_2, r23[3], r23_1, r23_2, cosa, a, m, d12, d23;
  distance(r2, r1, r12, r12_1, r12_2);
  distance(r3, r2, r23, r23_1, r23_2);
  cosa = -dotProduct(r12, r23) / r12_1 / r23_1;
  a = acos(cosa);
  u += k * (a - aeq) * (a - aeq);
  m = -2 * k * (a - aeq);          // moment
  m /= -1 * sqrt(1 - cosa * cosa); //*d(a)/d(cos a)
  for (int i = 0; i < 3; ++i) {
    d12 = -cosa * r12[i] / r12_2 - r23[i] / r12_1 / r23_1; // d(cos a)/d(r12)
    d23 = -cosa * r23[i] / r23_2 - r12[i] / r12_1 / r23_1; // d(cos a)/d(r23)
    f1[i] -= m * d12;
    f2[i] += m * (d12 - d23);
    f3[i] += m * d23;
  };
}

const double PotentialBase::ONE_OVER_4_PI_EPSILON0 =
    14.399644532010862; // eV e^2 / Angstrom

/** Calculate centre of two points.
Calculate the centre of @a r1 and @a r2 and store the result in @a rc.
@note The centre is calculated for the image or @a r2 the closed to @a r1.
*/
void PotentialBase::calculateCentre(double const r1[], double const r2[],
                                    double rc[]) {
  for (int i = 0; i < 3; ++i) {
    rc[i] = (r1[i] + unBreak1(r2[i], r1[i], i)) * 0.5;
  };
}

/// Calculate centre of three points. @see calculateCentre().
void PotentialBase::calculateCentre(double const r1[], double const r2[],
                                    double const r3[], double rc[]) {
  for (int i = 0; i < 3; ++i) {
    rc[i] =
        (r1[i] + unBreak1(r2[i], r1[i], i) + unBreak1(r3[i], r1[i], i)) / 3.0;
  };
}

/** Calculate barycentre.
The inverse of this functions the forces is spreadWeightedForce()
@pre Sum of weigths <em> w1, w2, w3 </em> must be 1.
*/
void PotentialBase::calculateWeightedCentre(double const w1, double const w2,
                                            double const w3, double const r1[],
                                            double const r2[],
                                            double const r3[], double rc[]) {
  assert(fabs(w1 + w2 + w3 - 1.0) < 1e-9);
  for (int i = 0; i < 3; ++i) {
    rc[i] = w1 * r1[i] + w2 * unBreak1(r2[i], r1[i], i) +
            w3 * unBreak1(r3[i], r1[i], i);
  };
}

/** Coulomb interaction with single charge based cutoff.
@f$ U = \frac{q_1q_2}{|\mathbf r_1- \mathbf r_2|} @f$
@param[in]              r1	Positions of atom 1.
@param[in]              r2	Positions of atom 2.
@param[in,out]       f1	Incremented by force on atom 1.
@param[in,out]       f2	Incremented by force on atom 2.
@param[in,out]       energy          Incremented by energy.
@param[in]              qq        Product of the charges.
@warning	@a f1, @a f2 and @a u are incremented.
@see getCutoff() and getSwitchingWidth().*/
void PotentialBase::coulombWithCutoff(const double r1[], const double r2[],
                                      double f1[], double f2[], double &energy,
                                      double qq) {
  double d1, d[3];
  // d: vector distance (r1-r2), d1 distance and d2 squared distance
  distance(r1, r2, d, d1);
  if (d1 < cutoff_) {
    double f, en;
    coulomb(d1, f, en, qq);
    switching(d1, f, en);
    for (int l = 0; l < 3; ++l) {
      f1[l] += f * d[l] / d1;
      f2[l] -= f * d[l] / d1;
    };
    energy += en;
  }
}

/** Compute Coulomb interaction between two charges.
@f$ U = \frac{q_1q_2}{|\mathbf r_1- \mathbf r_2|} @f$
@param[in]              r1	Positions of atom 1.
@param[in]              r2	Positions of atom 2.
@param[in,out]       f1	Incremented by force on atom 1.
@param[in,out]       f2	Incremented by force on atom 2.
@param[in,out]       energy           Incremented by energy.
@param[in]              qq        Product of the charges.
@warning	@a f1, @a f2 and @a u are incremented.*/
void PotentialBase::coulomb(const double r1[], const double r2[], double f1[],
                            double f2[], double &energy, double qq) {
  double d1, d[3], f, en;
  // d: vector distance (r1-r2), d1 distance and d2 squared distance
  distance(r1, r2, d, d1);
  coulomb(d1, f, en, qq);
  for (int l = 0; l < 3; ++l) {
    f1[l] += f * d[l] / d1;
    f2[l] -= f * d[l] / d1;
  };
  energy += en;
}

/** Compute Coulomb interaction between two charges.
@f$ U = \frac{q_1q_2}{|\mathbf r_1- \mathbf r_2|} @f$
@param[in]       distance between two charges
@param[out]       force	Force Resulting from the coulom interaction.
@param[out]       energy           Potential energy.
@param[in]         qq        Product of the charges.*/
void PotentialBase::coulomb(double distance, double &force, double &energy,
                            double const qq) {
  energy = ONE_OVER_4_PI_EPSILON0 * qq / distance;
  force = energy / distance;
}

/** Spread force on centre to other points.
Spread force @a fc to @a f1, @a f2. Each receives a half of the force.
*/
void PotentialBase::spreadForce(double f1[], double f2[], double const fc[]) {
  for (int i = 0; i < 3; ++i) {
    double const f = fc[i] / 2;
    f1[i] += f;
    f2[i] += f;
  };
}

/** Spread force on centre to other points.
Spread force @a fc to @a f1, @a f2, @a f3. Each receives 1/3 of the force.
*/
void PotentialBase::spreadForce(double f1[], double f2[], double f3[],
                                double const fc[]) {
  for (int i = 0; i < 3; ++i) {
    double const f = fc[i] / 3;
    f1[i] += f;
    f2[i] += f;
    f3[i] += f;
  };
}

/** Spread force on barycentre to atoms.
This is the inverse functions calculateWeightedCentre()
@pre Sum of weigths <em> w1, w2, w3 </em> must be 1.
@param[in] "w1, w2, w3" weights.
@param[in,out] "f1, f2, f3" vector to spread @a fc to.
@param[in,out] fc force to spread.
*/
void PotentialBase::spreadWeightedForce(double const w1, double const w2,
                                        double const w3, double f1[],
                                        double f2[], double f3[],
                                        double const fc[]) {
  for (int i = 0; i < 3; ++i) {
    double const f = fc[i];
    f1[i] += f * w1;
    f2[i] += f * w2;
    f3[i] += f * w3;
  };
}

/** Conversion for Lennard-Jones.
The Lennard-Jones is sometimes expressed as:
@f[ \frac{A}{r^{12}}-\frac{B}{r^6} @f]
whereas the functions lennardJones() use the expression:
@f[
    4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^6\right]
    @f]
The function calculates the parameter @a epsilon.
*/
double PotentialBase::epsilon(double const A, double const B) {
  return B * B / A / 4.0;
}

/**  Lennard-Jones 12-6 between two atoms.
Return the @a energy and scalar @a force resulting from a Lennard-Jones
potential.
@f$ E=4*\epsilon*(x^{12}-x^6) @f$ with @f$ x=\frac{\sigma}{r} @f$ and @f$
F=-\frac{dE}{dr} @f$ .
@param[in]  distance     Distance between two atoms.
@param[out] force         Force.
@param[out] energy      Potential energy.
@param[in] "epsilon, sigma"   Lennard-Jones parameters.
@see sigma() and epsilon()*/
void PotentialBase::lennardJones(double distance, double &force, double &energy,
                                 double const epsilon, double const sigma) {
  double x = sigma * sigma / (distance * distance);
  x *= x * x;
  energy = 4 * epsilon * (x - 1) * x;
  force = 24 * epsilon * (2 * x - 1) * x / distance;
}

/** Lennard-Jones 12-6 between two atoms.
@param[in]              "R1, R2"	Positions of atom 1 and 2.
@param[in,out]       "F1, F2"	Force on atom 1 and 2.
@param[in,out]       E                    Energy.
@param[in]  "epsilon, sigma"  Parameters of the Lennard-Jones potential.
@warning	@a F1, @a F2 and @a u are incremented.
@see @ref lennardJones(double const, double &, double &, double const, double
const) "lennardJones(distance, ...)".*/
void PotentialBase::lennardJones(const double R1[], const double R2[],
                                 double F1[], double F2[], double &E,
                                 double const epsilon, double const sigma) {
  // WARNING: F1 and F2 are incremented.
  double F, Fx, en;
  double R12[3];
  double L, L2;
  distance(R2, R1, R12, L, L2);
  lennardJones(L, F, en, epsilon, sigma);
  E += en;
  Fx = F * R12[0] / L;
  F1[0] -= Fx;
  F2[0] += Fx;
  Fx = F * R12[1] / L;
  F1[1] -= Fx;
  F2[1] += Fx;
  Fx = F * R12[2] / L;
  F1[2] -= Fx;
  F2[2] += Fx;
}

/** Lennard Jones 12-6 Potential with cutoff.
Interaction between two atoms.
Cutoff and switching width should be set with setCutoff() and
setSwitchingWidth().
@param[in]  "r1, r2"     Positions of the two atoms.
@param[in,out] "f1, f2" Forces on the two atoms.
@param[in,out]    energy      Potential energy.
@param[in]  "epsilon, sigma" See @ref lennardJones(double const distance, double
& force, double & energy, double const epsilon, double const sigma)
"lennardJones()" for definition of @a sigma and @a esplion.
@see @ref lennardJones(double const, double &, double &, double const, double
const) "lennardJones()".
*/
void PotentialBase::lennardJonesWithCutoff(double const r1[], double const r2[],
                                           double f1[], double f2[],
                                           double &energy, double const epsilon,
                                           double const sigma) {
  double r12[3] = {0}, d = 0.0, f, en; // distance
  distance(r2, r1, r12, d);
  if (d < cutoff_) {
    lennardJones(d, f, en, epsilon, sigma);
    switching(d, f, en);
    for (int i = 0; i < 3; ++i) {
      f1[i] -= f * r12[i] / d;
      f2[i] += f * r12[i] / d;
    };
    energy += en;
  };
}

/** Conversion for Lennard-Jones.
The Lennard-Jones is sometimes expressed as:
@f[ \frac{A}{r^{12}}-\frac{B}{r^6} @f]
whereas the functions lennardJones() use the expression:
@f[
    4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^6\right]
    @f]
The function calculates the parameter @a sigma.
*/
double PotentialBase::sigma(double const A, double const B) {
  return pow(A / B, 1.0 / 6.0);
}

/** Smith and Kong combination rules.
@f[
    \epsilon^{SK}_{12}=
    \frac{2^{13}\epsilon_1\sigma_1^6\epsilon_2\sigma_2^6}{\left[(\epsilon_1\sigma_1^{12})^\frac{1}{13}+(\epsilon_1\sigma_1^{12})^\frac{1}{13}\right]^{13}}
    @f]
@see smithKongSigma()*/
double PotentialBase::smithKongEpsilon(double sigma1, double epsilon1,
                                       double sigma2, double epsilon2) {
  double sa6, sb6, n, d;
  sa6 = pow(sigma1, 6.0);
  sb6 = pow(sigma2, 6.0);
  n = pow(2.0, 13.0) * epsilon1 * sa6 * epsilon2 * sb6;
  d = pow(epsilon1 * sa6 * sa6, 1.0 / 13.0) +
      pow(epsilon2 * sb6 * sb6, 1.0 / 13.0);
  d = pow(d, 13);
  return n / d;
}

/** Smith and Kong combination rules.
Combination of unlike Lennard-Jones parameters
@f[
    \sigma^{SK}=
    \left[
        \frac{(\epsilon_1\sigma_1^{12})^{\frac{1}{13}}+(\epsilon_2\sigma_2^{12})^{\frac{1}{13}}}
        {2^{13}\sqrt{\epsilon_1\sigma_1^6\epsilon_2\sigma_2^6}}
        \right]^\frac{1}{6}
    @f]
*/
double PotentialBase::smithKongSigma(double sigma1, double epsilon1,
                                     double sigma2, double epsilon2) {
  double const sa6 = pow(sigma1, 6.0);
  double const sb6 = pow(sigma2, 6.0);
  double n, d;
  n = pow(epsilon1 * sa6 * sa6, 1.0 / 13.0) +
      pow(epsilon2 * sb6 * sb6, 1.0 / 13);
  n = pow(n, 13.0);
  d = pow(2.0, 13.0) * sqrt(epsilon1 * sa6 * epsilon2 * sb6);
  return pow(n / d, 1.0 / 6.0);
}

/** Smooth cutoff switch off.
Should be called when the distance between the two atoms is in the switching
zone (see getSwitchinWidth() for more explanation).
@param[in] distance Distance between two atoms.
@param[in,out] force On atoms. The input values should be the full force for the
interaction between the two atoms. The function returns the corrected force.
@param[in,out] energy Potential energy. The input should be the full energy for
the interaction between the two atoms. The functions returns the corrected
energy.
@note The function is for atom based cutoff only.
@see getSwitchinWidth().
*/
void PotentialBase::switching(double const distance, double &force,
                              double &energy) {
  double x = (distance - cutoff_ + switchingWidth_) / switchingWidth_;
  if (x >= 1.0) { // if more: cutoff
    force = 0.0;
    energy = 0.0;
  }
  if (x > 0.0) { // if more: switching zone
    double S = (2 * x - 3) * x * x + 1;
    double dS = 6 * x * (x - 1);
    double const u = energy;
    force = force * S - u * dS / switchingWidth_;
    energy *= S;
  };
  // else: full interaction, leave force and energy as the are.
}

/** Smooth cutoff switch off.
Should be called when the distance between the two atoms is in the switching
zone (see getSwitchinWidth() for more explanation).
@param[in] "r1, r2" Positions of atoms.
@param[in,out] "f1, f2" Forces on atoms. The input values should be the full
forces for the interaction between the two atoms. The function returns the
corrected forces.
@param[in,out] energy Potential energy. The input should be the full energy for
the interaction between the two atoms. The functions returns the corrected
energy.
@note The function is for atom based cutoff only.
@see getSwitchinWidth().
*/
void PotentialBase::switching(double const r1[], double const r2[], double f1[],
                              double f2[], double &energy) {
  double z[3], z1;
  distance(r1, r2, z, z1);
  double x = (z1 - cutoff_ + switchingWidth_) / switchingWidth_;
  assert(x > 0.0);
  assert(x < 1.0);
  double S = (2 * x - 3) * x * x + 1;
  double dS = 6 * x * (x - 1);
  double const u = energy;
  for (int k = 0; k < 3; ++k) {
    f1[k] = f1[k] * S - u * dS / switchingWidth_ * z[k] / z1;
    f2[k] = f2[k] * S + u * dS / switchingWidth_ * z[k] / z1;
  };
  energy *= S;
}

/** Minimum image representation.
@return Minimum image representation of @a r. The value returned is @f$
\frac{-period}{2} \le r \le \frac{period}{2} @f$. The periods should be set with
setPeriodicity();
@param[in] r      Three-dimension vector.
@see setPeriodicity().*/
void PotentialBase::applyPeriodicity1(double r[]) {
  applyPeriodicity0(r, periods_);
}

/** Minimum image representation.
@return Minimum image representation of @a r. The value returned is @f$
\frac{-period}{2} \le r \le \frac{period}{2} @f$. The periods should be set with
setPeriodicity();
@param[in] r Cyclic coordinate.
@param[in] axis   Either 0, 1 or 2.
@see setPeriodicity().*/
double PotentialBase::applyPeriodicity1(double r, int const axis) {
  return applyPeriodicity0(r, periods_[axis]);
}

/** Distance vector, norm and norm square.
@param[in] "x, y" Two three dimension vector.
@param[out] z     @f$ \mathbf x -\mathbf  y @f$ with pbc applied.
@param[out] z1    @f$ |\mathbf z | @f$
@param[out] z2    @f$ |\mathbf z |^2 @f$
@see setPeriodicity().*/
void PotentialBase::distance(const double x[], const double y[], double z[],
                             double &z1, double &z2) {
  assert(x);
  assert(y);
  assert(z);
  for (int i = 0; i < 3; ++i)
    z[i] = x[i] - y[i];
  applyPeriodicity1(z);
  z2 = z[0] * z[0] + z[1] * z[1] + z[2] * z[2];
  z1 = sqrt(z2);
}

/** Distance.
@param[in] "x, y" Two three dimension vector.
@param[out] z[]
@param[out] z1     @f$| \mathbf x -\mathbf y | @f$ with pbc applied.
@see setPeriodicity().*/
void PotentialBase::distance(const double x[], const double y[], double z[],
                             double &z1) {
  double z2;
  distance(x, y, z, z1, z2);
}

/** Distance vector, norm.
@see distance(const double x[], const double y[], double z[], double & z1,
double & z2) and setPeriodicity().*/
void PotentialBase::distance(const double x[], const double y[], double z[]) {
  double z1;
  distance(x, y, z, z1);
}
/// Distance. @see distance
void PotentialBase::distance(const double x[], const double y[], double &z1) {
  double z[3];
  distance(x, y, z, z1);
}

/** Distance vector, norm and norm square.
@param[in] "x, y" Two three dimension vector.
@param[out] z     Return x-y, |x-y| and |x-y|^2.
@see setPeriodicity().*/
void PotentialBase::distance(const double x[], const double y[], Vector3 &z) {
  distance(x, y, z.v, z._1, z._2);
}

/** Set periodicity.
Must be used before calling certain functions: distance(), applyPeriodicity(),
etc...
*/
void PotentialBase::setPeriodicity(const double periods[]) {
  periods_[0] = periods[0];
  periods_[1] = periods[1];
  periods_[2] = periods[2];
}

/** Undo the separation of two atoms created by the periodic boundaries.
Example: two atoms H and O are at respectively x = -9.5 and x = 9.5 and the
period along @em x is 20. The real distance between the two atoms is 1. We use
unBreak() on the coordinates of these two atoms:
@code
unBreak(rh, ro);
@endcode
The coordinates @em rh of H is now 10.5 thus making O and H appears at a
distance of 1 of each other.\n Apply periodic boundaries condition to vector @f$
r - r_{ef} @f$.
@param[in]    r    Returns @f$ \textrm{applyPeriodicity}(r - r_{ef}) +r_{ef} @f$
@param[in] ref
@param[in] period
*/
double PotentialBase::unBreak0(double const r, double const ref,
                               double const period) {
  return applyPeriodicity0(r - ref, period) + ref;
}
/// Undo the separation of two atoms created by the periodic boundaries. @see
/// unBreak0()
void PotentialBase::unBreak0(double r[], double const ref[],
                             double const periods[]) {
  for (int a = 0; a < 3; ++a)
    r[a] = unBreak0(r[a], ref[a], periods[a]);
}
/// Undo the separation of two atoms created by the periodic boundaries. @see
/// unBreak0()
void PotentialBase::unBreak1(double r[], double const ref[]) {
  unBreak0(r, ref, periods_);
}
/// Undo the separation of two atoms created by the periodic boundaries. @see
/// unBreak0()
double PotentialBase::unBreak1(double const r, double const ref,
                               int const axis) {
  return unBreak0(r, ref, periods_[axis]);
}

/** Potential for Platinum.
Simple 12-6 Lennard-Jones potential for Platinum.\n
Parameters from: C. Kittel, Introduction to Solid State Physics, (Wiley, New
York, 1986).
@param[in]  nAtoms      Number of atoms.
@param[in]  positions     Positions of atoms.
@param[out] forces     Forces on atoms.
@param[out] energy      Potential energy.
@param[in]  periods     Periods for periodic boundary conditions.
@param[in] fixed  Array of boolean of length equal to @a nAtoms. The array shall
contain @em true is the atom is fixed, @em false is movable.
*/
void PotentialBase::computePt(int const nAtoms, double positions[],
                              double forces[], double &energy,
                              double const periods[], bool const fixed[]) {
  int const nCoord = 3 * nAtoms;
  for (int i = 0; i < nCoord; ++i)
    forces[i] = 0.0;
  energy = 0.0;
  double const(*r)[3] = reinterpret_cast<double const(*)[3]>(positions);
  double(*f)[3] = reinterpret_cast<double(*)[3]>(forces);
  setPeriodicity(periods); // Set Periodic boundaries. Essential in order for
                           // some functions to work.
  for (int i = nAtoms - 1; i > 0; --i) {
    for (int j = i - 1; j >= 0; --j) {
      if (not fixed[i] and not fixed[j]) {
        lennardJonesWithCutoff(r[i], r[j], f[i], f[j], energy, EPSILON_PT,
                               SIGMA_PT);
      };
    };
  };
}

/** Platinum Lennard-Jones.
Lennard-Jones parameters for platinum.
@see computePt().
@{ */
const double PotentialBase::EPSILON_PT = 0.68165797577788501; // eV
const double PotentialBase::SIGMA_PT = 2.54;                  // Angstrom
/// @}
/*
 const double PotentialBase::ONE_OVER_4_PI_EPSILON0 =
 unit_system::ONE_OVER_4_PI_EPSILON0; const double
 PotentialBase::EPSILON_PT=65.77*KJ_PER_MOL; const double
 PotentialBase::SIGMA_PT=0.254*NM;
 */

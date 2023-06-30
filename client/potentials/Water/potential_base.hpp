#pragma once
#ifndef FORCEFIELDS_POTENTIAL_BASE_HPP
#define FORCEFIELDS_POTENTIAL_BASE_HPP

#include <cassert>
#include <cmath>

/** @file
Basic  tools to write potentials.
@author Jean-Claude C. Berthet
@date 2006-2008
University of Iceland
*/
namespace forcefields {
class PotentialBase {
public:
  PotentialBase();
  PotentialBase(double cutoff, double switchingWidth);
  virtual ~PotentialBase(){};
  double getCutoff() const;
  void setCutoff(double cutoff);
  double getSwitchingWidth() const;
  void setSwitchingWidth(double width);
  static double applyPeriodicity0(double r, double const period);
  static void applyPeriodicity0(double r[], double const periods[]);

protected:
  //----------------------------------------Bond Pair
  // interactions--------------------------------------------------
  /// @name Bond Interactions
  /// @{
  static double isoscelesBase(double length, double angle);
  void restrainLength(const double R1[], const double R2[], double F1[],
                      double F2[], double &u, const double k, const double r0);
  /// @}
  /// @name Angle
  //@{
  void restrainAngle(double const r1[], double const r2[], double const r3[],
                     double f1[], double f2[], double f3[], double &u,
                     const double k, const double aeq);
  //@}
  //----------------------------------------coulomb--------------------------------------------------
  /// @name Coulomb
  /// @{
  static const double ONE_OVER_4_PI_EPSILON0;
  void calculateCentre(double const r1[], double const r2[], double rc[]);
  void calculateCentre(double const r1[], double const r2[], double const r3[],
                       double rc[]);
  void calculateWeightedCentre(double const w1, double const w2,
                               double const w3, double const r1[],
                               double const r2[], double const r3[],
                               double rc[]);
  void coulombWithCutoff(const double r1[], const double r2[], double f1[],
                         double f2[], double &u, double const qq);
  void coulomb(const double r1[], const double r2[], double f1[], double f2[],
               double &u, double const qq);
  void coulomb(double distance, double &force, double &energy, double const qq);
  void spreadForce(double f1[], double f2[], double const fc[]);
  void spreadForce(double f1[], double f2[], double f3[], double const fc[]);
  void spreadWeightedForce(double const w1, double const w2, double const w3,
                           double f1[], double f2[], double f3[],
                           double const fc[]);
  template <int N, class R = double (*const)[N], class F = double (*const)[N]>
  struct ChargeGroup {
    double const *const centre_;
    R position_;
    F force_;
  };
  template <int N, class R, class F>
  void addForces(ChargeGroup<N, R, F> const &g1, ChargeGroup<N, R, F> &g2);
  template <int N, class R, class F>
  void switching(ChargeGroup<N, R, F> &g1, ChargeGroup<N, R, F> &g2,
                 double &energy, double cutoff, double switchingWidth);
  //-------------------------------------------Lennard-Jones---------------------------------------------------------
  /// @name Lennard-Jones
  /// @{
  static double epsilon(double const A, double const B);
  void lennardJones(double const distance, double &force, double &energy,
                    double const epsilon, double const sigma);
  void lennardJones(const double R1[], const double R2[], double F1[],
                    double F2[], double &E, double const epsilon,
                    double const sigma);
  void lennardJonesWithCutoff(double const r1[], double const r2[], double f1[],
                              double f2[], double &energy, double const epsilon,
                              double const sigma);
  static double sigma(double const A, double const B);
  static double smithKongEpsilon(double sigma1, double epsilon1, double sigma2,
                                 double epsilon2);
  static double smithKongSigma(double sigma1, double epsilon1, double sigma2,
                               double epsilon2);
  void switching(double const distance, double &force, double &energy);
  void switching(double const r1[], double const r2[], double f1[], double f2[],
                 double &energy);
  /// @}
  //------------------------------------Periodicity and
  // distances---------------------------------------------
  /// @name Periodic Boundaries, Distance...
  /// @{
  /// 3-D vector.
  struct Vector3 {
    double v[3]; ///< Vector.
    double _1;   ///< Vector's norm.
    double _2;   ///< Vector's norm square.
  };
  double applyPeriodicity1(double r, int const axis);
  void applyPeriodicity1(double r[]);
  void distance(const double x[], const double y[], double z[], double &z1,
                double &z2);
  void distance(const double x[], const double y[], double z[], double &z1);
  void distance(const double x[], const double y[], double z[]);
  void distance(const double x[], const double y[], double &z1);
  void distance(const double x[], const double y[], Vector3 &z);
  void setPeriodicity(const double periods[]);
  static double unBreak0(double const r, double const ref, double const period);
  static void unBreak0(double r[], double const ref[], double const periods[]);
  double unBreak1(double const r, double const ref, int const axis);
  void unBreak1(double r[], double const ref[]);
  /// @}
  //----------------------------------Common potential for platinum
  //--------------------------------------
  void computePt(int const nAtoms, double positions[], double forces[],
                 double &energy, double const periods[], bool const fixed[]);
  static const double EPSILON_PT;
  static const double SIGMA_PT;
  //-----------------------------------------3D vector
  // functions---------------------------------------------------
  /** @defgroup vector3d Common operations 3D vectors.
      Implementation is in header file for better optimisaition.
      @{ */
  static inline void divide(double v[], double const divisor);
  static inline double *crossProduct(const double v[], const double w[],
                                     double e[]);
  static inline double dotProduct(double const v[], double const w[]);
  static inline void multiply(double v[], double const factor);
  static inline double norm(double const v[]);
  static inline void normalise(double v[]);
  /// @}
  double cutoff_;
  double periods_[3];
  double switchingWidth_;
};
} // namespace forcefields

/** Increment forces.
Increment forces in @a g2 by forces in @a g1.
*/
template <int N, class R, class F>
void forcefields::PotentialBase::addForces(ChargeGroup<N, R, F> const &g1,
                                           ChargeGroup<N, R, F> &g2) {
  increment(g2.force_, g1.force_);
}

template <int N, class R, class F>
void forcefields::PotentialBase::switching(ChargeGroup<N, R, F> &g1,
                                           ChargeGroup<N, R, F> &g2,
                                           double &energy, double cutoff,
                                           double switchingWidth) {
  double z[3], z1;
  distance(g1.centre_, g2.centre_, z, z1);
  assert(z1 <= cutoff);
  assert(z1 >= cutoff - switchingWidth);
  double x = (z1 - cutoff + switchingWidth) / switchingWidth;
  double S = (2 * x - 3) * x * x + 1;
  double dS = 6 * x * (x - 1);
  double const u = energy;
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < 3; ++k) {
      g1.force_[i][k] =
          g1.force_[i][k] * S - u * dS / switchingWidth * z[k] / z1 / N;
      g2.force_[i][k] =
          g2.force_[i][k] * S + u * dS / switchingWidth * z[k] / z1 / N;
    }
  };
  energy *= S;
}

/** Cross product.
@f$ \mathbf e=\mathbf v \times \mathbf w @f$
@param[in] v      3D vector.
@param[in]  w     3D vector.
@param[out] e     Result.
*/
double *forcefields::PotentialBase::crossProduct(const double v[],
                                                 const double w[], double e[]) {
  e[0] = v[1] * w[2] - v[2] * w[1];
  e[1] = v[2] * w[0] - v[0] * w[2];
  e[2] = v[0] * w[1] - v[1] * w[0];
  return e;
}

/** Divide vector.
@param[in,out] v Three-dimension vector to divide.
@param[in] divisor Divide @a v by @a divisor.
*/
void forcefields::PotentialBase::divide(double v[], double const divisor) {
  v[0] /= divisor;
  v[1] /= divisor;
  v[2] /= divisor;
}

/** Dot product.
@f$ \mathbf v . \mathbf w @f$
@param[in] v      3D vector.
@param[in]  w     3D vector.
@return Dot product.
*/
double forcefields::PotentialBase::dotProduct(double const v[],
                                              double const w[]) {
  return v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
}

/** Multiply vector.
@param[in,out] v Three-dimension vector to multiply.
@param[in] factor Multiply @a v by @a factor.
*/
void forcefields::PotentialBase::multiply(double v[], double const factor) {
  v[0] *= factor;
  v[1] *= factor;
  v[2] *= factor;
}

/** Norm of 3D vector.
@f$ |\mathbf v| @f$
@param[in] v      3D vector.
@return Norm.
*/
double forcefields::PotentialBase::norm(double const v[]) {
  return sqrt(dotProduct(v, v));
}

/** Normalise 3D vector
@f$ \frac{\mathbf v}{ |\mathbf v|} @f$
@param[in,out] v  3D vector to normalise.
*/
void forcefields::PotentialBase::normalise(double v[]) {
  double const n = norm(v);
  divide(v, n);
}

#endif

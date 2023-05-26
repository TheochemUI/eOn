/** @file
Potential for water and Platinum
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
*/
#include "zhu_philpott.hpp"
#include <cassert>
#include <cmath>

/** @class forcefields::ZhuPhilpott
@brief Forcefield for water and platinum interactions.
This forcefield is the A2 water-platinum potential invented by @ref zhu1994 "Zhu
and Philpott". It includes the SPC/E+CCL potential for interaction between the
molecules of water. The SPC/E is a potential for constrained water. This
implementation includes restraints so it may be used without constraint (see
class SpceCcl).\n The potential uses Kong's rules to combine the Lennard-Jones
parameters of platinum with oxygen and hydrogen. A @ref combinationReview
"review of combination rules" including Kong's rule can be find on the web.\n
The system of unit used by the class is (eV, Angstrom, fs, e).
@section references References
@anchor zhu1994
Interaction of water with metal surfaces, S.-B. Zhu and M.R. Philpott J. Chem.
Phys. (1994) vol. 100, No 9,  p. 6961.\n
@anchor kimura2003
A Molecular Dynamics Simulation of Water Droplet in Contact with a Platinum
Surface, T. Kimura and S. Maruyama, University of Tokyo.\n
@anchor combinationReview
Unlike Lennard-Jones Parameters for Vapor-Liquid Equilibria, Thorsten Schnabel,
Jadran Vrabec , Hans Hasse, Institut fur Technische Thermodynamik und Thermische
Verfahrenstechnik, Universitat Stuttgart, D-70550 Stuttgart, Germany,
http://www.itt.uni-stuttgart.de/~schnabel/CR.pdf.
@internal The class contains many functions with some of them calculating only
one type of interaction. Those functions often have a parameter called @a force
and one @a energy . Some of those functions simply return the force and energy
into the variables provided while others add the force and energy to the values
already stored in the variables. Look at the documentation of the function to
know how it operates. The parameter is marked as both input and ouput then it
adds the values, otherwise it is simply marked as output.
*/

// #define DEBUG_LEVEL 5
#ifdef DEBUG_LEVEL
#warning zhu_philpott debug: DEBUG_LEVEL defined
#define DEBUG_LEVEL_RETURN(x)                                                  \
  if (x == DEBUG_LEVEL)                                                        \
    continue;
#else
#define DEBUG_LEVEL_RETURN(x)
#endif

using namespace std;
using namespace forcefields;

template <class P> ZhuPhilpott<P>::ZhuPhilpott() : SpceCcl() {
  nPlatinum_ = 0;
  positions_ = 0;
  forces_ = 0;
}

template <class P>
ZhuPhilpott<P>::ZhuPhilpott(double cutoff, double switchingWidth)
    : SpceCcl(cutoff, switchingWidth) {
  nPlatinum_ = 0;
  positions_ = 0;
  forces_ = 0;
}

template <class P> ZhuPhilpott<P>::ZhuPhilpott(ZhuPhilpott const &zhuPhilpott) {
  operator=(zhuPhilpott);
}

template <class P>
void ZhuPhilpott<P>::operator=(ZhuPhilpott const &zhuPhilpott) {
  setPlatinum(zhuPhilpott.nPlatinum_, zhuPhilpott.positions_);
}

template <class P> ZhuPhilpott<P>::~ZhuPhilpott() {
  delete[] positions_;
  delete[] forces_;
}

/** Compute water-platinum forcefield.
@param[in] nWater Number of molecules of water.
@param[in] nPt    Number of platinum atoms.
@param[in] r      Positions of the atoms.
@param[out] f     Forces on the atoms.
@param[out] energy Potential energy.
@param[in]  b     Periodic boundaries. Length: 3
@param[in]  fixed True for fixed atoms, false otherwise.
The order of the atoms is all hydrogens (hydrogens belonging to the same
molecule are next to each others), all oxygen, all platinum. The length of
arrays @a r and @a f is @f$ 3\times(3 \times nWater+nPt) @f$. The length of
array @a fixed is @f$ 3 \times nWater+nPt @f$. When @a fixed is provided the
interaction between two fixed atoms may be skipped to speed up the calculation.
The parameter is optional.
@warning  The conductivity is simulated by interaction a mirror images of
charges. The mirror surface is set at z=0, which should correspond to the
location of the highest layer of platinum atoms.
*/
template <class P>
void ZhuPhilpott<P>::computeHH_O_Pt_(const int nWater, const int nPt,
                                     const double r[], double f[],
                                     double &energy, double const b[],
                                     bool const fixed[]) {
  const double(*const rh1)[6] = reinterpret_cast<const double(*)[6]>(r);
  const double(*const rh2)[6] = reinterpret_cast<const double(*)[6]>(&r[3]);
  const double(*const ro)[3] =
      reinterpret_cast<const double(*)[3]>(&r[nWater * 6]);
  const double(*const rPt)[3] =
      reinterpret_cast<const double(*)[3]>(&r[nWater * 9]);

  double(*const fh1)[6] = reinterpret_cast<double(*)[6]>(f);
  double(*const fh2)[6] = reinterpret_cast<double(*)[6]>(&f[3]);
  double(*const fo)[3] = reinterpret_cast<double(*)[3]>(&f[nWater * 6]);
  double(*const fPt)[3] = reinterpret_cast<double(*)[3]>(&f[nWater * 9]);

  if (fixed) {
    bool const(*const xh1)[2] = reinterpret_cast<bool const(*)[2]>(fixed);
    bool const(*const xh2)[2] = reinterpret_cast<bool const(*)[2]>(&fixed[1]);
    bool const(*const xo)[1] =
        reinterpret_cast<bool const(*)[1]>(&fixed[nWater * 2]);
    bool const *const xPt = &fixed[nWater * 3];
    computeTemplate(nWater, rh1, rh2, ro, fh1, fh2, fo, nPt, rPt, fPt, energy,
                    b, xh1, xh2, xo, xPt);
  } else {
    bool const(*const xh1)[2] = 0;
    bool const(*const xh2)[2] = 0;
    bool const(*const xo)[1] = 0;
    computeTemplate(nWater, rh1, rh2, ro, fh1, fh2, fo, nPt, rPt, fPt, energy,
                    b, xh1, xh2, xo, 0);
  };
}

/** Compute water-platinum interactions, call with water's positions only.
This function is called with the coordinates of the atoms of water only.
The positions of the platinum atoms must have been previously provided through
function setPlatinum().
@param[in] nWater Number of molecules of water.
@param[in] r      Positions of atoms in water.
@param[out] f     Forces on the atoms in water.
@param[out] energy Potential energy.
@param[in]  b     Periodic boundaries. Length: 3
@param[in]  fixed True for fixed atoms, false otherwise.
The length of arrays @a r and @a f is @f$ 9\times nWater @f$.
The length of array @a fixed is @f$ 3 \times nWater @f$.
@note The function does not return the forces on platinum atoms as they as
assumed to be fixed.
@see setPlatinum() and computeHH_O_Pt_().
*/
template <class P>
void ZhuPhilpott<P>::computeHH_O_(const int nWater, const double r[],
                                  double f[], double &energy, double const b[],
                                  bool const fixed[]) {
  const double(*const rh1)[6] = reinterpret_cast<const double(*)[6]>(r);
  const double(*const rh2)[6] = reinterpret_cast<const double(*)[6]>(&r[3]);
  const double(*const ro)[3] =
      reinterpret_cast<const double(*)[3]>(&r[nWater * 6]);
  const double(*const rPt)[3] =
      reinterpret_cast<const double(*)[3]>(positions_);

  double(*const fh1)[6] = reinterpret_cast<double(*)[6]>(f);
  double(*const fh2)[6] = reinterpret_cast<double(*)[6]>(&f[3]);
  double(*const fo)[3] = reinterpret_cast<double(*)[3]>(&f[nWater * 6]);
  double(*const fPt)[3] = reinterpret_cast<double(*)[3]>(forces_);

  if (fixed) {
    bool const(*const xh1)[2] = reinterpret_cast<bool const(*)[2]>(fixed);
    bool const(*const xh2)[2] = reinterpret_cast<bool const(*)[2]>(&fixed[1]);
    bool const(*const xo)[1] =
        reinterpret_cast<bool const(*)[1]>(&fixed[nWater * 2]);
    computeTemplate(nWater, rh1, rh2, ro, fh1, fh2, fo, nPlatinum_, rPt, fPt,
                    energy, b, xh1, xh2, xo, 0);
  } else {
    bool const(*const xh1)[2] = 0;
    bool const(*const xh2)[2] = 0;
    bool const(*const xo)[1] = 0;
    computeTemplate(nWater, rh1, rh2, ro, fh1, fh2, fo, nPlatinum_, rPt, fPt,
                    energy, b, xh1, xh2, xo, 0);
  };
}

/// Number of platinum atoms.
/// @see setPlatinum().
template <class P> int ZhuPhilpott<P>::nPlatinum() const { return nPlatinum_; }

/** Initialises the positions of the atoms of platinum.
Use before calling computeHH_O_().
*/
template <class P>
void ZhuPhilpott<P>::setPlatinum(int nPlatinum, double const positions[]) {
  assert(nPlatinum >= 0);
  if (nPlatinum == 0) {
    nPlatinum_ = 0;
    delete[] positions_;
    positions_ = 0;
    delete[] forces_;
    forces_ = 0;
    return;
  };
  int const n = nPlatinum * 3;
  if (nPlatinum_ != nPlatinum) {
    nPlatinum_ = nPlatinum;
    delete[] positions_;
    positions_ = new double[n];
    delete[] forces_;
    forces_ = new double[n];
  };
  for (int i = 0; i < n; ++i) {
    positions_[i] = positions[i];
  };
}

/// Name of the potential.
template <class P> char const *ZhuPhilpott<P>::getName() {
  return "ZhuPhilpott";
}

template <class P>
template <int H, int O, int H3, int O3>
void ZhuPhilpott<P>::computeTemplate(
    const int nWater, const double (*const rh1)[H3],
    const double (*const rh2)[H3], const double (*const ro)[O3],
    double (*const fh1)[H3], double (*const fh2)[H3], double (*const fo)[O3],
    const int nPt, const double rPt[][3], double fPt[][3], double &energy,
    double const b[], bool const (*const xh1)[H], bool const (*const xh2)[H],
    bool const (*const xo)[O], bool const xPt[]) {
  for (int i = 0; i < nWater; ++i) {
    for (int a = 0; a < 3; ++a) {
      fh1[i][a] = 0.0;
      fh2[i][a] = 0.0;
      fo[i][a] = 0.0;
    };
  };
  for (int i = 0; i < nPt; ++i) {
    for (int a = 0; a < 3; ++a)
      fPt[i][a] = 0.0;
  };
  energy = 0.0;
  setPeriodicity(b);
  for (int i = nWater - 1; i >= 0; --i) {
    double centre1[3] = {0};
    calculateCentre(rh1[i], rh2[i], ro[i], centre1);
    Water w1(rh1[i], rh2[i], ro[i], centre1, fh1[i], fh2[i], fo[i]);
    intramolecular(w1, energy);
    assert(not isnan(energy));
    DEBUG_LEVEL_RETURN(0)
    // Two next lines for interaction with platinum
    interactWithCorePt(w1, nPt, rPt, fPt,
                       energy); // called Vw-core in Zhu and Philpott
    assert(not isnan(energy));
    DEBUG_LEVEL_RETURN(1)
    interactWithImage(w1, w1, energy); // called Vw-cond in Zhu and Philpott
    DEBUG_LEVEL_RETURN(2)
    assert(not isnan(energy));
    for (int j = i - 1; j >= 0; --j) {
      bool areFixed = false;
      if (xh1 and xh2 and xo) {
        areFixed = xh1[i][0] and xh2[i][0] and xo[i][0];
        // check if all the atoms of molecule j are fixed.
        areFixed &= xh1[j][0] and xh2[j][0] and xo[j][0];
      };
      // if both molecules are fixed skip force calculation
      if (not areFixed) {
        double centre2[3] = {0};
        calculateCentre(rh1[j], rh2[j], ro[j], centre2);
        Water w2(rh1[j], rh2[j], ro[j], centre2, fh1[j], fh2[j], fo[j]);
        lennardJonesWithCutoff(w1, w2, energy);
        DEBUG_LEVEL_RETURN(3)
        coulombWithCutoff(w1, w2, energy, 1.0);
        DEBUG_LEVEL_RETURN(4)
        // Two next lines for interaction with platinum
        // called Vw-cond in Zhu and Philpott
        interactWithImage(w1, w2, energy); // w1 interacts with image of w2
        interactWithImage(w2, w1, energy); // w2 interacts with image of w1
        DEBUG_LEVEL_RETURN(5)
      };
    };
#ifndef NDEBUG
    for (int a = 0; a < 3; a++) {
      assert(not isnan(fh1[i][a]));
      assert(not isnan(fh2[i][a]));
      assert(not isnan(fo[i][a]));
    };
#endif
  };
  assert(not isnan(energy) and not isinf(energy));
}

/** Interaction of one molecule of water with the whole platinum.
@param[in, out]   water    Molecule of water. The forces are incremented.
@param[in]  nPt   Number of Platinum atom.
@param[in]  rPt   Positions of the platinum.
@param[in, out] fPt     Forces on the platinum. The forces are incremented
@param[in, out] energy Add the energy.
*/
template <class P>
void ZhuPhilpott<P>::interactWithCorePt(Water &water, int const nPt,
                                        double const rPt[][3], double fPt[][3],
                                        double &energy) {
  for (int i = 0; i < nPt; ++i) {
    double r[3], r1; // for distance between two atoms
                     // Pt - O
    distance(water.ro_, rPt[i], r, r1);
    if (r1 <= cutoff_ - switchingWidth_) {
      interactionPtO(water.ro_, rPt[i], water.fo_, fPt[i], energy);
    } else if (r1 < cutoff_) {
      double en = 0;
      double fo[3] = {0}, fPtTmp[3] = {0};
      interactionPtO(water.ro_, rPt[i], fo, fPtTmp, en);
      switching(water.ro_, rPt[i], fo, fPtTmp, en);
      energy += en;
      for (int j = 0; j < 3; j++) {
        water.fo_[j] += fo[j];
        fPt[i][j] += fPtTmp[j];
      };
    };
    // Pt - H1
    distance(water.rh1_, rPt[i], r, r1);
    if (r1 <= cutoff_ - switchingWidth_) {
      interactionPtH(water.rh1_, rPt[i], water.fh1_, fPt[i], energy);
    } else if (r1 < cutoff_) {
      double en = 0;
      double fh[3] = {0}, fPtTmp[3] = {0};
      interactionPtH(water.rh1_, rPt[i], fh, fPtTmp, en);
      switching(water.rh1_, rPt[i], fh, fPtTmp, en);
      energy += en;
      for (int j = 0; j < 3; j++) {
        water.fh1_[j] += fh[j];
        fPt[i][j] += fPtTmp[j];
      };
    };
    // Pt - H2
    distance(water.rh2_, rPt[i], r, r1);
    if (r1 <= cutoff_ - switchingWidth_) {
      interactionPtH(water.rh2_, rPt[i], water.fh2_, fPt[i], energy);
    } else if (r1 < cutoff_) {
      double en = 0;
      double fh[3] = {0}, fPtTmp[3] = {0};
      interactionPtH(water.rh2_, rPt[i], fh, fPtTmp, en);
      switching(water.rh2_, rPt[i], fh, fPtTmp, en);
      energy += en;
      for (int j = 0; j < 3; j++) {
        water.fh2_[j] += fh[j];
        fPt[i][j] += fPtTmp[j];
      };
    };
  };
}

template <class P>
void ZhuPhilpott<P>::interactionPtO(double const R1[], double const R2[],
                                    double F1[], double F2[], double &energy) {
  double R12[3], F12[3] = {0}, r1, f1 = 0;
  distance(R1, R2, R12, r1);
  anisotropic(R12, F12, energy, P::epsilonOPt_, P::sigmaOPt_, P::alpha_);
  isotropic10(r1, f1, energy, P::epsilonOPt_, P::sigmaOPt_, P::C10_O_);
  for (int j = 0; j < 3; j++) {
    F1[j] += F12[j] + f1 * R12[j] / r1;
    F2[j] -= F12[j] + f1 * R12[j] / r1;
  };
}

template <class P>
void ZhuPhilpott<P>::interactionPtH(double const R1[], double const R2[],
                                    double F1[], double F2[], double &energy) {
  double R12[3], F12[3] = {0}, r1, f1 = 0;
  distance(R1, R2, R12, r1);
  anisotropic(R12, F12, energy, P::epsilonHPt_, P::sigmaHPt_, P::alpha_);
  isotropic10(r1, f1, energy, P::epsilonHPt_, P::sigmaHPt_, P::C10_H_);
  for (int j = 0; j < 3; j++) {
    F1[j] += F12[j] + f1 * R12[j] / r1;
    F2[j] -= F12[j] + f1 * R12[j] / r1;
  };
}

/** Anisotropic interaction between water and platinum.
The isotropic iteraction between a atom Pt and an atom O or H of a molecule of
water has the function form:
@f[
    E=4\epsilon\left[
        \left(\frac{\sigma^2}{\alpha^2\rho^2+z^2}\right)^6-\left(\frac{\sigma^2}{\rho^2/\alpha^2+z^2}\right)^3
        \right]
    @f]
@param[in]      distance   Vector distance between to atom.
@param[in, out]      force         Vector force created by the interaction
between the two atoms.
@param[in, out]      energy      Energy created by the interaction.
@param[in] epsilon See equation.
@param[in] sigma See equation.
@param[in] alpha See equation.
@internal The subroutine has not been optimised.*/
template <class P>
void ZhuPhilpott<P>::anisotropic(const double distance[], double force[],
                                 double &energy, double const epsilon,
                                 double const sigma, double const alpha) {
  double z2, z, a, b, A, B, dE_da, dE_db, rho, rho2, alpha2;
  // WARNING: F1 and F2 are incremented.
  alpha2 = alpha * alpha;
  double f;
  rho2 = distance[0] * distance[0] + distance[1] * distance[1];
  rho = sqrt(rho2);
  z = distance[2];
  z2 = z * z;
  a = pow(sigma, 2.0) / (rho2 * alpha2 + z2);
  b = pow(sigma, 2.0) / (rho2 / alpha2 + z2);
  A = pow(a, 6.0);
  B = pow(b, 3.0);
  energy += 4.0 * epsilon * (A - B);
  dE_da = 24.0 * epsilon * pow(a, 5.0);
  dE_db = -12.0 * epsilon * pow(b, 2.0);
  f = 2.0 * (dE_da * pow(a, 2.0) * alpha2 + dE_db * pow(b, 2.0) / alpha2) /
      pow(sigma, 2.0);
  force[0] += f * distance[0];
  force[1] += f * distance[1];
  f = 2 * (dE_da * pow(a, 2.0) + dE_db * pow(b, 2.0)) / pow(sigma, 2.0);
  force[2] += f * distance[2];
}

/** Isotropic interaction between water and platinum.
The isotropic iteraction between a atom Pt and an atom O or H of a molecule of
water has the functional form:
@f[ E=-4\epsilon\ C_{10}\left(\frac{\sigma}{r}\right)^{10} @f]
@param[in]      distance    Distance between two atoms.
@param[in,out]     force   Force create by the interaction between the two
atoms.
@param[in,out]      energy    Energy created by the interaction.
@param[in] epsilon See equation.
@param[in] sigma See equation.
@param[in] C10 See equation.
@internal The subroutine has not been optimised.*/
template <class P>
void ZhuPhilpott<P>::isotropic10(double const distance, double &force,
                                 double &energy, double const epsilon,
                                 double const sigma, double const C10) {
  double const s = sigma / distance;
  double s10 = s * s * s * s * s;
  s10 *= s10;
  energy += -4.0 * epsilon * C10 * s10;
  force -= 40.0 * epsilon * C10 * s10 / distance;
}

/** Interactions of a molecule with an image.
Coulomb interaction between @a w1 and the image of @a w2. Includes cutoff.
@param[in, out] w1 Molecule 1.
@param[in, out] w2 Molecule 2.
@param[in, out] U Add the potential energy to @a U.
@note Does not include the interaction of @a w2 with the image of @a w1. For
this you must explicitely call interactWithImage(w2, w1, U).
*/
template <class P>
void ZhuPhilpott<P>::interactWithImage(Water &w1, Water &w2, double &U) {
  // f is used to store the force on the atoms' images. These forces have no use
  // and are at the end discarded. we use the same vector to store of the forces
  // even when they are on different atoms' images.
  double r[4][3], f[3][3] = {{0}};
  for (int k = 0; k < 2; k++) {
    r[0][k] = w2.rh1_[k];
    r[1][k] = w2.rh2_[k];
    r[2][k] = w2.ro_[k];
    r[3][k] = w2.rc_[k];
  };
  // mirror about plane at z=0.
  r[0][2] = -w2.rh1_[2];
  r[1][2] = -w2.rh2_[2];
  r[2][2] = -w2.ro_[2];
  r[3][2] = -w2.rc_[2];
  Water w2m(r[0], r[1], r[2], r[3], f[0], f[1], f[2]);
  // -2 is given as the relative permittivity. This value is not the
  // permittivity of the metal which is infinite in the case of a perfect metal.
  // It is the value (1+epsilon_r)/(1-epsilon_r) which is equal to -1 for a
  // perfect metal. This value is the 'effective' permittivity for the
  // interaction of a charge with a mirror image.
  coulombWithCutoff(w1, w2m, U, -2.0);
  for (int k = 0; k < 2; k++) {
    w2.fh1_[k] += f[0][k];
    w2.fh2_[k] += f[1][k];
    w2.fo_[k] += f[2][k];
  };
  // mirror about plane
  w2.fh1_[2] -= f[0][2];
  w2.fh2_[2] -= f[1][2];
  w2.fo_[2] -= f[2][2];
}

/** Coulomb interaction between two molecules within the swithcing zone.
The precondition for using the function is that the distance between @a w1 and
@a w2 defined as @f$ |\mathbf r_{n1}-\mathbf r_{n2}| @f$ verifies the equation:
@f[
    C_{utoff}-S_{witchingWidth} < |\mathbf r_{n1}-\mathbf r_{n2}| < C_{utoff}
    @f]
@param[in, out] w1 Molecule 1.
@param[in, out] w2 Molecule 2.
@param[in, out] U Add the potential energy to @a U.
@param[in] relativePermittivity Relative permittivity.
@see getCutOff() and getSwitchingWidth().
*/
template <class P>
void ZhuPhilpott<P>::coulombWithCutoff(Water &w1, Water &w2, double &U,
                                       double const relativePermittivity) {
  double z[3], z1, z2;
  distance(w1.rc_, w2.rc_, z, z1, z2);
  if (z1 <= cutoff_ - switchingWidth_) {
    coulombFull(w1, w2, U, relativePermittivity);
  } else if (z1 < cutoff_) {
    double f1[3][3] = {{0}}, f2[3][3] = {{0}};
    Water v1(w1, f1[0], f1[1], f1[2]);
    Water v2(w2, f2[0], f2[1], f2[2]);
    double energy = 0.0;
    coulombFull(v1, v2, energy, relativePermittivity);
    ChargeGroup<3> g1 = {v1.rc_, 0, f1};
    ChargeGroup<3> g2 = {v2.rc_, 0, f2};
    switching(g1, g2, energy, cutoff_, switchingWidth_);
    U += energy;
    w1.addForces(v1);
    w2.addForces(v2);
  };
}

/** Coulomb interaction between two molecules with cutoff.
Coulomb interaction between two molecules of water.
@param[in, out] "w1, w2" Water molecules.
@param[in, out] U Add the potential energy to @a U.
@param[in] relativePermittivity
@see getCutOff() and getSwitchingWidth().
*/
template <class P>
void ZhuPhilpott<P>::coulombFull(Water &w1, Water &w2, double &U,
                                 double const relativePermittivity) {
  // q^2/relativePermittivity;
  const double qq2overEr = charge2_ / relativePermittivity;
  // Coulomb interactions between hydrogens
  coulomb(w1.rh1_, w2.rh1_, w1.fh1_, w2.fh1_, U, qq2overEr);
  coulomb(w1.rh1_, w2.rh2_, w1.fh1_, w2.fh2_, U, qq2overEr);
  coulomb(w1.rh2_, w2.rh1_, w1.fh2_, w2.fh1_, U, qq2overEr);
  coulomb(w1.rh2_, w2.rh2_, w1.fh2_, w2.fh2_, U, qq2overEr);
  // interactions between H and O.
  coulomb(w1.ro_, w2.rh1_, w1.fo_, w2.fh1_, U, -2.0 * qq2overEr);
  coulomb(w1.ro_, w2.rh2_, w1.fo_, w2.fh2_, U, -2.0 * qq2overEr);
  coulomb(w1.rh1_, w2.ro_, w1.fh1_, w2.fo_, U, -2.0 * qq2overEr);
  coulomb(w1.rh2_, w2.ro_, w1.fh2_, w2.fo_, U, -2.0 * qq2overEr);
  // interactions between O1, O2
  coulomb(w1.ro_, w2.ro_, w1.fo_, w2.fo_, U, 4.0 * qq2overEr);
}

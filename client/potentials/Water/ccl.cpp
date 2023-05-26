/** @file
Potential CCL Table II for intramolecular interactions in water.
@author Jean-Claude C. Berthet
@date 2008
University of Iceland
*/

#include "ccl.hpp"
#include <cassert>
#include <cmath>
// #include "unit_system.hpp"

using namespace std;
using namespace forcefields;
namespace {
#if defined(FORCEFIELDS_UNIT_SYSTEM_HPP) &&                                    \
    (FORCEFIELDS_UNIT_SYSTEM_HPP !=                                            \
     FORCEFIELDS_UNIT_SYSTEM_ELECTRONVOLT_ANGSTROM_FEMTOSECOND_ECHARGE)
using namespace unit_system;

double const re_ = 0.9572 * ANGSTROM;
double const thetae_ = 104.52 * DEGREE;

double const re2_ = re_ * re_;
//    double const ERGS_PER_ANGSTROM2=ERGS/ANGSTROM2;

// ------------------------ Quadratic ---------------------------
double const ro_2_ = 84.54e-12 * ERGS_PER_ANGSTROM2;
double const ro1_ro2_ = -1.01e-12 * ERGS_PER_ANGSTROM2;
double const ro_theta_ = 2.288e-12 * ERGS_PER_ANGSTROM2;
double const theta_2_ = 7.607e-12 * ERGS_PER_ANGSTROM2;

// --------------------------------------------- Cubic
// ----------------------------------------
double const ro_3_ = -10.168e-12 * ERGS_PER_ANGSTROM2;
double const ro_ro1_ro2_ = 0.201e-12 * ERGS_PER_ANGSTROM2;
double const ro_2_theta_ = 4.308e-12 * ERGS_PER_ANGSTROM2;
double const ro1_ro2_theta_ = -4.020e-12 * ERGS_PER_ANGSTROM2;
double const ro_theta_2_ = -1.175e-12 * ERGS_PER_ANGSTROM2;
double const thetat_3_ = -1.595e-12 * ERGS_PER_ANGSTROM2;

// --------------------------------------------- Quartic
// ----------------------------------------
double const ro_4_ = -10.684e-12 * ERGS_PER_ANGSTROM2;
double const ro1_ro2_ro_2_ = -6.162e-12 * ERGS_PER_ANGSTROM2;
double const ro1_2_ro2_2_ = 2.717e-12 * ERGS_PER_ANGSTROM2;

double const ro_3_theta_ = 6.328e-12 * ERGS_PER_ANGSTROM2;
double const ro_ro1_ro2_theta_ = -4.020e-12 * ERGS_PER_ANGSTROM2;

double const ro_2_theta_2_ = -4.70e-12 * ERGS_PER_ANGSTROM2;
double const ro1_ro2_theta_2_ = 3.05e-12 * ERGS_PER_ANGSTROM2;

// ro_theta_3 = 0
double const theta_4_ = -0.0318e-12 * ERGS_PER_ANGSTROM2;
#else
double const re_ = 0.9572;                  // ANGSTROM
double const thetae_ = 1.82421813418447321; // RADIANS

double const re2_ = 0.91623184; // ANGSTROM^2
//    double const ERGS_PER_ANGSTROM2 = 624150947960.771851; // eV / Angstrom^2

// ------------------------ Quadratic ---------------------------
double const ro_2_ = 52.7657211406036524;      // eV / Angstrom^2
double const ro1_ro2_ = -0.630392457440379528; // eV / Angstrom^2
double const ro_theta_ = 1.42805736893424595;  // eV / Angstrom^2
double const theta_2_ = 4.74791626113759158;   // eV / Angstrom^2

// --------------------------------------------- Cubic
// ----------------------------------------
double const ro_3_ = -6.3463668388651282;         // eV / Angstrom^2
double const ro_ro1_ro2_ = 0.125454340540115145;  // eV / Angstrom^2
double const ro_2_theta_ = 2.68884228381500501;   // eV / Angstrom^2
double const ro1_ro2_theta_ = -2.509086810802303; // eV / Angstrom^2
double const ro_theta_2_ = -0.733377363853906838; // eV / Angstrom^2
double const thetat_3_ = -0.995520761997431114;   // eV / Angstrom^2

// --------------------------------------------- Quartic
// ----------------------------------------
double const ro_4_ = -6.668428728012886;             // eV / Angstrom^2
double const ro1_ro2_ro_2_ = -3.84601814133427622;   // eV / Angstrom^2
double const ro1_2_ro2_2_ = 1.69581812560941692;     // eV / Angstrom^2
double const ro_3_theta_ = 3.94962719869576429;      // eV / Angstrom^2
double const ro_ro1_ro2_theta_ = -2.509086810802303; // eV / Angstrom^2
double const ro_2_theta_2_ = -2.93350945541562735;   // eV / Angstrom^2
double const ro1_ro2_theta_2_ = 1.90366039128035425; // eV / Angstrom^2
double const theta_4_ = -0.0198480001451525438;      // eV / Angstrom^2
#endif
} // namespace

/** @class forcefields::Ccl
@brief Quartic intramolecular potential for water.
This is a potential to intramolecular interaction in water published by @ref
carney1976 "Carney, Curtiss and Langhoff". The functions and parameters for the
potential are in Table II of the publication.
@section references References
@anchor carney1976
Improved Potential Functions for Bent AB2 Molecules: Water and Ozone. GD Carney,
LA Curtiss, SR Langhoff, J. Mol. Spectroscopy @b 1976, vol. 61, p. 371-381
*/

/// Distance OH at equilibrium
double const Ccl::re_ = ::re_;
/// Distance OH at equilibrium
double const Ccl::thetae_ = ::thetae_;

/** Compute the forces and the energy.
The order of the atoms is very important. For this function the order is H1, H1,
H2, H2, etc ... , O1, O2, etc ... The numbers are for the molecules. The suffix
HH_O_ is to reminds the ordering of the atoms. It indicates the position of the
atoms of the first molecule, while the underscore marks the locations of other
atoms.
@param[in]  nAtoms Number of Atoms.
@param[in]  R      coordinates.
@param[out] F     Forces.
@param[out] U     Potential energy.
@param[in]  b      Periodic boundaries.
@warning Be careful with the order of the atoms.
*/
void Ccl::computeHH_O_(const int nAtoms, const double R[], double F[],
                       double &U, const double b[]) {
  computeHH_O_(nAtoms, R, F, U, b, 0);
}

/** Compute the forces and the energy (with fixed atom optimisation).
In a simulation, some atoms may be fixed (not allowed to move). Computing the
interaction between two fixed atoms is useless. When the potential knows which
atoms are fixed, it can avoid these useless computations.
@param[in]  nAtoms Number of Atoms.
@param[in]  R      coordinates.
@param[out] F     Forces.
@param[out] U     Potential energy.
@param[in]  b      Periodic boundaries.
@param[in] fixed The length of the array must be equal to @a nAtoms in
compute(). True when atom is fixed, false otherwise.
@see Check compute() for the order of the atoms.
*/
void Ccl::computeHH_O_(const int nAtoms, const double R[], double F[],
                       double &U, const double b[], const bool fixed[]) {
  int const nMolecules = nAtoms / 3;
  const double(*const rh1)[6] = reinterpret_cast<const double(*)[6]>(R);
  const double(*const rh2)[6] = reinterpret_cast<const double(*)[6]>(&R[3]);
  const double(*const ro)[3] =
      reinterpret_cast<const double(*)[3]>(&R[nMolecules * 6]);
  double(*const fh1)[6] = reinterpret_cast<double(*)[6]>(F);
  double(*const fh2)[6] = reinterpret_cast<double(*)[6]>(&F[3]);
  double(*const fo)[3] = reinterpret_cast<double(*)[3]>(&F[nMolecules * 6]);
  bool const(*const xh1)[2] = reinterpret_cast<bool const(*)[2]>(fixed);
  bool const(*const xh2)[2] = reinterpret_cast<bool const(*)[2]>(&fixed[1]);
  bool const(*const xo)[1] =
      reinterpret_cast<bool const(*)[1]>(&fixed[nMolecules * 2]);
  computeTemplate(nMolecules, rh1, rh2, ro, fh1, fh2, fo, U, b, xh1, xh2, xo);
}

/// Name of the potential.
char const *Ccl::getName() const { return "Ccl"; }

/** Constructor with cutoff for derived classes.
Ccl does not use any cutoff (bond interaction only). The constructor is for
derived class that may require it.
*/
Ccl::Ccl(double cutoff, double switchingWidth)
    : PotentialBase(cutoff, switchingWidth) {}

/** @var forcefields::Ccl::Rho::n
@brief Convert derivative to force.
Vector @a n is by definition the vector which satisfy the equation:
@f[   \mathbf{F}=-\frac{dE}{d\rho}\mathbf{n}    @f]
It is equal to:
@f[\mathbf n=\frac{\mathbf{r}}{|\mathbf{r}|}\frac{r_e}{\rho^2} @f]
*/

/** @var forcefields::Ccl::Dtheta::n1
@brief Convert derivative to force.
The vector has should have direction of the force applied on hydrogen one. The
norm of the vector should be @f$ \frac{1}{r_1} @f$. This way multiplying the
vector by the moment gives the force on hydrogen 1.\n Reminder:\n The moment is
the equivalent of the force for angular coordinates:
@f[  M=-\frac{dE}{d\theta}   @f]
The moment is related to the force:
@f[   \mathbf M=\mathbf r \times \mathbf F      @f]
The vector @a n is defined as:
@f[   \mathbf n=\frac{\mathbf F}{M}       @f]
and
@f[     |\mathbf n|=\frac{1}{r}    @f]
*/

/// Initialise Rho.
void Ccl::initialiseRho(Vector3 const &v, Rho &ro) {
  ro._1 = (v._1 - re_) / v._1;
  ro._2 = ro._1 * ro._1;
  ro._3 = ro._2 * ro._1;
  double const d = re_ / v._2 / v._1;
  for (int i = 0; i < 3; ++i)
    ro.n[i] = v.v[i] * d;
}

/// Initialise Dtheta.
void Ccl::initialiseDtheta(Vector3 const &v1, Vector3 const &v2, Dtheta &dth,
                           double const thetaEquilibrium) {
  double const cos_theta = dotProduct(v1.v, v2.v) / v1._1 / v2._1;
  double const theta = acos(cos_theta);
  dth._1 = theta - thetaEquilibrium;
  dth._2 = dth._1 * dth._1;
  dth._3 = dth._2 * dth._1;
  double const d_theta = -1.0 / sqrt(1.0 - cos_theta * cos_theta);
  for (int k = 0; k < 3; ++k) {
    dth.n1[k] = 0.0;
    dth.n2[k] = 0.0;
    for (int j = 0; j < 3; ++j) {
      dth.n1[k] -= v2.v[j] / v2._1 / v1._1 * v1.v[j] * v1.v[k] / v1._2;
      dth.n2[k] -= v1.v[j] / v1._1 / v2._1 * v2.v[j] * v2.v[k] / v2._2;
    };
    dth.n1[k] += v2.v[k] / v2._1 / v1._1;
    dth.n2[k] += v1.v[k] / v1._1 / v2._1;
    dth.n1[k] *= d_theta;
    dth.n2[k] *= d_theta;
  };
}

/** Interactions inside one molecules.
@param[in] "rh1, rh2, ro"   Positions of H1, H2, O.
@param[in,out] "fh1, fh2, fo"   Forces on atoms.
@param[in, out] energy
@note The function is for one molecules so rh1[3], rh2[3], etc... @a energy, @a
fh1, etc, are incremented.
*/
void Ccl::intramolecular(double const rh1[], double const rh2[],
                         double const ro[], double fh1[], double fh2[],
                         double fo[], double &energy) {
  Vector3 v1, v2;
  distance(rh1, ro, v1);
  distance(rh2, ro, v2);
  // ------prepare _
  Rho ro1, ro2;
  initialiseRho(v1, ro1);
  initialiseRho(v2, ro2);

  // ------- prepare Delta theta
  Dtheta dth;
  initialiseDtheta(v1, v2, dth, thetae_);
  intramolecular(ro1, ro2, dth, energy, fh1, fh2, fo);
}

/** Compute intramolecular energy and forces.
@param[in]  "ro1, ro2, dth"      Potential coordinates properly initialised.
@param[in,out]    energy      Add the energy.
@param[in,out]    "fh1, fh2, fo"    Add the forces on atoms.
*/
void Ccl::intramolecular(Rho const &ro1, Rho const &ro2, Dtheta const &dth,
                         double &energy, double fh1[], double fh2[],
                         double fo[]) {
  double d1 = 0, d2 = 0,
         d3 = 0; // d1=d(E)/d(rho1), d2=d(E)/d(rho2), d(E)/d(theta)

  // -------- quadratic interaction
  ro_2(ro1, ro2, energy, d1, d2);
  ro1_ro2(ro1, ro2, energy, d1, d2);
  ro_theta(ro1, ro2, dth, energy, d1, d2, d3);
  theta_2(dth, energy, d3);
  //*/

  // -------- Cubic interaction
  ro_3(ro1, ro2, energy, d1, d2);
  ro_ro1_ro2(ro1, ro2, energy, d1, d2);
  ro_2_theta(ro1, ro2, dth, energy, d1, d2, d3);
  ro1_ro2_theta(ro1, ro2, dth, energy, d1, d2, d3);
  ro_theta_2(ro1, ro2, dth, energy, d1, d2, d3);
  theta_3(dth, energy, d3);
  //*/

  // --------------------------------------------- Quartic
  // -------------------------------------
  ro_4(ro1, ro2, energy, d1, d2);
  ro1_ro2_ro_2(ro1, ro2, energy, d1, d2);
  ro1_2_ro2_2(ro1, ro2, energy, d1, d2);

  ro_3_theta(ro1, ro2, dth, energy, d1, d2, d3);
  ro_ro1_ro2_theta(ro1, ro2, dth, energy, d1, d2, d3);

  ro_2_theta_2(ro1, ro2, dth, energy, d1, d2, d3);
  ro1_ro2_theta_2(ro1, ro2, dth, energy, d1, d2, d3);

  // ro_theta_3 = 0
  theta_4(dth, energy, d3);
  //*/
  for (int i = 0; i < 3; ++i) {
    fh1[i] -= d1 * ro1.n[i] + d3 * dth.n1[i];
    fh2[i] -= d2 * ro2.n[i] + d3 * dth.n2[i];
    fo[i] += d1 * ro1.n[i] + d2 * ro2.n[i] + d3 * dth.n1[i] + d3 * dth.n2[i];
  };
}

/** Compute forces and energy.
The function calculates and returns the energy and forces applied on each atom.
The function has preconditions with which the user must comply before calling
this function. The periodic boundaries must set with (setPeriodicity()). The
content of and arrays @a fh1 , @a fh2 , @a fo must be zero for all elements. \n
Arrays ending in @em h1 are for the hydrogen one of the molecules and those
ending in @em h2 for hydrogen 2.
@param[in] nMolecules Number of molecules.
@param[in] "rh1, rh2, ro" Positions of hydrogens and oxygens.
@param[in,out] "fh1, fh2, fo" Forces on hydrogens and oxygens.
@param[in,out] energy      Potential energy.
@param[in] b      Periodic boundaries.
@param[in] "xh1, xh2, xo" Tell which atoms are fixed and which are movable. This
parameter are optional. When provided, interaction between fixed atoms are
skipped. The three arrays @a xh1, @a xh2, @a xo must be provided for the
optimisation to work.
@warning Remember the preconditions.
*/
template <int H, int O>
void Ccl::computeTemplate(
    const int nMolecules, const double (*const rh1)[H * 3],
    const double (*const rh2)[H * 3], const double (*const ro)[O * 3],
    double (*const fh1)[H * 3], double (*const fh2)[H * 3],
    double (*const fo)[O * 3], double &energy, double const b[],
    bool const (*const xh1)[H], bool const (*const xh2)[H],
    bool const (*const xo)[O]) {
  for (int i = 0; i < nMolecules; ++i) {
    for (int a = 0; a < 3; a++) {
      fh1[i][a] = 0.0;
      fh2[i][a] = 0.0;
      fo[i][a] = 0.0;
    };
  };
  energy = 0.0;
  setPeriodicity(b);

  for (int i = nMolecules - 1; i >= 0; --i) {
    bool const areFixed = xh1[i][0] and xh2[i][0] and xo[i][0];
    if (not areFixed) {
      intramolecular(rh1[i], rh2[i], ro[i], fh1[i], fh2[i], fo[i], energy);
    };
  };
  assert(not isnan(energy) and not isinf(energy));
}

// --------------------------------------------- Square
// ----------------------------------------

/** Energy and derivatives for Term 1.
Compute energy for term : @f$ r_e^2(\rho_1^2+\rho_2^2)/2 @f$ and derivatives:
@f$ \frac{dE}{d\rho_1} @f$ and @f$ \frac{dE}{d\rho_2} @f$. We remind that @f$
\rho_1=\frac{r_1-r_{eq}}{\rho_1} @f$.
@param[in] "ro1, ro2"      Two sets containing parameters rho and rho^2.
@param[in,out] energy          Energy. Add energy to @a e.
@param[in,out] "d1, d2"   Derivatives. Add derivatives to parameters @a d1 and
@a d2.
*/
void Ccl::ro_2(Rho const &ro1, Rho const &ro2, double &energy, double &d1,
               double &d2) {
  energy += ro_2_ * re2_ * (ro1._2 + ro2._2) / 2.0;
  d1 += ro_2_ * re2_ * ro1._1;
  d2 += ro_2_ * re2_ * ro2._1;
}

void Ccl::ro1_ro2(Rho const &ro1, Rho const &ro2, double &energy, double &d1,
                  double &d2) {
  energy += ro1_ro2_ * re2_ * ro1._1 * ro2._1;
  d1 += ro1_ro2_ * re2_ * ro2._1;
  d2 += ro1_ro2_ * re2_ * ro1._1;
}

void Ccl::ro_theta(Rho const &ro1, Rho const &ro2, Dtheta const &dth,
                   double &energy, double &d1, double &d2, double &d3) {
  energy += ro_theta_ * re2_ * (ro1._1 + ro2._1) * dth._1;
  d1 += ro_theta_ * re2_ * dth._1;            // dE/drho1
  d2 += ro_theta_ * re2_ * dth._1;            // dE/drho2
  d3 += ro_theta_ * re2_ * (ro1._1 + ro2._1); // dE/d(theta)
}

/** Energy and derivatives for Term 4.
Compute energy for term : @f$ (r_e\Delta \theta)^2 @f$ and derivatives: @f$
\frac{dE}{d\theta} @f$.
@param[in] dth    Two sets containing parameters @f$ (\Delta \theta)^2@f$ and
@f$ \Delta \theta @f$.
@param[in,out] energy          Energy. Add energy to @a e.
@param[in,out] d3   Derivatives. Add derivatives to parameters @a d3.
@note Parameter @a d3 is a moment of force.
*/
void Ccl::theta_2(Dtheta const &dth, double &energy, double &d3) {
  energy += theta_2_ * re2_ * dth._2 / 2.0;
  d3 += theta_2_ * re2_ * dth._1; // dE/d(theta)= - moment
}

// --------------------------------------------- Cubic
// ----------------------------------------

void Ccl::ro_3(Rho const &ro1, Rho const &ro2, double &energy, double &d1,
               double &d2) {
  energy += ro_3_ * re2_ * (ro1._3 + ro2._3);
  d1 += 3.0 * ro_3_ * re2_ * ro1._2;
  d2 += 3.0 * ro_3_ * re2_ * ro2._2;
}

void Ccl::ro_ro1_ro2(Rho const &ro1, Rho const &ro2, double &energy, double &d1,
                     double &d2) {
  energy += ro_ro1_ro2_ * re2_ * (ro1._1 + ro2._1) * ro1._1 * ro2._1;
  d1 += ro_ro1_ro2_ * re2_ * (2.0 * ro1._1 * ro2._1 + ro2._2);
  d2 += ro_ro1_ro2_ * re2_ * (2.0 * ro1._1 * ro2._1 + ro1._2);
}

void Ccl::ro_2_theta(Rho const &ro1, Rho const &ro2, Dtheta const &dth,
                     double &energy, double &d1, double &d2, double &d3) {
  energy += ro_2_theta_ * re2_ * (ro1._2 + ro2._2) * dth._1;
  d1 += ro_2_theta_ * re2_ * 2.0 * ro1._1 * dth._1; // dE/drho1
  d2 += ro_2_theta_ * re2_ * 2.0 * ro2._1 * dth._1; // dE/drho2
  d3 += ro_2_theta_ * re2_ * (ro1._2 + ro2._2);     // dE/d(theta)
}

void Ccl::ro1_ro2_theta(Rho const &ro1, Rho const &ro2, Dtheta const &dth,
                        double &energy, double &d1, double &d2, double &d3) {
  energy += ro1_ro2_theta_ * re2_ * ro1._1 * ro2._1 * dth._1;
  d1 += ro1_ro2_theta_ * re2_ * ro2._1 * dth._1; // dE/drho1
  d2 += ro1_ro2_theta_ * re2_ * ro1._1 * dth._1; // dE/drho2
  d3 += ro1_ro2_theta_ * re2_ * ro1._1 * ro2._1; // dE/d(theta)
}

void Ccl::ro_theta_2(Rho const &ro1, Rho const &ro2, Dtheta const &dth,
                     double &energy, double &d1, double &d2, double &d3) {
  energy += ro_theta_2_ * re2_ * (ro1._1 + ro2._1) * dth._2;
  d1 += ro_theta_2_ * re2_ * dth._2;                           // dE/drho1
  d2 += ro_theta_2_ * re2_ * dth._2;                           // dE/drho2
  d3 += ro_theta_2_ * re2_ * (ro1._1 + ro2._1) * 2.0 * dth._1; // dE/d(theta)
}

void Ccl::theta_3(Dtheta const &dth, double &energy, double &d3) {
  energy += thetat_3_ * re2_ * dth._3;
  d3 += thetat_3_ * re2_ * 3.0 * dth._2; // dE/d(theta)
}

// --------------------------------------------- Quartic
// ----------------------------------------

void Ccl::ro_4(Rho const &ro1, Rho const &ro2, double &energy, double &d1,
               double &d2) {
  energy += ro_4_ * re2_ * (ro1._2 * ro1._2 + ro2._2 * ro2._2);
  d1 += ro_4_ * re2_ * 4.0 * ro1._3; // dE/drho1
  d2 += ro_4_ * re2_ * 4.0 * ro2._3; // dE/drho2
}

void Ccl::ro1_ro2_ro_2(Rho const &ro1, Rho const &ro2, double &energy,
                       double &d1, double &d2) {
  energy += ro1_ro2_ro_2_ * re2_ * ro1._1 * ro2._1 * (ro1._2 + ro2._2);
  d1 += ro1_ro2_ro_2_ * re2_ * (3.0 * ro2._1 * ro1._2 + ro2._3); // dE/drho1
  d2 += ro1_ro2_ro_2_ * re2_ * (3.0 * ro1._1 * ro2._2 + ro1._3); // dE/drho2
}

void Ccl::ro1_2_ro2_2(Rho const &ro1, Rho const &ro2, double &energy,
                      double &d1, double &d2) {
  energy += ro1_2_ro2_2_ * re2_ * ro1._2 * ro2._2;
  d1 += ro1_2_ro2_2_ * re2_ * 2.0 * ro1._1 * ro2._2; // dE/drho1
  d2 += ro1_2_ro2_2_ * re2_ * 2.0 * ro1._2 * ro2._1; // dE/drho2
}

void Ccl::ro_3_theta(Rho const &ro1, Rho const &ro2, Dtheta const &dth,
                     double &energy, double &d1, double &d2, double &d3) {
  energy += ro_3_theta_ * re2_ * (ro1._3 + ro2._3) * dth._1;
  d1 += ro_3_theta_ * re2_ * 3.0 * ro1._2 * dth._1; // dE/drho1
  d2 += ro_3_theta_ * re2_ * 3.0 * ro2._2 * dth._1; // dE/drho2
  d3 += ro_3_theta_ * re2_ * (ro1._3 + ro2._3);     // dE/d(theta)
}

void Ccl::ro_ro1_ro2_theta(Rho const &ro1, Rho const &ro2, Dtheta const &dth,
                           double &energy, double &d1, double &d2, double &d3) {
  energy +=
      ro_ro1_ro2_theta_ * re2_ * (ro1._1 + ro2._1) * ro1._1 * ro2._1 * dth._1;
  d1 += ro_ro1_ro2_theta_ * re2_ * (2.0 * ro1._1 * ro2._1 + ro2._2) *
        dth._1; // dE/drho1
  d2 += ro_ro1_ro2_theta_ * re2_ * (2.0 * ro1._1 * ro2._1 + ro1._2) *
        dth._1; // dE/drho2
  d3 += ro_ro1_ro2_theta_ * re2_ * (ro1._1 + ro2._1) * ro1._1 *
        ro2._1; // dE/d(theta)
}

void Ccl::ro_2_theta_2(Rho const &ro1, Rho const &ro2, Dtheta const &dth,
                       double &energy, double &d1, double &d2, double &d3) {
  energy += ro_2_theta_2_ * re2_ * (ro1._2 + ro2._2) * dth._2;
  d1 += ro_2_theta_2_ * re2_ * 2.0 * ro1._1 * dth._2;            // dE/drho1
  d2 += ro_2_theta_2_ * re2_ * 2.0 * ro2._1 * dth._2;            // dE/drho2
  d3 += ro_2_theta_2_ * re2_ * 2.0 * (ro1._2 + ro2._2) * dth._1; // dE/d(theta)
}

void Ccl::ro1_ro2_theta_2(Rho const &ro1, Rho const &ro2, Dtheta const &dth,
                          double &energy, double &d1, double &d2, double &d3) {
  energy += ro1_ro2_theta_2_ * re2_ * ro1._1 * ro2._1 * dth._2;
  d1 += ro1_ro2_theta_2_ * re2_ * ro2._1 * dth._2;                // dE/drho1
  d2 += ro1_ro2_theta_2_ * re2_ * ro1._1 * dth._2;                // dE/drho2
  d3 += ro1_ro2_theta_2_ * re2_ * ro1._1 * ro2._1 * 2.0 * dth._1; // dE/d(theta)
}

// rCcl::o_theta_3 = 0
void Ccl::theta_4(Dtheta const &dth, double &energy, double &d3) {
  energy += theta_4_ * re2_ * dth._2 * dth._2;
  d3 += theta_4_ * re2_ * 4.0 * dth._3; // dE/d(theta)
}

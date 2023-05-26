/** @file
TIP4P potential for water
@author Jean-Claude C. Berthet
@date 2006-2007
University of Iceland
*/
#include "tip4p_ccl.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
// If unit_system.hpp is removed the system of unit will be  (Energy: eV,
// distance: Angstrom, time: fs, charge: e).
#include "tip4p_unit_system.hpp"

using namespace std;
using namespace forcefields;
/* defined(FORCEFIELDS_UNIT_SYSTEM_HPP) && \
( FORCEFIELDS_UNIT_SYSTEM_HPP !=
FORCEFIELDS_UNIT_SYSTEM_ELECTRONVOLT_ANGSTROM_FEMTOSECOND_ECHARGE) //*/
using namespace forcefields::unit_system;
namespace {
double const re_ = 0.9572 * ANGSTROM;
double const thetae_ = 104.52 * DEGREE;
double const charge_ = 0.520 * ECHARGE;    ///< Charge on one hydrogen
double const charge2_ = charge_ * charge_; ///< Square of # charge_
double const sigma_ =
    3.154 * ANGSTROM; ///< Lennard-Jones sigma between oxygen atoms.
double const epsilon_ =
    78.0 * KELVIN; ///< Lennard-Jones epsilon between oxygen atoms.
double const ron_ =
    0.150 * ANGSTROM; ///< Distance between oxygen and the middle charge N.
double const rok_ =
    ::re_ * cos(::thetae_ / 2.0); ///< Distance between oxygen and the centre of
                                  ///< the two hydrogen atoms (point K).
double const wh_ = ron_ / rok_ * 0.5;
double const wo_ = (1.0 - wh_ * 2.0);
} // namespace

/** @class forcefields::Tip4p
@brief TIP4P forcefield with fexible molecules and group based cutoff.
The function computing the forces and energy is compute(). Parameters for the
forcefield were taken from @ref abascal2005 "Abascal et al.". The parameters are
those of TIP4P (original version). TIP4P is a potential for constrained
molecules of water. In this implementation quadratic restraints were added to
bonds OH and HH inside molecules of water. So the potential can be used without
constraints. The potential has a cutoff for long range interaction. The cutoff
can be changed by setCutoff(). There is also and switching zone at the edge of
the cutoff to cut off the interactions smoothly. The width of this switching
zone is controlled by setSwitchingWidth().\n The unit system for this potential
is eV (electron volt), Angstrom, e (e charge).
@section references References
@anchor abascal2005
A general purpose model for the condensed phases of water: TIP4P/2005, J.L.F.
Abascal and C. Vega, J. Chem. Phys. (2005) vol. 123, p. 234505.*/

Tip4p::Tip4p() : Ccl() {}

Tip4p::Tip4p(double cutoff, double switchingWidth)
    : Ccl(cutoff, switchingWidth) {}

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
@warning Be careful with the order of the atoms.*/
void Tip4p::computeHH_O_(const int nAtoms, const double R[], double F[],
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
void Tip4p::computeHH_O_(const int nAtoms, const double R[], double F[],
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
char const *Tip4p::getName() { return "Tip4p"; }

// -----------------------------------------------------private----------------------------------------

/// Pointers for one molecule of water.
struct Tip4p::Water {
  const double *const rh1_;
  const double *const rh2_;
  const double *const ro_;
  const double *const rn_;
  const double *const rc_;
  double *const fh1_;
  double *const fh2_;
  double *const fo_;
  double *const fn_;
};

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
void Tip4p::computeTemplate(
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
    double rc1[3], rn1[3], fn1[3] = {0}; // rn position charge N, rc centre of
                                         // charge for groupbased cutoff.
    calculateWeightedCentre(wo_, wh_, wh_, ro[i], rh1[i], rh2[i], rn1);
    calculateCentre(rn1, rh1[i], rh2[i], rc1);
    Water w1 = {rh1[i], rh2[i], ro[i], rn1, rc1, fh1[i], fh2[i], fo[i], fn1};
    intramolecular(rh1[i], rh2[i], ro[i], fh1[i], fh2[i], fo[i], energy);
    for (int j = i - 1; j >= 0; --j) {
      bool areFixed = false;
      if (xh1 and xh2 and xo) {
        areFixed = xh1[i][0] and xh2[i][0] and xo[i][0];
        // check if all the atoms of molecule j are fixed.
        areFixed &= xh1[j][0] and xh2[j][0] and xo[j][0];
      };
      // if both molecules are fixed skip force calculation
      if (not areFixed) {
        double rc2[3], rn2[3], fn2[3] = {0};
        calculateWeightedCentre(wo_, wh_, wh_, ro[j], rh1[j], rh2[j], rn2);
        calculateCentre(rn2, rh1[j], rh2[j], rc2);
        Water w2 = {rh1[j], rh2[j], ro[j], rn2, rc2,
                    fh1[j], fh2[j], fo[j], fn2};
        coulombWithCutoff(w1, w2, energy);
        spreadWeightedForce(wo_, wh_, wh_, w2.fo_, w2.fh1_, w2.fh2_, w2.fn_);
        lennardJonesWithCutoff(w1, w2, energy);
      };
    };
    spreadWeightedForce(wo_, wh_, wh_, w1.fo_, w1.fh1_, w1.fh2_, w1.fn_);
  };
  assert(not isnan(energy) and not isinf(energy));
}

/** Interactions between two molecules.
Coulomb interaction between two molecules of water with molecules-based cutoff
@param[in, out] w1 Molecule 1.
@param[in, out] w2 Molecule 2.
@param[in, out] U Add the potential energy to @a U.
*/
void Tip4p::coulombWithCutoff(Water &w1, Water &w2, double &U) {
  double z[3], z1, z2;
  distance(w1.rc_, w2.rc_, z, z1, z2);
  if (z1 <= cutoff_ - switchingWidth_) {
    coulombFull(w1, w2, U);
  } else if (z1 < cutoff_) {
    double f1[3][3] = {{0}}, f2[3][3] = {{0}};
    // f1[0], f1[1], 0, f1[2] are respectively, H1, H2, O, N. There is no charge
    // on O.
    Water v1 = {w1.rh1_, w1.rh2_, w1.ro_, w1.rn_, w1.rc_,
                f1[0],   f1[1],   0,      f1[2]};
    Water v2 = {w2.rh1_, w2.rh2_, w2.ro_, w2.rn_, w2.rc_,
                f2[0],   f2[1],   0,      f2[2]};
    double energy = 0.0;
    coulombFull(v1, v2, energy);
    ChargeGroup<3> g1 = {v1.rc_, 0, f1};
    ChargeGroup<3> g2 = {v2.rc_, 0, f2};
    // Calculate the weakened forces and energy.
    switching(g1, g2, energy, cutoff_, switchingWidth_);
    // add weakened force and energy to those of other interactions.
    U += energy;
    for (int i = 0; i < 3; ++i) {
      w1.fh1_[i] += v1.fh1_[i];
      w1.fh2_[i] += v1.fh2_[i];
      w1.fn_[i] += v1.fn_[i];
      w2.fh1_[i] += v2.fh1_[i];
      w2.fh2_[i] += v2.fh2_[i];
      w2.fn_[i] += v2.fn_[i];
    };
  }
}

/** Interactions between two molecules.
Coulomb interaction between two molecules. Full interaction (i.e. no cutoff).
@param[in, out] w1 Molecule 1.
@param[in, out] w2 Molecule 2.
@param[in, out] U Add the potential energy to @a U.
@see intermolecularFull() and intermolecularSwitching().
@note The function add the force on charge N to Water::fn_ but does not spread
it to Water::fh1_, etc ...
@see spreadN().
*/
void Tip4p::coulombFull(Water &w1, Water &w2, double &U) {
  // Coulomb interactions between hydrogens
  coulomb(w1.rh1_, w2.rh1_, w1.fh1_, w2.fh1_, U, charge2_);
  coulomb(w1.rh1_, w2.rh2_, w1.fh1_, w2.fh2_, U, charge2_);
  coulomb(w1.rh2_, w2.rh1_, w1.fh2_, w2.fh1_, U, charge2_);
  coulomb(w1.rh2_, w2.rh2_, w1.fh2_, w2.fh2_, U, charge2_);
  // interactions between H and N.
  coulomb(w1.rn_, w2.rh1_, w1.fn_, w2.fh1_, U, -2.0 * charge2_);
  coulomb(w1.rn_, w2.rh2_, w1.fn_, w2.fh2_, U, -2.0 * charge2_);
  coulomb(w1.rh1_, w2.rn_, w1.fh1_, w2.fn_, U, -2.0 * charge2_);
  coulomb(w1.rh2_, w2.rn_, w1.fh2_, w2.fn_, U, -2.0 * charge2_);
  // interactions between N1, N2
  coulomb(w1.rn_, w2.rn_, w1.fn_, w2.fn_, U, 4.0 * charge2_);
}

/** Interactions between two molecules.
Lennard Jones between oxygen only with cutoff
@param[in, out] w1 Molecule of water 1.
@param[in, out] w2 Molecule of water 2.
@param[in, out] U Add the potential energy to @a U.
@see intermolecularFull() and intermolecularSwitching().
*/
void Tip4p::lennardJonesWithCutoff(Water &w1, Water &w2, double &U) {
  double z[3], z1;
  distance(w1.ro_, w2.ro_, z, z1);
  if (z1 <= cutoff_ - switchingWidth_) {
    lennardJones(w1.ro_, w2.ro_, w1.fo_, w2.fo_, U, epsilon_, sigma_);
  } else if (z1 < cutoff_) {
    double f1[3] = {0}, f2[3] = {0};
    double energy = 0.0;
    lennardJones(w1.ro_, w2.ro_, f1, f2, energy, epsilon_, sigma_);
    switching(w1.ro_, w2.ro_, f1, f2, energy);
    U += energy;
    for (int i = 0; i < 3; i++) {
      w1.fo_[i] += f1[i];
      w2.fo_[i] += f2[i];
    };
  };
}

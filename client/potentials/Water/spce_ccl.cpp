/** @file
SPC/E potential for water.
@author Jean-Claude C. Berthet
@date 2006-2007
University of Iceland
*/
#include "spce_ccl.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
// #include "unit_system.hpp"

using namespace std;
using namespace forcefields;

#if defined(FORCEFIELDS_UNIT_SYSTEM_HPP) &&                                    \
    (FORCEFIELDS_UNIT_SYSTEM_HPP !=                                            \
     FORCEFIELDS_UNIT_SYSTEM_ELECTRONVOLT_ANGSTROM_FEMTOSECOND_ECHARGE)
using namespace unit_system;
const double SpceCcl::roh_ = 1.0 * ANGSTROM;
const double SpceCcl::theta_ =
    acos(-1.0 / 3.0); // tetrahedron about 109.47*DEGREE
const double SpceCcl::rhh_ = isoscelesBase(roh_, theta_);
const double SpceCcl::charge_ = 0.4238 * ECHARGE;
const double SpceCcl::charge2_ = charge_ * charge_;
// Definition of A and B are inverted compared to original publication by
// Berendsen
const double SpceCcl::A_ = pow(0.3428 * NM, 12.0) * KJ_PER_MOL;
const double SpceCcl::B_ = pow(0.37122 * NM, 6.0) * KJ_PER_MOL;
const double SpceCcl::sigma_ = sigma(A_, B_);
const double SpceCcl::epsilon_ = epsilon(A_, B_);
const double SpceCcl::polarisationEnergy_ = 5.22 * KJ_PER_MOL;
#else
const double SpceCcl::roh_ = 1.0;                      // Angstrom
const double SpceCcl::theta_ = 1.91063;                // radians
const double SpceCcl::rhh_ = 1.63299;                  // Angstrom
const double SpceCcl::charge_ = 0.4238;                // e
const double SpceCcl::charge2_ = 0.179606;             // e
const double SpceCcl::A_ = 27291.6;                    // Angstrom^-1
const double SpceCcl::B_ = 27.1223;                    // Angstrom^-1
const double SpceCcl::sigma_ = 3.16556;                // Angstrom
const double SpceCcl::epsilon_ = 0.00673853;           // eV
const double SpceCcl::polarisationEnergy_ = 0.0541015; // eV
#endif

SpceCcl::SpceCcl() : Ccl() {}

SpceCcl::SpceCcl(double cutoff, double switchingWidth)
    : Ccl(cutoff, switchingWidth) {}

void SpceCcl::computeHH_O_(const int nAtoms, const double R[], double F[],
                           double &U, const double b[]) {
  computeHH_O_(nAtoms, R, F, U, b, 0);
}

void SpceCcl::computeHH_O_(const int nAtoms, const double R[], double F[],
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

char const *SpceCcl::getName() const { return "SpceCcl"; }

void SpceCcl::lennardJonesWithCutoff(Water &w1, Water &w2, double &U) {
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

void SpceCcl::coulombWithCutoff(Water &w1, Water &w2, double &U) {
  double z[3], z1, z2;
  distance(w1.rc_, w2.rc_, z, z1, z2);
  if (z1 <= cutoff_ - switchingWidth_) {
    coulombFull(w1, w2, U);
  } else if (z1 < cutoff_) {
    // store forces and energy of the full interaction in temporaries
    double f1[3][3] = {{0}}, f2[3][3] = {{0}};
    Water v1(w1.rh1_, w1.rh2_, w1.ro_, w1.rc_, f1[0], f1[1], f1[2]);
    Water v2(w2.rh1_, w2.rh2_, w2.ro_, w2.rc_, f2[0], f2[1], f2[2]);
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
      w1.fo_[i] += v1.fo_[i];
      w2.fh1_[i] += v2.fh1_[i];
      w2.fh2_[i] += v2.fh2_[i];
      w2.fo_[i] += v2.fo_[i];
    };
  };
}

void SpceCcl::coulombFull(Water &w1, Water &w2, double &U) {
  // Coulomb interactions between hydrogens
  coulomb(w1.rh1_, w2.rh1_, w1.fh1_, w2.fh1_, U, charge2_);
  coulomb(w1.rh1_, w2.rh2_, w1.fh1_, w2.fh2_, U, charge2_);
  coulomb(w1.rh2_, w2.rh1_, w1.fh2_, w2.fh1_, U, charge2_);
  coulomb(w1.rh2_, w2.rh2_, w1.fh2_, w2.fh2_, U, charge2_);
  // interactions between H and O.
  coulomb(w1.ro_, w2.rh1_, w1.fo_, w2.fh1_, U, -2.0 * charge2_);
  coulomb(w1.ro_, w2.rh2_, w1.fo_, w2.fh2_, U, -2.0 * charge2_);
  coulomb(w1.rh1_, w2.ro_, w1.fh1_, w2.fo_, U, -2.0 * charge2_);
  coulomb(w1.rh2_, w2.ro_, w1.fh2_, w2.fo_, U, -2.0 * charge2_);
  // interactions between O1, O2
  coulomb(w1.ro_, w2.ro_, w1.fo_, w2.fo_, U, 4.0 * charge2_);
}

void SpceCcl::initialiseRho(Vector3 const &v, Rho &ro) {
  // re_ = r_CCL and roh_ = r_SPC
  ro._1 = (v._1 - roh_) / (v._1 - roh_ + re_);
  ro._2 = ro._1 * ro._1;
  ro._3 = ro._2 * ro._1;
  double const a = v._1 - roh_ + re_;
  double const d = re_ / a / a / v._1;
  for (int i = 0; i < 3; ++i)
    ro.n[i] = v.v[i] * d;
}

void SpceCcl::intramolecular(Water &w, double &energy) {
  Vector3 v1, v2;
  distance(w.rh1_, w.ro_, v1);
  distance(w.rh2_, w.ro_, v2);
  // ------ prepare --------------------
  Rho ro1, ro2;
  initialiseRho(v1, ro1);
  initialiseRho(v2, ro2);

  // ------- prepare Delta theta
  Dtheta dth;
  initialiseDtheta(v1, v2, dth, theta_); // theta_ SPC equilibrium angle.
  Ccl::intramolecular(ro1, ro2, dth, energy, w.fh1_, w.fh2_, w.fo_);
}

template <int H, int O>
void SpceCcl::computeTemplate(
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
    double rc1[3] = {0};
    calculateCentre(rh1[i], rh2[i], ro[i], rc1);
    Water w1(rh1[i], rh2[i], ro[i], rc1, fh1[i], fh2[i], fo[i]);
    intramolecular(w1, energy);
    energy += polarisationEnergy_;
    for (int j = i - 1; j >= 0; --j) {
      bool areFixed = false;
      if (xh1 and xh2 and xo) {
        areFixed = xh1[i][0] and xh2[i][0] and xo[i][0];
        // check if all the atoms of molecule j are fixed.
        areFixed &= xh1[j][0] and xh2[j][0] and xo[j][0];
      };
      // if both molecules are fixed skip force calculation
      if (not areFixed) {
        double rc2[3] = {0};
        calculateCentre(rh1[j], rh2[j], ro[j], rc2);
        Water w2(rh1[j], rh2[j], ro[j], rc2, fh1[j], fh2[j], fo[j]);
        lennardJonesWithCutoff(w1, w2, energy);
        coulombWithCutoff(w1, w2, energy);
      };
    };
  };
  assert(not isnan(energy) and not isinf(energy));
}

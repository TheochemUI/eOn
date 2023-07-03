#pragma once
#ifndef FORCEFIELDS_SPCE_CCL_HPP
#define FORCEFIELDS_SPCE_CCL_HPP
#include "ccl.hpp"
/** @file
SPC/E potential for water.
@author Jean-Claude C. Berthet
@date 2006-2007
University of Iceland
*/

namespace forcefields {
/** SPC/E+CCl water potential.
This potential combines potential SPC/E (for intermolecular interactions) with
potential CCL (for intramolecular interactions). The functions and parameters
for the SPC/E are from @ref berendsen1987 "Berendsen, et al.". References for
potential CCL can be found in in class Ccl 's documentation. The rules used to
merge these two potentials were taken from @ref amira2004 "Amira, et al".
@section references References
@anchor berendsen1987
The Missing Term in Effective Pair Potentials. HJC Berendsen, JR Grigera et al.
J. Phys. Chem. @b 1987, vol. 91, p. 6269-6271 \n
@anchor amira2004
Derivation and evaluation of a flexible SPC model for liquid water. S Amira, D
Spangberg, K Hermansson, Chemical Physics, @b 2004, vol. 303, p. 327-334
*/
class SpceCcl : public Ccl {
public:
  SpceCcl();
  SpceCcl(double cutoff, double switchingWidth);

  /** Compute the forces and the energy.
      The order of the atoms is very important. For this function the order is
     H1, H1, H2, H2, etc ... , O1, O2, etc ... The numbers are for the
     molecules. The suffix HH_O_ is to reminds the ordering of the atoms. It
     indicates the position of the atoms of the first molecule, while the
     underscore marks the locations of other atoms.
      @param[in]  nAtoms Number of Atoms.
      @param[in]  R      coordinates.
      @param[out] F     Forces.
      @param[out] U     Potential energy.
      @param[in]  b      Periodic boundaries.
      @warning Be careful with the order of the atoms.*/
  void computeHH_O_(const int nAtoms, const double R[], double F[], double &U,
                    const double b[]);

  /** Compute the forces and the energy (with fixed atom optimisation).
      In a simulation, some atoms may be fixed (not allowed to move). Computing
     the interaction between two fixed atoms is useless. When the potential
     knows which atoms are fixed, it can avoid these useless computations.
      @param[in]  nAtoms Number of Atoms.
      @param[in]  R      coordinates.
      @param[out] F     Forces.
      @param[out] U     Potential energy.
      @param[in]  b      Periodic boundaries.
      @param[in] fixed The length of the array must be equal to @a nAtoms in
     compute(). True when atom is fixed, false otherwise.
      @see Check compute() for the order of the atoms.
      */
  void computeHH_O_(const int nAtoms, const double R[], double F[], double &U,
                    const double b[], const bool fixed[]);

  char const *getName() const; ///< Name of the potential.

protected:
  /** Pointers to molecule of water.
  Pointers to coordinates and forces of a molecule of water.
  */
  struct Water {
    Water(double const rh1[], double const rh2[], double const ro[],
          double const rc[], double fh1[], double fh2[], double fo[])
        : rh1_(rh1), rh2_(rh2), ro_(ro), rc_(rc), fh1_(fh1), fh2_(fh2),
          fo_(fo) {}
    Water(Water const &w, double fh1[], double fh2[], double fo[])
        : rh1_(w.rh1_), rh2_(w.rh2_), ro_(w.ro_), rc_(w.rc_), fh1_(fh1),
          fh2_(fh2), fo_(fo) {}
    void addForces(Water const &w) {
      for (int a = 0; a < 3; ++a) {
        fh1_[a] += w.fh1_[a];
        fh2_[a] += w.fh2_[a];
        fo_[a] += w.fo_[a];
      };
    }
    const double *const rh1_;
    const double *const rh2_;
    const double *const ro_;
    const double *const rc_;
    double *const fh1_;
    double *const fh2_;
    double *const fo_;
  };

  /** Interactions within a molecules.
      CCL interaction with @ref amira2004 "Amira's grafting rules".
      @param[in, out] water   Add forces.
      @param[in, out] U Add energy to @a U.
      */
  void intramolecular(Water &water, double &U);

  /** Interactions between two molecules.
      Lennard Jones between oxygen only with cutoff
      @param[in, out] w1 Molecule of water 1.
      @param[in, out] w2 Molecule of water 2.
      @param[in, out] U Add the potential energy to @a U.
      @see intermolecularFull() and intermolecularSwitching().
      */
  void lennardJonesWithCutoff(Water &w1, Water &w2, double &U);

  /** Interactions between two molecules.
      Coulomb interaction between two molecules of water with molecules-based
     cutoff
      @param[in, out] w1 Molecule 1.
      @param[in, out] w2 Molecule 2.
      @param[in, out] U Add the potential energy to @a U.
      @see intermolecularFull() and intermolecularSwitching().
      */
  void coulombWithCutoff(Water &w1, Water &w2, double &U);

  /** Interactions between two molecules.
      Coulomb interaction between two molecules. Full interaction (i.e. no
     cutoff).
      @param[in, out] w1 Molecule 1.
      @param[in, out] w2 Molecule 2.
      @param[in, out] U Add the potential energy to @a U.
      @see intermolecularFull() and intermolecularSwitching().
      */
  void coulombFull(Water &w1, Water &w2, double &U);

  /// Distance OH
  static const double roh_;
  /// Angle HOH
  static const double theta_;
  /// Distance HH
  static const double rhh_;
  /// Charge on one hydrogen
  static const double charge_;
  /// Square of # charge_
  static const double charge2_;
  /// Lennard-Jones. See PotentialBase::sigma() and PotentialBase::epsilon() for
  /// definition
  static const double A_;
  /// Lennard-Jones. See PotentialBase::sigma() and PotentialBase::epsilon() for
  /// definition
  static const double B_;
  /// Lennard-Jones. See PotentialBase::lennardJones() for definition.
  static const double sigma_;
  /// Lennard-Jones. See PotentialBase::lennardJones() for definition.
  static const double epsilon_;
  /// Polarisation correction. Energy added per molecule (see @ref berendsen1987
  /// "original publication" for definition).
  static double const polarisationEnergy_;

private:
  /** Initialise Rho.
  Compute Rho coordinates according to @ref amira2004 "Amira's grafting rules":
  @f[
      \rho^{SPC+CCL}=\frac{r-r^{SPC}_e}{r-r^{SPC}_e+r^{CCL}_e}
      @f]
  @f[
      \frac{d\rho}{d\mathbf
  r}=\frac{-r_e^{CCL}}{\left(r-r^{SPC}_e+r^{CCL}_e\right)^2}\frac{\mathbf r}{r}
      @f]
  */
  void initialiseRho(Vector3 const &v, Rho &r);

  /** Compute forces and energy.
      The function calculates and returns the energy and forces applied on each
     atom. The function has preconditions with which the user must comply before
     calling this function. The periodic boundaries must set with
     (setPeriodicity()). The content of and arrays @a fh1 , @a fh2 , @a fo must
     be zero for all elements. \n Arrays ending in @em h1 are for the hydrogen
     one of the molecules and those ending in @em h2 for hydrogen 2.
      @param[in] nMolecules Number of molecules.
      @param[in] "rh1, rh2, ro" Positions of hydrogens and oxygens.
      @param[in,out] "fh1, fh2, fo" Forces on hydrogens and oxygens.
      @param[in,out] energy      Potential energy.
      @param[in] b      Periodic boundaries.
      @param[in] "xh1, xh2, xo" Tell which atoms are fixed and which are
     movable. This parameter are optional. When provided, interaction between
     fixed atoms are skipped. The three arrays @a xh1, @a xh2, @a xo must be
     provided for the optimisation to work.
      @warning Remember the preconditions.*/
  template <int H, int O>
  void computeTemplate(const int nMolecules, const double (*const rh1)[H * 3],
                       const double (*const rh2)[H * 3],
                       const double (*const ro)[O * 3],
                       double (*const fh1)[H * 3], double (*const fh2)[H * 3],
                       double (*const fo)[O * 3], double &energy,
                       double const b[], bool const (*const xh1)[H] = 0,
                       bool const (*const xh2)[H] = 0,
                       bool const (*const xo)[O] = 0);
};
} // namespace forcefields
#endif

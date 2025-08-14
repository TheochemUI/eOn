#ifndef ZBLPOT_H
#define ZBLPOT_H

#include "../../Potential.h"
#include <map>
#include <memory>
#include <vector>

// This is modified from the GNU GPL licensed LAMMPS codebase
// LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
// https://www.lammps.org/, Sandia National Laboratories
// LAMMPS development team: developers@lammps.org

// From J.F. Zeigler, J. P. Biersack and U. Littmark,
// "The Stopping and Range of Ions in Matter" volume 1, Pergamon, 1985.

class ZBLPot : public Potential {
public:
  /**
   * @brief Constructor for the LAMMPS-style ZBL Potential.
   * @param p A shared pointer to the eOn Parameters object
   */
  ZBLPot(std::shared_ptr<Parameters> p);

  /**
   * @brief Default virtual destructor.
   */
  ~ZBLPot() override = default;

  /**
   * @brief Calculates the total potential energy and forces for all atoms in
   * the system. This is the main computational method called repeatedly during
   * a simulation.
   */
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  // --- Private Methods ---

  /**
   * @brief Performs just-in-time initialization on the first force call.
   * Extracts unique elements from the system's atomic numbers and calculates
   * all pair coefficients needed for the simulation.
   */
  void initialize_for_system(long N, const int *atomicNrs);

  /**
   * @brief Calculates and stores the coefficients for a single pair of atom
   * types (i, j). This is a direct port of the logic from the LAMMPS
   * `set_coeff` function.
   */
  void set_pair_coeffs(int i, int j, double zi, double zj);

  // --- Direct ports of LAMMPS ZBL calculation functions ---
  double e_zbl(double r, int i, int j) const;
  double dzbldr(double r, int i, int j) const;
  double d2zbldr2(double r, int i, int j) const;

  // --- Member Variables ---
  // Flag for lazy initialization
  bool is_initialized;

  // Cutoffs
  double cut_inner;
  double cut_global;
  double cut_inner_sq;
  double cut_global_sq;

  // These are lazy initialized on the first force call
  // Lookup table from atomic number (e.g., 79) to internal type index (e.g., 0)
  std::map<int, int> atomic_nr_to_type;

  // Storage for parameters, indexed by internal type index
  int n_types;
  std::vector<double> z_values;

  // Pre-calculated coefficient tables, indexed by internal type indices
  std::vector<std::vector<double>> d1a, d2a, d3a, d4a;
  std::vector<std::vector<double>> zze;
  std::vector<std::vector<double>> sw1, sw2, sw3, sw4, sw5;

  // --- Universal ZBL Constants (from pair_zbl_const.h and LAMMPS units) ---
  static constexpr double C1 = 0.02817;
  static constexpr double C2 = 0.28022;
  static constexpr double C3 = 0.50986;
  static constexpr double C4 = 0.18175;

  // Slightly off wrt the equations, but matched with the LAMMPS source
  static constexpr double D1 = 3.19980;
  static constexpr double D2 = 0.94229;
  static constexpr double D3 = 0.40290;
  static constexpr double D4 = 0.20162;

  static constexpr double PZBL = 0.23;
  // Screening factor
  static constexpr double A0 = 0.46850;
  // Conversion factor in LAMMPS 'metal' units (eV*Angstrom)
  static constexpr double QQR2E = 14.399645;
};

#endif // ZBLPOT_H

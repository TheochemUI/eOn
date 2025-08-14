#ifndef ZBLPOT_H
#define ZBLPOT_H

#include "../../Potential.h"
#include <map>
#include <memory>
#include <vector>

/**
 * @brief LAMMPS-style ZBL Potential
 *        Based on: J.F. Ziegler, J.P. Biersack, U. Littmark,
 *        "Stopping and Range of Ions in Matter", Pergamon (1985).
 *
 * Directly adapted from the GNU GPL licensed LAMMPS codebase.
 */
class ZBLPot : public Potential {
public:
  ZBLPot(std::shared_ptr<Parameters> p);

  ~ZBLPot() override = default;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  // --- Private Methods ---

  /// Just-in-time initialization for the current system
  void initialize_for_system(long N, const int *atomicNrs);

  /// Set coefficients for a pair of atom types
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

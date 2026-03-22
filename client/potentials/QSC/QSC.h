/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once

#include <vector>

#ifndef QSC_STANDALONE
#include "../../Potential.h"
#endif

/// Quantum Sutton-Chen (QSC) potential -- EAM-type with:
///   F_i(rho_i) = c * sqrt(rho_i)
///   rho_i = sum_(j!=i) (a/r_ij)^m
///   V(r_ij) = (a/r_ij)^n
class QSC
#ifndef QSC_STANDALONE
    : public Potential
#endif
{
public:
  explicit QSC(const Parameters &params)
      : Potential(PotType::QSC, params) {
    int i = 0;
    while (qsc_default_params[i].Z != -1) {
      qsc_params_.push_back(qsc_default_params[i]);
      i++;
    }
  }
  QSC() : QSC(Parameters{}) {}
  ~QSC() override = default;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
  void set_verlet_skin(double dr);
  void set_cutoff(double c);
  [[nodiscard]] double get_cutoff() const;
  void set_qsc_parameter(int Z, double n, double m, double epsilon, double c,
                         double a);

  long vlist_updates{0};

private:
  struct qsc_parameters {
    int Z;
    double n;
    double m;
    double epsilon;
    double c;
    double a;
  };

  struct distance {
    double d[3];
    double r;
  };

  bool init_{false};
  long N_{0};
  double cutoff_{8.0};
  double verlet_skin_{0.5};

  // Flattened N*N 2D arrays (row-major)
  std::vector<distance> distances_; // [N*N]
  std::vector<int> vlist_;          // [N*N] neighbor indices
  std::vector<int> nlist_;          // [N] neighbor counts
  std::vector<double> oldR_;        // [3*N]
  std::vector<double> rho_;         // [N]
  std::vector<double> sqrtrho_;     // [N]
  std::vector<double> V_;           // [N*N]
  std::vector<double> phi_;         // [N*N]

  void initialize(long N);
  void energy(long N, const double *R, const int *atomicNrs, double *U,
              const double *box);
  void new_vlist(long N, const double *R, const double *box);
  void update_distances(long N, const double *R, const double *box);
  [[nodiscard]] bool verlet_needs_update(long N, const double *R,
                                         const double *box) const;
  void calc_distance(const double *box, const double *R1, int i,
                     const double *R2, int j, distance *d) const;

  static const qsc_parameters qsc_default_params[];
  std::vector<qsc_parameters> qsc_params_;

  [[nodiscard]] qsc_parameters get_qsc_parameters(int a, int b) const;
  [[nodiscard]] static double dpowi(double x, unsigned n);
  [[nodiscard]] static double pair_potential(double r, double a, double n);
};

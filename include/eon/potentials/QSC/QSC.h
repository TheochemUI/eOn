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
#include "eon/Potential.h"
#endif
#include "eon/VesinNeighbors.h"

/// Quantum Sutton-Chen (QSC) potential -- EAM-type with:
///   F_i(rho_i) = c * sqrt(rho_i)
///   rho_i = sum_(j!=i) (a/r_ij)^m
///   V(r_ij) = (a/r_ij)^n
///
/// Neighbor pairs come from the thread-local ``eonc::PairListCache`` (Verlet
/// skin on top of vesin), so repeated force evaluations on nearby geometries
/// skip the list rebuild without racing under NEB.
class QSC
#ifndef QSC_STANDALONE
    : public Potential
#endif
{
public:
  explicit QSC(const Parameters &params) : Potential(PotType::QSC, params) {
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
  /// rho_/sqrtrho_ are instance state written during force.
  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return true;
  }
  /// Verlet skin for the cached pair list (Angstrom).
  void set_verlet_skin(double dr);
  void set_cutoff(double c);
  [[nodiscard]] double get_cutoff() const;
  void set_qsc_parameter(int Z, double n, double m, double epsilon, double c,
                         double a);

  /// Number of force evaluations that consulted the pair-list cache.
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

  double cutoff_{8.0};
  double skin_{1.0};

  std::vector<double> rho_;     // [N]
  std::vector<double> sqrtrho_; // [N]

  /// Energy from ``nl`` (caller must have computed it on this thread).
  void energy_from_nl(long N, const double *R, const int *atomicNrs, double *U,
                      const eonc::CachedPairList &nl);

  static const qsc_parameters qsc_default_params[];
  std::vector<qsc_parameters> qsc_params_;

  [[nodiscard]] qsc_parameters get_qsc_parameters(int a, int b) const;
  [[nodiscard]] static double dpowi(double x, unsigned n);
  [[nodiscard]] static double pair_potential(double r, double a, double n);
};

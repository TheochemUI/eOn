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

/// @file QSC.cpp
/// @brief Quantum Sutton-Chen potential implementation.
///
/// EAM-type potential with density = (a/r)^m, pair = (a/r)^n,
/// embedding = c * epsilon * sqrt(rho). Neighbor pairs via a thread_local
/// eonc::VesinNeighbors (allocation reuse; safe under NEB parallel).

#include "eon/potentials/QSC/QSC.h"
#include "eon/potentials/QSC/Parameters.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

void QSC::energy_from_nl(long N, const int *atomicNrs, double *U,
                         const eonc::VesinNeighbors &nl) {
  *U = 0.0;
  rho_.assign(static_cast<std::size_t>(N), 0.0);
  sqrtrho_.assign(static_cast<std::size_t>(N), 0.0);
  if (N < 2 || nl.size() == 0) {
    return;
  }

  std::vector<double> pair_term(static_cast<std::size_t>(N), 0.0);

  for (std::size_t p = 0; p < nl.size(); ++p) {
    const long i = static_cast<long>(nl.i(p));
    const long j = static_cast<long>(nl.j(p));
    if (i == j) {
      continue;
    }
    const double r_ij = nl.distance(p);
    if (r_ij <= 0.0 || r_ij > cutoff_) {
      continue;
    }

    const auto p_ii = get_qsc_parameters(atomicNrs[i], atomicNrs[i]);
    const auto p_ij = get_qsc_parameters(atomicNrs[i], atomicNrs[j]);
    const auto p_jj = get_qsc_parameters(atomicNrs[j], atomicNrs[j]);

    const double phi_ij = pair_potential(r_ij, p_jj.a, p_jj.m);
    const double phi_ji = pair_potential(r_ij, p_ii.a, p_ii.m);
    rho_[static_cast<std::size_t>(i)] += phi_ij;
    rho_[static_cast<std::size_t>(j)] += phi_ji;

    const double V = p_ij.epsilon * pair_potential(r_ij, p_ij.a, p_ij.n);
    pair_term[static_cast<std::size_t>(i)] += V;
  }

  for (long i = 0; i < N; i++) {
    const auto p_ii_e = get_qsc_parameters(atomicNrs[i], atomicNrs[i]);
    sqrtrho_[static_cast<std::size_t>(i)] =
        std::sqrt(rho_[static_cast<std::size_t>(i)]);
    const double embedding =
        p_ii_e.c * p_ii_e.epsilon * sqrtrho_[static_cast<std::size_t>(i)];
    *U += pair_term[static_cast<std::size_t>(i)] - embedding;
  }
}

void QSC::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, double *variance, const double *box) {
  variance = nullptr;
  for (long i = 0; i < 3 * N; i++) {
    F[i] = 0.0;
  }
  *U = 0.0;
  if (N < 2 || cutoff_ <= 0.0) {
    return;
  }

  thread_local eonc::VesinNeighbors nl;
  eonc::VesinNeighbors::Options opt;
  opt.cutoff = cutoff_;
  opt.full = false;
  opt.return_distances = true;
  opt.return_vectors = true;
  nl.compute(R, static_cast<std::size_t>(N), box, opt);
  ++vlist_updates;

  energy_from_nl(N, atomicNrs, U, nl);

  for (std::size_t p = 0; p < nl.size(); ++p) {
    const long i = static_cast<long>(nl.i(p));
    const long j = static_cast<long>(nl.j(p));
    if (i == j) {
      continue;
    }
    const double r_ij = nl.distance(p);
    if (r_ij <= 0.0 || r_ij > cutoff_) {
      continue;
    }

    const auto p_ii = get_qsc_parameters(atomicNrs[i], atomicNrs[i]);
    const auto p_ij = get_qsc_parameters(atomicNrs[i], atomicNrs[j]);
    const auto p_jj = get_qsc_parameters(atomicNrs[j], atomicNrs[j]);

    const double phi_ij = pair_potential(r_ij, p_jj.a, p_jj.m);
    const double phi_ji = pair_potential(r_ij, p_ii.a, p_ii.m);
    const double V = p_ij.epsilon * pair_potential(r_ij, p_ij.a, p_ij.n);

    double Fij = p_ij.n * V;
    Fij -= p_ii.epsilon * p_ii.c * p_jj.m * 0.5 *
           (1.0 / sqrtrho_[static_cast<std::size_t>(i)]) * phi_ij;
    Fij -= p_jj.epsilon * p_jj.c * p_ii.m * 0.5 *
           (1.0 / sqrtrho_[static_cast<std::size_t>(j)]) * phi_ji;
    Fij /= r_ij;

    const double *v = nl.vector(p);
    const double fscale = Fij / r_ij;
    const double fx = fscale * (-v[0]);
    const double fy = fscale * (-v[1]);
    const double fz = fscale * (-v[2]);

    F[3 * i] += fx;
    F[3 * i + 1] += fy;
    F[3 * i + 2] += fz;
    F[3 * j] -= fx;
    F[3 * j + 1] -= fy;
    F[3 * j + 2] -= fz;
  }
}

double QSC::dpowi(double x, unsigned n) {
  double p = x;
  double r = 1.0;
  while (n > 0) {
    if (n % 2 == 1)
      r *= p;
    p *= p;
    n /= 2;
  }
  return r;
}

double QSC::pair_potential(double r, double a, double n) {
  double x = a / r;
  if ((n - std::floor(n)) == 0.0 && n > 0) {
    return dpowi(x, static_cast<unsigned>(n));
  }
  return std::pow(x, n);
}

void QSC::set_verlet_skin(double dr) {
  assert(dr > 0.0);
  (void)dr;
}

void QSC::set_cutoff(double c) {
  assert(c > 0.0);
  cutoff_ = c;
}

double QSC::get_cutoff() const { return cutoff_; }

void QSC::set_qsc_parameter(int Z, double n, double m, double epsilon, double c,
                            double a) {
  qsc_parameters p{Z, n, m, epsilon, c, a};
  for (auto &existing : qsc_params_) {
    if (existing.Z == Z) {
      existing = p;
      return;
    }
  }
  qsc_params_.push_back(p);
}

QSC::qsc_parameters QSC::get_qsc_parameters(int element_a,
                                            int element_b) const {
  int ia = -1, ib = -1;
  for (size_t i = 0; i < qsc_params_.size(); i++) {
    if (element_a == qsc_params_[i].Z)
      ia = static_cast<int>(i);
    if (element_b == qsc_params_[i].Z)
      ib = static_cast<int>(i);
    if (ia != -1 && ib != -1)
      break;
  }

  if (ia == -1) {
    std::fprintf(stderr, "ERROR: QSC lacks parameters for element %i\n",
                 element_a);
    throw 1;
  }
  if (ib == -1) {
    std::fprintf(stderr, "ERROR: QSC lacks parameters for element %i\n",
                 element_b);
    throw 1;
  }

  if (ia == ib)
    return qsc_params_[ia];

  return qsc_parameters{
      0,
      0.5 * (qsc_params_[ia].n + qsc_params_[ib].n),
      0.5 * (qsc_params_[ia].m + qsc_params_[ib].m),
      std::sqrt(qsc_params_[ia].epsilon * qsc_params_[ib].epsilon),
      qsc_params_[ia].c,
      0.5 * (qsc_params_[ia].a + qsc_params_[ib].a),
  };
}

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
/// embedding = c * epsilon * sqrt(rho). Uses Verlet neighbor list
/// for O(N) force evaluation.

#include "QSC.h"
#include "Parameters.h"

#include <cassert>
#include <cmath>
#include <cstdio>

void QSC::initialize(long N) {
  N_ = N;
  rho_.assign(N, 0.0);
  sqrtrho_.assign(N, 0.0);
  nlist_.assign(N, 0);
  oldR_.assign(3 * N, 0.0);
  // Flattened N*N arrays
  distances_.resize(N * N);
  vlist_.assign(N * N, 0);
  V_.assign(N * N, 0.0);
  phi_.assign(N * N, 0.0);
  init_ = true;
}

void QSC::new_vlist(long N, const double *R, const double *box) {
  double rv = cutoff_ + verlet_skin_;

  for (long i = 0; i < N; i++) {
    nlist_[i] = 0;
    for (long j = i + 1; j < N; j++) {
      calc_distance(box, R, i, R, j, &distances_[i * N + j]);
      if (distances_[i * N + j].r <= rv) {
        vlist_[i * N + nlist_[i]] = static_cast<int>(j);
        nlist_[i]++;
      }
    }
  }

  for (long i = 0; i < 3 * N; i++) {
    oldR_[i] = R[i];
  }
  vlist_updates++;
}

void QSC::update_distances(long N, const double *R, const double *box) {
  for (long i = 0; i < N; i++) {
    for (int k = 0; k < nlist_[i]; k++) {
      int j = vlist_[i * N + k];
      calc_distance(box, R, i, R, j, &distances_[i * N + j]);
    }
  }
}

bool QSC::verlet_needs_update(long N, const double *R,
                              const double *box) const {
  distance diff;
  double dist_max1 = 0.0;
  double dist_max2 = 0.0;

  for (long i = 0; i < N; i++) {
    calc_distance(box, oldR_.data(), i, R, i, &diff);
    if (diff.r > dist_max1) {
      dist_max2 = dist_max1;
      dist_max1 = diff.r;
    } else if (diff.r > dist_max2) {
      dist_max2 = diff.r;
    }
    if (dist_max1 + dist_max2 > verlet_skin_) {
      return true;
    }
  }
  return false;
}

void QSC::energy(long N, const double *R, const int *atomicNrs, double *U,
                 const double *box) {
  *U = 0.0;
  for (long i = 0; i < N; i++) {
    rho_[i] = 0.0;
  }

  int prev_i_Z = -1, prev_j_Z = -1;
  for (long i = 0; i < N; i++) {
    double pair_term = 0.0;
    qsc_parameters p_ii{};

    if (prev_i_Z != atomicNrs[i]) {
      p_ii = get_qsc_parameters(atomicNrs[i], atomicNrs[i]);
    }
    prev_i_Z = atomicNrs[i];

    for (int k = 0; k < nlist_[i]; k++) {
      qsc_parameters p_ij{}, p_jj{};
      int j = vlist_[i * N + k];
      double r_ij = distances_[i * N + j].r;
      if (r_ij > cutoff_) continue;

      if (prev_j_Z != atomicNrs[j]) {
        p_ij = get_qsc_parameters(atomicNrs[i], atomicNrs[j]);
        p_jj = get_qsc_parameters(atomicNrs[j], atomicNrs[j]);
      }
      prev_j_Z = atomicNrs[j];

      // Density contributions
      phi_[i * N + j] = pair_potential(r_ij, p_jj.a, p_jj.m);
      rho_[i] += phi_[i * N + j];
      phi_[j * N + i] = pair_potential(r_ij, p_ii.a, p_ii.m);
      rho_[j] += phi_[j * N + i];

      // Repulsive pair term
      V_[i * N + j] = p_ij.epsilon * pair_potential(r_ij, p_ij.a, p_ij.n);
      pair_term += V_[i * N + j];
    }
    sqrtrho_[i] = std::sqrt(rho_[i]);
    double embedding = p_ii.c * p_ii.epsilon * sqrtrho_[i];
    *U += pair_term - embedding;
  }
}

void QSC::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, double *variance, const double *box) {
  variance = nullptr;
  if (!init_) {
    initialize(N);
    new_vlist(N, R, box);
  } else if (N != N_) {
    initialize(N);
    new_vlist(N, R, box);
  } else if (verlet_needs_update(N, R, box)) {
    new_vlist(N, R, box);
  } else {
    update_distances(N, R, box);
  }

  energy(N, R, atomicNrs, U, box);

  for (long i = 0; i < 3 * N; i++) {
    F[i] = 0.0;
  }

  int prev_i_Z = -1, prev_j_Z = -1;
  for (long i = 0; i < N; i++) {
    qsc_parameters p_ii{};
    if (prev_i_Z != atomicNrs[i]) {
      p_ii = get_qsc_parameters(atomicNrs[i], atomicNrs[i]);
    }
    prev_i_Z = atomicNrs[i];

    for (int k = 0; k < nlist_[i]; k++) {
      qsc_parameters p_ij{}, p_jj{};
      int j = vlist_[i * N + k];
      double r_ij = distances_[i * N + j].r;
      if (r_ij > cutoff_) continue;

      if (prev_j_Z != atomicNrs[j]) {
        p_ij = get_qsc_parameters(atomicNrs[i], atomicNrs[j]);
        p_jj = get_qsc_parameters(atomicNrs[j], atomicNrs[j]);
      }
      prev_j_Z = atomicNrs[j];

      double Fij = p_ij.n * V_[i * N + j];
      Fij -= p_ii.epsilon * p_ii.c * p_jj.m * 0.5 * (1.0 / sqrtrho_[i]) *
             phi_[i * N + j];
      Fij -= p_jj.epsilon * p_jj.c * p_ii.m * 0.5 * (1.0 / sqrtrho_[j]) *
             phi_[j * N + i];
      Fij /= r_ij;

      const auto &dist = distances_[i * N + j];
      double fscale = Fij / r_ij;
      double fx = fscale * dist.d[0];
      double fy = fscale * dist.d[1];
      double fz = fscale * dist.d[2];

      F[3 * i] += fx;
      F[3 * i + 1] += fy;
      F[3 * i + 2] += fz;
      F[3 * j] -= fx;
      F[3 * j + 1] -= fy;
      F[3 * j + 2] -= fz;
    }
  }
}

double QSC::dpowi(double x, unsigned n) {
  double p = x;
  double r = 1.0;
  while (n > 0) {
    if (n % 2 == 1) r *= p;
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

void QSC::calc_distance(const double *box, const double *R1, int i,
                        const double *R2, int j, distance *d) const {
  double dx = R1[3 * i] - R2[3 * j];
  double dy = R1[3 * i + 1] - R2[3 * j + 1];
  double dz = R1[3 * i + 2] - R2[3 * j + 2];

  // Orthogonal PBC
  dx -= box[0] * std::floor(dx / box[0] + 0.5);
  dy -= box[4] * std::floor(dy / box[4] + 0.5);
  dz -= box[8] * std::floor(dz / box[8] + 0.5);

  d->r = std::sqrt(dx * dx + dy * dy + dz * dz);
  d->d[0] = dx;
  d->d[1] = dy;
  d->d[2] = dz;
}

void QSC::set_verlet_skin(double dr) {
  assert(dr > 0.0);
  verlet_skin_ = dr;
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
    if (element_a == qsc_params_[i].Z) ia = static_cast<int>(i);
    if (element_b == qsc_params_[i].Z) ib = static_cast<int>(i);
    if (ia != -1 && ib != -1) break;
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

  if (ia == ib) return qsc_params_[ia];

  // Mixing rules
  return qsc_parameters{
      0,
      0.5 * (qsc_params_[ia].n + qsc_params_[ib].n),
      0.5 * (qsc_params_[ia].m + qsc_params_[ib].m),
      std::sqrt(qsc_params_[ia].epsilon * qsc_params_[ib].epsilon),
      qsc_params_[ia].c,
      0.5 * (qsc_params_[ia].a + qsc_params_[ib].a),
  };
}

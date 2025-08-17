#include "ZBLPot.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

ZBLPot::ZBLPot(std::shared_ptr<Parameters> p)
    : Potential(p),
      is_initialized(false) {
  cut_inner = p->zbl_options.cut_inner;
  cut_global = p->zbl_options.cut_global;

  if (cut_inner <= 0.0 || cut_inner >= cut_global) {
    throw std::runtime_error(
        "Invalid ZBL cutoffs: require 0.0 < cut_inner < cut_global.");
  }
  cut_inner_sq = cut_inner * cut_inner;
  cut_global_sq = cut_global * cut_global;
}

void ZBLPot::initialize_for_system(long N, const int *atomicNrs) {
  std::set<int> unique_z(atomicNrs, atomicNrs + N);
  n_types = static_cast<int>(unique_z.size());
  if (n_types == 0)
    return;

  z_values.resize(n_types);
  int idx = 0;
  for (int z : unique_z) {
    atomic_nr_to_type[z] = idx;
    z_values[idx++] = static_cast<double>(z);
  }

  auto resize_tbl = [&](auto &tbl) {
    tbl.assign(n_types, std::vector<double>(n_types, 0.0));
  };
  resize_tbl(d1a);
  resize_tbl(d2a);
  resize_tbl(d3a);
  resize_tbl(d4a);
  resize_tbl(zze);
  resize_tbl(sw1);
  resize_tbl(sw2);
  resize_tbl(sw3);
  resize_tbl(sw4);
  resize_tbl(sw5);

  for (int i = 0; i < n_types; ++i)
    for (int j = i; j < n_types; ++j)
      set_pair_coeffs(i, j, z_values[i], z_values[j]);
}

void ZBLPot::set_pair_coeffs(int i, int j, double zi, double zj) {
  const double a_inv = (std::pow(zi, PZBL) + std::pow(zj, PZBL)) / A0;

  d1a[i][j] = d1a[j][i] = D1 * a_inv;
  d2a[i][j] = d2a[j][i] = D2 * a_inv;
  d3a[i][j] = d3a[j][i] = D3 * a_inv;
  d4a[i][j] = d4a[j][i] = D4 * a_inv;

  zze[i][j] = zze[j][i] = zi * zj * QQR2E;

  const double tc = cut_global - cut_inner;
  const double fc = e_zbl(cut_global, i, j);
  const double fcp = dzbldr(cut_global, i, j);
  const double fcpp = d2zbldr2(cut_global, i, j);
  const double swa = (-3.0 * fcp + tc * fcpp) / (tc * tc);
  const double swb = (2.0 * fcp - tc * fcpp) / (tc * tc * tc);
  const double swc = -fc + (tc / 2.0) * fcp - (tc * tc / 12.0) * fcpp;

  sw1[i][j] = sw1[j][i] = swa;
  sw2[i][j] = sw2[j][i] = swb;
  sw3[i][j] = sw3[j][i] = swa / 3.0;
  sw4[i][j] = sw4[j][i] = swb / 4.0;
  sw5[i][j] = sw5[j][i] = swc;
}

void ZBLPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double * /*box*/) {
  if (!is_initialized) {
    initialize_for_system(N, atomicNrs);
    is_initialized = true;
  }

  if (variance)
    *variance = 0.0;
  *U = 0.0;
  std::fill(F, F + 3 * N, 0.0);

  for (long i = 0; i < N; ++i) {
    for (long j = i + 1; j < N; ++j) {
      const double dx = R[3 * j] - R[3 * i];
      const double dy = R[3 * j + 1] - R[3 * i + 1];
      const double dz = R[3 * j + 2] - R[3 * i + 2];
      const double r_sq = dx * dx + dy * dy + dz * dz;

      if (r_sq >= cut_global_sq)
        continue;

      const double r = std::sqrt(r_sq);
      const int type_i = atomic_nr_to_type.at(atomicNrs[i]);
      const int type_j = atomic_nr_to_type.at(atomicNrs[j]);

      double energy_pair = e_zbl(r, type_i, type_j) + sw5[type_i][type_j];
      if (r_sq > cut_inner_sq) {
        const double t = r - cut_inner;
        energy_pair +=
            t * t * t * (sw3[type_i][type_j] + sw4[type_i][type_j] * t);
      }
      *U += energy_pair;

      double fpair_term = dzbldr(r, type_i, type_j);
      if (r_sq > cut_inner_sq) {
        const double t = r - cut_inner;
        fpair_term += t * t * (sw1[type_i][type_j] + sw2[type_i][type_j] * t);
      }
      const double f_div_r = -fpair_term / r;

      const double fx = dx * f_div_r;
      const double fy = dy * f_div_r;
      const double fz = dz * f_div_r;

      F[3 * i] -= fx;
      F[3 * i + 1] -= fy;
      F[3 * i + 2] -= fz;
      F[3 * j] += fx;
      F[3 * j + 1] += fy;
      F[3 * j + 2] += fz;
    }
  }
}

double ZBLPot::e_zbl(double r, int i, int j) const {
  const double r_inv = 1.0 / r;
  const double sum =
      C4 * std::exp(-d1a[i][j] * r) + C3 * std::exp(-d2a[i][j] * r) +
      C2 * std::exp(-d3a[i][j] * r) + C1 * std::exp(-d4a[i][j] * r);

  return zze[i][j] * sum * r_inv;
}

double ZBLPot::dzbldr(double r, int i, int j) const {
  const double r_inv = 1.0 / r;
  const double e1 = std::exp(-d1a[i][j] * r);
  const double e2 = std::exp(-d2a[i][j] * r);
  const double e3 = std::exp(-d3a[i][j] * r);
  const double e4 = std::exp(-d4a[i][j] * r);

  const double sum = C4 * e1 + C3 * e2 + C2 * e3 + C1 * e4;
  const double sum_p = -C4 * d1a[i][j] * e1 - C3 * d2a[i][j] * e2 -
                       C2 * d3a[i][j] * e3 - C1 * d4a[i][j] * e4;

  return zze[i][j] * (sum_p - sum * r_inv) * r_inv;
}

double ZBLPot::d2zbldr2(double r, int i, int j) const {
  const double r_inv = 1.0 / r;
  const double d1 = d1a[i][j], d2 = d2a[i][j], d3 = d3a[i][j], d4 = d4a[i][j];

  const double e1 = std::exp(-d1 * r);
  const double e2 = std::exp(-d2 * r);
  const double e3 = std::exp(-d3 * r);
  const double e4 = std::exp(-d4 * r);

  const double sum = C4 * e1 + C3 * e2 + C2 * e3 + C1 * e4;
  const double sum_p =
      C4 * e1 * d1 + C3 * e2 * d2 + C2 * e3 * d3 + C1 * e4 * d4;
  const double sum_pp = C4 * e1 * d1 * d1 + C3 * e2 * d2 * d2 +
                        C2 * e3 * d3 * d3 + C1 * e4 * d4 * d4;

  return zze[i][j] *
         (sum_pp + 2.0 * sum_p * r_inv + 2.0 * sum * r_inv * r_inv) * r_inv;
}

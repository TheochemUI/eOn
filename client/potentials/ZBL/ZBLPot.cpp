#include "ZBLPot.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

ZBLPot::ZBLPot(std::shared_ptr<Parameters> p)
    : Potential(p),
      is_initialized(false) {
  cut_inner = p->cut_inner;
  cut_global = p->cut_global;

  if (cut_inner <= 0.0 || cut_inner >= cut_global) {
    throw std::runtime_error(
        "Invalid ZBL cutoffs: must be 0.0 < cut_inner < cut_global.");
  }
  cut_inner_sq = cut_inner * cut_inner;
  cut_global_sq = cut_global * cut_global;
}

void ZBLPot::initialize_for_system(long N, const int *atomicNrs) {
  // 1. Find the unique atomic numbers from the provided system
  std::set<int> unique_elements_set(atomicNrs, atomicNrs + N);
  std::vector<int> unique_elements(unique_elements_set.begin(),
                                   unique_elements_set.end());

  if (unique_elements.empty()) {
    n_types = 0;
    return;
  }
  n_types = unique_elements.size();

  // 2. Populate the lookup table and Z-value vector
  z_values.resize(n_types);
  int type_index = 0;
  for (const auto &atomic_nr : unique_elements) {
    atomic_nr_to_type[atomic_nr] = type_index;
    z_values[type_index] = static_cast<double>(atomic_nr);
    type_index++;
  }

  // 3. Pre-calculate all pair coefficients
  auto resize_vec = [&](std::vector<std::vector<double>> &vec) {
    vec.assign(n_types, std::vector<double>(n_types, 0.0));
  };
  resize_vec(d1a);
  resize_vec(d2a);
  resize_vec(d3a);
  resize_vec(d4a);
  resize_vec(zze);
  resize_vec(sw1);
  resize_vec(sw2);
  resize_vec(sw3);
  resize_vec(sw4);
  resize_vec(sw5);

  for (int i = 0; i < n_types; ++i) {
    for (int j = i; j < n_types; ++j) {
      set_pair_coeffs(i, j, z_values[i], z_values[j]);
    }
  }
}

void ZBLPot::set_pair_coeffs(int i, int j, double zi, double zj) {
  // This logic is ported directly from LAMMPS's set_coeff function
  double a_inv = (pow(zi, PZBL) + pow(zj, PZBL)) / A0;

  d1a[i][j] = D1 * a_inv;
  d2a[i][j] = D2 * a_inv;
  d3a[i][j] = D3 * a_inv;
  d4a[i][j] = D4 * a_inv;
  zze[i][j] = zi * zj * QQR2E;

  // Ensure symmetry for easier lookup
  d1a[j][i] = d1a[i][j];
  d2a[j][i] = d2a[i][j];
  d3a[j][i] = d3a[i][j];
  d4a[j][i] = d4a[i][j];
  zze[j][i] = zze[i][j];

  // Calculate switching function coefficients
  double tc = cut_global - cut_inner;
  double fc = e_zbl(cut_global, i, j);
  double fcp = dzbldr(cut_global, i, j);
  double fcpp = d2zbldr2(cut_global, i, j);

  double swa = (-3.0 * fcp + tc * fcpp) / (tc * tc);
  double swb = (2.0 * fcp - tc * fcpp) / (tc * tc * tc);
  double swc = -fc + (tc / 2.0) * fcp - (tc * tc / 12.0) * fcpp;

  sw1[i][j] = swa;
  sw2[i][j] = swb;
  sw3[i][j] = swa / 3.0;
  sw4[i][j] = swb / 4.0;
  sw5[i][j] = swc;

  // Ensure symmetry
  sw1[j][i] = sw1[i][j];
  sw2[j][i] = sw2[i][j];
  sw3[j][i] = sw3[i][j];
  sw4[j][i] = sw4[i][j];
  sw5[j][i] = sw5[i][j];
}

void ZBLPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  // --- Just-In-Time Initialization ---
  if (!is_initialized) {
    initialize_for_system(N, atomicNrs);
    is_initialized = true;
  }
  if (variance) {
    *variance = 0.0;
  }
  *U = 0.0;
  std::fill(F, F + 3 * N, 0.0);

  // Main pairwise loop
  for (long i = 0; i < N; ++i) {
    for (long j = i + 1; j < N; ++j) {
      double dx = R[3 * j] - R[3 * i];
      double dy = R[3 * j + 1] - R[3 * i + 1];
      double dz = R[3 * j + 2] - R[3 * i + 2];

      // NOTE: Assuming minimum image convention is handled by the caller (e.g.,
      // Matter class)

      double r_sq = dx * dx + dy * dy + dz * dz;

      if (r_sq < cut_global_sq) {
        double r = sqrt(r_sq);
        int type_i = atomic_nr_to_type.at(atomicNrs[i]);
        int type_j = atomic_nr_to_type.at(atomicNrs[j]);

        // --- Energy Calculation ---
        double energy_pair = e_zbl(r, type_i, type_j) + sw5[type_i][type_j];
        if (r_sq > cut_inner_sq) {
          double t = r - cut_inner;
          energy_pair +=
              t * t * t * (sw3[type_i][type_j] + sw4[type_i][type_j] * t);
        }
        *U += energy_pair;

        // --- Force Calculation ---
        double fpair_term = dzbldr(r, type_i, type_j);
        if (r_sq > cut_inner_sq) {
          double t = r - cut_inner;
          fpair_term += t * t * (sw1[type_i][type_j] + sw2[type_i][type_j] * t);
        }

        double f_div_r = -fpair_term / r;

        double fx = dx * f_div_r;
        double fy = dy * f_div_r;
        double fz = dz * f_div_r;

        F[3 * i] -= fx;
        F[3 * i + 1] -= fy;
        F[3 * i + 2] -= fz;
        F[3 * j] += fx;
        F[3 * j + 1] += fy;
        F[3 * j + 2] += fz;
      }
    }
  }
}

// The following three functions are direct ports of the LAMMPS calculations
double ZBLPot::e_zbl(double r, int i, int j) const {
  double r_inv = 1.0 / r;

  // The correct pairing is: C4-d1a, C3-d2a, C2-d3a, C1-d4a
  double sum = C4 * exp(-d1a[i][j] * r) + C3 * exp(-d2a[i][j] * r) +
               C2 * exp(-d3a[i][j] * r) + C1 * exp(-d4a[i][j] * r);

  return zze[i][j] * sum * r_inv;
}

double ZBLPot::dzbldr(double r, int i, int j) const {
  double r_inv = 1.0 / r;

  // Pre-calculate exponentials for clarity
  double e1 = exp(-d1a[i][j] * r); // Term with D1
  double e2 = exp(-d2a[i][j] * r); // Term with D2
  double e3 = exp(-d3a[i][j] * r); // Term with D3
  double e4 = exp(-d4a[i][j] * r); // Term with D4

  // Correct pairing for the sum
  double sum = C4 * e1 + C3 * e2 + C2 * e3 + C1 * e4;

  // Correct pairing for the derivative of the sum
  double sum_p = -C4 * d1a[i][j] * e1 - C3 * d2a[i][j] * e2 -
                 C2 * d3a[i][j] * e3 - C1 * d4a[i][j] * e4;

  return zze[i][j] * (sum_p - sum * r_inv) * r_inv;
}

double ZBLPot::d2zbldr2(double r, int i, int j) const {
  double r_inv = 1.0 / r;

  // Use local variables for dXa coefficients for readability
  double d1 = d1a[i][j];
  double d2 = d2a[i][j];
  double d3 = d3a[i][j];
  double d4 = d4a[i][j];

  // Pre-calculate exponentials
  double e1 = exp(-d1 * r);
  double e2 = exp(-d2 * r);
  double e3 = exp(-d3 * r);
  double e4 = exp(-d4 * r);

  // Correct pairing for sum, first derivative term, and second derivative term
  double sum = C4 * e1 + C3 * e2 + C2 * e3 + C1 * e4;
  double sum_p = C4 * e1 * d1 + C3 * e2 * d2 + C2 * e3 * d3 + C1 * e4 * d4;
  double sum_pp = C4 * e1 * d1 * d1 + C3 * e2 * d2 * d2 + C2 * e3 * d3 * d3 +
                  C1 * e4 * d4 * d4;

  return zze[i][j] *
         (sum_pp + 2.0 * sum_p * r_inv + 2.0 * sum * r_inv * r_inv) * r_inv;
}

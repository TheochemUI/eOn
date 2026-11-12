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
#include "eon/SetflEam.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace eonc {
namespace {

struct Spline {
  double origin{0.0};
  double step{1.0};
  std::vector<double> y;
  std::vector<double> y2;

  void build() {
    const int n = static_cast<int>(y.size());
    y2.assign(static_cast<size_t>(n), 0.0);
    if (n < 3) {
      return;
    }
    std::vector<double> u(static_cast<size_t>(n), 0.0);
    y2[0] = 0.0;
    u[0] = 0.0;
    for (int i = 1; i < n - 1; ++i) {
      const double sig = 0.5;
      const double p = sig * y2[static_cast<size_t>(i - 1)] + 2.0;
      y2[static_cast<size_t>(i)] = (sig - 1.0) / p;
      const double d =
          (y[static_cast<size_t>(i + 1)] - y[static_cast<size_t>(i)]) / step -
          (y[static_cast<size_t>(i)] - y[static_cast<size_t>(i - 1)]) / step;
      u[static_cast<size_t>(i)] =
          (6.0 * d / (2.0 * step) - sig * u[static_cast<size_t>(i - 1)]) / p;
    }
    y2[static_cast<size_t>(n - 1)] = 0.0;
    for (int k = n - 2; k >= 0; --k) {
      y2[static_cast<size_t>(k)] =
          y2[static_cast<size_t>(k)] * y2[static_cast<size_t>(k + 1)] +
          u[static_cast<size_t>(k)];
    }
  }

  void eval(double x, double &value, double &derivative) const {
    const int n = static_cast<int>(y.size());
    if (n == 0 || !(step > 0.0)) {
      value = 0.0;
      derivative = 0.0;
      return;
    }
    double t = (x - origin) / step;
    if (t < 0.0) {
      t = 0.0;
    }
    if (t > static_cast<double>(n - 1)) {
      t = static_cast<double>(n - 1);
    }
    const int i = std::min(static_cast<int>(t), n - 2);
    const double a = t - static_cast<double>(i);
    const double b = 1.0 - a;
    const double yi = y[static_cast<size_t>(i)];
    const double yj = y[static_cast<size_t>(i + 1)];
    const double y2i = y2[static_cast<size_t>(i)];
    const double y2j = y2[static_cast<size_t>(i + 1)];
    value = b * yi + a * yj +
            ((b * b * b - b) * y2i + (a * a * a - a) * y2j) * (step * step) /
                6.0;
    derivative = (yj - yi) / step +
                 ((-3.0 * b * b + 1.0) * y2i + (3.0 * a * a - 1.0) * y2j) *
                     step / 6.0;
  }
};

std::vector<double> readValues(std::istream &in, int count, const char *what) {
  std::vector<double> values(static_cast<size_t>(count));
  for (int i = 0; i < count; ++i) {
    if (!(in >> values[static_cast<size_t>(i)])) {
      throw std::runtime_error(std::string("setfl EAM: short ") + what);
    }
  }
  return values;
}

} // namespace

struct SetflEam::Tables {
  double cutoff{0.0};
  double dr{1.0};
  Spline embedding;
  Spline density;
  Spline zphi;
  int atomicNumber{29};
};

SetflEam::SetflEam(const Parameters &params, std::string path)
    : Potential(PotType::EAM_CU, params) {
  std::ifstream file(path);
  if (!file) {
    throw std::runtime_error("setfl EAM: cannot open " + path);
  }
  // The comment may be several lines. The grid line is the first one
  // whose first and third tokens are large integers (Nrho, Nr).
  std::string line;
  int nElements = 0;
  std::string name;
  int nrho = 0;
  int nr = 0;
  double drho = 0.0;
  double dr = 0.0;
  double cutoff = 0.0;
  bool found = false;
  std::string pending;
  while (std::getline(file, line)) {
    std::istringstream row(line);
    int first = 0;
    double second = 0.0;
    int third = 0;
    double fourth = 0.0;
    double fifth = 0.0;
    if (row >> first >> second >> third >> fourth >> fifth && first > 100 &&
        third > 100 && second > 0.0 && fourth > 0.0 && fifth > 0.0) {
      nrho = first;
      drho = second;
      nr = third;
      dr = fourth;
      cutoff = fifth;
      found = true;
      break;
    }
    pending = line;
  }
  if (!found) {
    throw std::runtime_error("setfl EAM: no grid line in " + path);
  }
  {
    std::istringstream header(pending);
    header >> nElements >> name;
  }
  if (nElements != 1) {
    throw std::runtime_error("setfl EAM: only a one-element file is read");
  }
  int atomicNumber = 0;
  double mass = 0.0;
  double alat = 0.0;
  std::string lattice;
  file >> atomicNumber >> mass >> alat >> lattice;
  if (!file || nrho < 2 || nr < 2 || !(drho > 0.0) || !(dr > 0.0) ||
      !(cutoff > 0.0)) {
    throw std::runtime_error("setfl EAM: bad header in " + path);
  }
  auto tables = std::make_shared<Tables>();
  tables->cutoff = cutoff;
  tables->dr = dr;
  tables->atomicNumber = atomicNumber;
  tables->embedding.origin = 0.0;
  tables->embedding.step = drho;
  tables->embedding.y = readValues(file, nrho, "embedding table");
  tables->embedding.build();
  tables->density.origin = 0.0;
  tables->density.step = dr;
  tables->density.y = readValues(file, nr, "density table");
  tables->density.build();
  tables->zphi.origin = 0.0;
  tables->zphi.step = dr;
  tables->zphi.y = readValues(file, nr, "pair table");
  tables->zphi.build();
  tables_ = std::move(tables);
}

void SetflEam::force(long N, const double *R, const int *, double *F, double *U,
                     double *, const double *box) {
  const auto &tables = *tables_;
  Matrix3d cell = Matrix3d::Zero();
  for (int i = 0; i < 9; ++i) {
    cell(i / 3, i % 3) = box[i];
  }
  const double volume = std::abs(cell.determinant());
  const Matrix3d inverse = cell.inverse();
  const double cutoff2 = tables.cutoff * tables.cutoff;
  std::vector<double> rho(static_cast<size_t>(N), 0.0);
  for (long i = 0; i < N; ++i) {
    F[3 * i] = 0.0;
    F[3 * i + 1] = 0.0;
    F[3 * i + 2] = 0.0;
  }
  auto minimumImage = [&](long i, long j, Eigen::Vector3d &dr) {
    Eigen::RowVector3d row(R[3 * i] - R[3 * j], R[3 * i + 1] - R[3 * j + 1],
                           R[3 * i + 2] - R[3 * j + 2]);
    Eigen::RowVector3d frac = row * inverse;
    for (int axis = 0; axis < 3; ++axis) {
      frac[axis] -= std::floor(frac[axis] + 0.5);
    }
    row = frac * cell;
    dr = row.transpose();
  };
  for (long i = 0; i < N; ++i) {
    for (long j = i + 1; j < N; ++j) {
      Eigen::Vector3d dr;
      minimumImage(i, j, dr);
      const double r2 = dr.squaredNorm();
      if (!(r2 < cutoff2) || !(r2 > 0.0)) {
        continue;
      }
      const double r = std::sqrt(r2);
      double density = 0.0;
      double unused = 0.0;
      tables.density.eval(r, density, unused);
      rho[static_cast<size_t>(i)] += density;
      rho[static_cast<size_t>(j)] += density;
    }
  }
  double energy = 0.0;
  std::vector<double> dEmbedding(static_cast<size_t>(N), 0.0);
  for (long i = 0; i < N; ++i) {
    double value = 0.0;
    double derivative = 0.0;
    tables.embedding.eval(rho[static_cast<size_t>(i)], value, derivative);
    energy += value;
    dEmbedding[static_cast<size_t>(i)] = derivative;
  }
  Matrix3d virial = Matrix3d::Zero();
  for (long i = 0; i < N; ++i) {
    for (long j = i + 1; j < N; ++j) {
      Eigen::Vector3d dr;
      minimumImage(i, j, dr);
      const double r2 = dr.squaredNorm();
      if (!(r2 < cutoff2) || !(r2 > 0.0)) {
        continue;
      }
      const double r = std::sqrt(r2);
      double density = 0.0;
      double densityDerivative = 0.0;
      tables.density.eval(r, density, densityDerivative);
      double z = 0.0;
      double zDerivative = 0.0;
      tables.zphi.eval(r, z, zDerivative);
      const double phi = z / r;
      const double phiDerivative = zDerivative / r - z / (r * r);
      energy += phi;
      const double dEdr =
          (dEmbedding[static_cast<size_t>(i)] +
           dEmbedding[static_cast<size_t>(j)]) *
              densityDerivative +
          phiDerivative;
      const Eigen::Vector3d hat = dr / r;
      const Eigen::Vector3d fi = -dEdr * hat;
      F[3 * i] += fi[0];
      F[3 * i + 1] += fi[1];
      F[3 * i + 2] += fi[2];
      F[3 * j] -= fi[0];
      F[3 * j + 1] -= fi[1];
      F[3 * j + 2] -= fi[2];
      virial += dEdr * (dr * dr.transpose()) / r;
    }
  }
  *U = energy;
  stress_.setZero();
  if (volume > 0.0) {
    stress_ = virial / volume;
  }
}

} // namespace eonc

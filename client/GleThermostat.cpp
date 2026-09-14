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
#include <fstream>
#include <sstream>
#include <vector>

#include <unsupported/Eigen/MatrixFunctions>

#include "eon/GleThermostat.h"

namespace eonc {

MatrixXd GleThermostat::loadDriftMatrix(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    return MatrixXd();
  }
  std::vector<std::vector<double>> rows;
  std::string line;
  while (std::getline(in, line)) {
    const auto hash = line.find('#');
    if (hash != std::string::npos) {
      line = line.substr(0, hash);
    }
    std::istringstream ss(line);
    std::vector<double> row;
    double v;
    while (ss >> v) {
      row.push_back(v);
    }
    if (!row.empty()) {
      rows.push_back(std::move(row));
    }
  }
  const auto n = static_cast<long>(rows.size());
  if (n == 0) {
    return MatrixXd();
  }
  MatrixXd a(n, n);
  for (long i = 0; i < n; ++i) {
    if (static_cast<long>(rows[i].size()) != n) {
      return MatrixXd();
    }
    for (long j = 0; j < n; ++j) {
      a(i, j) = rows[i][j];
    }
  }
  return a;
}

GleThermostat::GleThermostat(const MatrixXd &a_drift, double kbt,
                             double dt_half, long n_dof) {
  const long dim = a_drift.rows();
  if (dim == 0 || a_drift.cols() != dim || !(kbt > 0.0) || !(dt_half > 0.0) ||
      n_dof <= 0) {
    return;
  }
  m_T = (-a_drift * dt_half).exp();
  // Fluctuation-dissipation: S S^T = kbt (I - T T^T). The right-hand
  // side is symmetric positive semi-definite; factor through its
  // eigendecomposition with negative rounding clamped so a
  // deterministic (dt -> 0) or overdamped block cannot poison the
  // Cholesky.
  const MatrixXd cov =
      kbt * (MatrixXd::Identity(dim, dim) - m_T * m_T.transpose());
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(cov);
  if (es.info() != Eigen::Success) {
    return;
  }
  m_S = es.eigenvectors() *
        es.eigenvalues().cwiseMax(0.0).cwiseSqrt().asDiagonal();
  m_Z = MatrixXd::Zero(dim, n_dof);
  m_valid = true;
}

void GleThermostat::apply(VectorXd &vel, const VectorXd &masses3N,
                          const std::function<double()> &gauss) {
  if (!m_valid || vel.size() != m_Z.cols()) {
    return;
  }
  const VectorXd sqrtM = masses3N.cwiseSqrt();
  m_Z.row(0) = (vel.cwiseProduct(sqrtM)).transpose();
  MatrixXd xi(m_Z.rows(), m_Z.cols());
  for (long i = 0; i < xi.rows(); ++i) {
    for (long j = 0; j < xi.cols(); ++j) {
      xi(i, j) = gauss();
    }
  }
  m_Z = m_T * m_Z + m_S * xi;
  vel = m_Z.row(0).transpose().cwiseQuotient(sqrtM);
}

} // namespace eonc

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
#include "eon/SolidState.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace eonc::solid {
namespace {

void projectTranslations(VectorXd &v, int nFree) {
  if (nFree < 2 || v.size() < 3L * nFree) {
    return;
  }
  for (int pass = 0; pass < 2; ++pass) {
    Eigen::Vector3d mean = Eigen::Vector3d::Zero();
    for (int i = 0; i < nFree; ++i) {
      mean += v.segment<3>(3 * i);
    }
    mean /= static_cast<double>(nFree);
    for (int i = 0; i < nFree; ++i) {
      v.segment<3>(3 * i) -= mean;
    }
  }
}

using GradientFn = std::function<VectorXd(const VectorXd &)>;
using EnergyFn = std::function<double(const VectorXd &)>;

VectorXd projectedAction(const VectorXd &x, const VectorXd &direction,
                         const GradientFn &gradient, int nFree, double step,
                         bool project, int &evaluations) {
  VectorXd dir = direction;
  if (project) {
    projectTranslations(dir, nFree);
  }
  const double norm = dir.norm();
  if (norm < 1e-14) {
    return VectorXd::Zero(x.size());
  }
  const VectorXd unit = dir / norm;
  const VectorXd gPlus = gradient(x + step * unit);
  const VectorXd gMinus = gradient(x - step * unit);
  evaluations += 2;
  if (gPlus.size() != x.size() || gMinus.size() != x.size()) {
    throw std::invalid_argument("index-check gradient length changed");
  }
  if (!gPlus.allFinite() || !gMinus.allFinite()) {
    throw std::runtime_error("index-check gradient is not finite");
  }
  VectorXd action = (gPlus - gMinus) * (norm / (2.0 * step));
  if (project) {
    projectTranslations(action, nFree);
  }
  return action;
}

std::optional<VectorXd> freshDirection(const std::vector<VectorXd> &basis,
                                       int nFree, int n, bool project) {
  for (int seed = 0; seed < n; ++seed) {
    VectorXd w = VectorXd::Zero(n);
    w[seed] = 1.0;
    if (project) {
      projectTranslations(w, nFree);
    }
    for (const VectorXd &q : basis) {
      w -= q * w.dot(q);
    }
    const double norm = w.norm();
    if (norm > 1e-8) {
      return w / norm;
    }
  }
  return std::nullopt;
}

std::pair<std::vector<double>, MatrixXd>
ritzFromActions(const std::vector<VectorXd> &basis,
                const std::vector<VectorXd> &actions) {
  const int m = static_cast<int>(std::min(actions.size(), basis.size()));
  MatrixXd projected = MatrixXd::Zero(m, m);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      projected(i, j) = basis[static_cast<size_t>(i)].dot(
          actions[static_cast<size_t>(j)]);
    }
  }
  for (int i = 0; i < m; ++i) {
    for (int j = i + 1; j < m; ++j) {
      const double mean = 0.5 * (projected(i, j) + projected(j, i));
      projected(i, j) = mean;
      projected(j, i) = mean;
    }
  }
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> sym =
      projected;
  Eigen::SelfAdjointEigenSolver<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
      solver(sym);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("index check could not diagonalise the Ritz matrix");
  }
  const int n = basis.empty() ? 0 : static_cast<int>(basis.front().size());
  MatrixXd vectors = MatrixXd::Zero(n, m);
  for (int column = 0; column < m; ++column) {
    for (int row = 0; row < m; ++row) {
      vectors.col(column) +=
          solver.eigenvectors()(row, column) * basis[static_cast<size_t>(row)];
    }
    const double norm = vectors.col(column).norm();
    if (norm > 0.0) {
      vectors.col(column) /= norm;
    }
  }
  std::vector<double> values(static_cast<size_t>(m));
  for (int i = 0; i < m; ++i) {
    values[static_cast<size_t>(i)] = solver.eigenvalues()(i);
  }
  return {std::move(values), std::move(vectors)};
}

IndexSpectrum ritzLowest(const VectorXd &x, const GradientFn &gradient,
                         int nFree, const IndexCheck &check) {
  const int n = static_cast<int>(x.size());
  const int k = std::min(check.lowest, n);
  VectorXd seed(n);
  for (int i = 0; i < n; ++i) {
    seed[i] = std::sin((i + 1) * 0.618033988749895);
  }
  if (check.project) {
    projectTranslations(seed, nFree);
  }
  const double seedNorm = seed.norm();
  IndexSpectrum spectrum;
  if (seedNorm < 1e-14) {
    spectrum.values.assign(static_cast<size_t>(k), 0.0);
    spectrum.vectors = MatrixXd::Zero(n, k);
    return spectrum;
  }
  seed /= seedNorm;
  std::vector<VectorXd> basis{seed};
  std::vector<VectorXd> actions;
  const int cap = std::min(n, std::max(64, k * 4));
  std::optional<std::vector<double>> previous;
  std::optional<std::pair<std::vector<double>, MatrixXd>> kept;
  bool confirmed = false;
  while (static_cast<int>(basis.size()) <= cap) {
    VectorXd hq = projectedAction(x, basis.back(), gradient, nFree, check.step,
                                  check.project, spectrum.evaluations);
    if (hq.size() != n) {
      throw std::invalid_argument("index-check Hessian action length changed");
    }
    VectorXd w = hq;
    for (int pass = 0; pass < 2; ++pass) {
      for (const VectorXd &q : basis) {
        w -= q * w.dot(q);
      }
    }
    actions.push_back(std::move(hq));
    const double off = w.norm();
    const int m = static_cast<int>(actions.size());
    const bool chainDone = off < 1e-8 || m == n;
    const bool ready = chainDone || (m >= k && m % 4 == 0);
    if (ready) {
      auto [values, vectors] = ritzFromActions(basis, actions);
      bool stable = false;
      if (previous.has_value()) {
        stable = true;
        const int take = std::min(
            {k, static_cast<int>(values.size()),
             static_cast<int>(previous->size())});
        for (int i = 0; i < take; ++i) {
          const double now = values[static_cast<size_t>(i)];
          const double then = (*previous)[static_cast<size_t>(i)];
          if (std::abs(now - then) >
              1e-5 * (1.0 + std::max(std::abs(now), std::abs(then)))) {
            stable = false;
            break;
          }
        }
      }
      previous = std::vector<double>(
          values.begin(),
          values.begin() + std::min(k, static_cast<int>(values.size())));
      kept = std::make_pair(values, vectors);
      if (m < n && (chainDone || (stable && !confirmed))) {
        if (auto next = freshDirection(
                std::vector<VectorXd>(basis.begin(), basis.begin() + m), nFree,
                n, check.project)) {
          confirmed = stable;
          basis.resize(static_cast<size_t>(m));
          basis.push_back(*next);
          continue;
        }
      }
      if (stable || chainDone || m == n) {
        break;
      }
    }
    if (m == n || !(off > 0.0) || !std::isfinite(off)) {
      break;
    }
    w /= off;
    basis.resize(static_cast<size_t>(m));
    basis.push_back(w);
  }
  if (!kept.has_value()) {
    throw std::runtime_error("index check did not build a Ritz pair");
  }
  const int take = std::min(k, static_cast<int>(kept->first.size()));
  spectrum.values.assign(kept->first.begin(), kept->first.begin() + take);
  spectrum.vectors = kept->second.leftCols(take);
  spectrum.negative = 0;
  for (double value : spectrum.values) {
    if (value < check.negative) {
      ++spectrum.negative;
    }
  }
  return spectrum;
}

void requireCheck(const VectorXd &x, const IndexCheck &check) {
  if (x.size() == 0 || !x.allFinite()) {
    throw std::invalid_argument("index check needs a finite coordinate");
  }
  if (check.lowest < 1 || !(check.step > 0.0) || !std::isfinite(check.step)) {
    throw std::invalid_argument(
        "index check needs a positive difference step and at least one eigenvalue");
  }
  if (!std::isfinite(check.negative) || !(check.escape >= 0.0) ||
      !std::isfinite(check.escape)) {
    throw std::invalid_argument(
        "index check cutoff and escape step must be finite");
  }
}

double rmsToNearest(const VectorXd &point,
                    const std::vector<VectorXd> &previous) {
  double nearest = std::numeric_limits<double>::max();
  for (const VectorXd &earlier : previous) {
    const double rms =
        (point - earlier).norm() / std::sqrt(static_cast<double>(point.size()));
    nearest = std::min(nearest, rms);
  }
  return nearest;
}

double rayleigh(const VectorXd &q, const VectorXd &direction,
                const GradientFn &gradient) {
  constexpr double step = 1e-4;
  const VectorXd gPlus = gradient(q + step * direction);
  const VectorXd gMinus = gradient(q - step * direction);
  if (gPlus.size() != q.size() || gMinus.size() != q.size()) {
    throw std::invalid_argument("index-check Rayleigh gradient length changed");
  }
  return direction.dot(gPlus - gMinus) / (2.0 * step);
}

bool walkOff(VectorXd &x, VectorXd &hardMode, const EnergyFn &energy,
             const GradientFn &gradient, const IndexSpectrum &spectrum,
             const IndexCheck &check, std::vector<VectorXd> &rejectedAt) {
  if (spectrum.negative < 2 ||
      spectrum.vectors.cols() < spectrum.negative ||
      spectrum.vectors.rows() != x.size()) {
    return false;
  }
  const VectorXd start = x;
  rejectedAt.push_back(start);
  const VectorXd soft = spectrum.vectors.col(spectrum.negative - 1);
  hardMode = spectrum.vectors.col(0);
  const double reference = energy(start);
  double sign = 0.0;
  double lowest = reference;
  for (double candidate : {1.0, -1.0}) {
    const VectorXd trial = start + soft * (candidate * 0.05);
    if (!trial.allFinite()) {
      continue;
    }
    const double value = energy(trial);
    if (value < lowest) {
      lowest = value;
      sign = candidate;
    }
  }
  if (sign == 0.0) {
    return false;
  }
  const double step = std::max(check.escape, 0.01);
  VectorXd q = start;
  for (int i = 0; i < 12; ++i) {
    q += soft * (sign * step);
    const double drift = (q - start).dot(hardMode);
    q -= hardMode * drift;
    if (rayleigh(q, soft, gradient) >= check.negative) {
      break;
    }
  }
  if (rmsToNearest(q, rejectedAt) <= 1e-4) {
    return false;
  }
  x = q;
  return true;
}

} // namespace

IndexSpectrum indexAt(Matter &matter, const Chart &chart,
                      const IndexCheck &check) {
  const AtomMatrix savedPos = matter.getPositions();
  const Matrix3d savedCell = matter.getCell();
  const VectorXd x = chart.coordinates(matter);
  requireCheck(x, check);
  const int nFree = static_cast<int>(chart.freeAtoms().size());
  auto gradient = [&](const VectorXd &q) {
    chart.setCoordinates(matter, q);
    return chart.energyGradient(matter);
  };
  IndexSpectrum spectrum;
  try {
    spectrum = ritzLowest(x, gradient, nFree, check);
  } catch (...) {
    matter.setPositions(savedPos);
    matter.setCell(savedCell);
    throw;
  }
  matter.setPositions(savedPos);
  matter.setCell(savedCell);
  return spectrum;
}

IndexSpectrum indexAtoms(Matter &matter, const IndexCheck &check) {
  const VectorXd x = matter.getPositionsV();
  requireCheck(x, check);
  const AtomMatrix saved = matter.getPositions();
  // The frozen-cell coordinate is every atom. A fixed atom is still a
  // row of that vector; the dimer does not move it, and the translation
  // to remove is the one shared by the free rows. All-free systems, which
  // are the copper cell and the analytic checks, project every row.
  const int nFree = static_cast<int>(matter.numberOfAtoms());
  auto gradient = [&](const VectorXd &q) {
    matter.setPositionsV(q);
    const AtomMatrix forces = matter.getForces();
    // Evaluate before `forces` dies. A bare Eigen expression would keep a
    // pointer into that temporary and read freed memory.
    const VectorXd gradient =
        -Eigen::Map<const VectorXd>(forces.data(), forces.size());
    return gradient;
  };
  IndexSpectrum spectrum;
  try {
    spectrum = ritzLowest(x, gradient, nFree, check);
  } catch (...) {
    matter.setPositions(saved);
    throw;
  }
  matter.setPositions(saved);
  return spectrum;
}

bool leaveHigherIndex(Matter &matter, const Chart &chart, VectorXd &hardMode,
                      const IndexSpectrum &spectrum, const IndexCheck &check,
                      std::vector<VectorXd> &rejectedAt) {
  const AtomMatrix savedPos = matter.getPositions();
  const Matrix3d savedCell = matter.getCell();
  VectorXd x = chart.coordinates(matter);
  auto energy = [&](const VectorXd &q) {
    chart.setCoordinates(matter, q);
    return chart.enthalpy(matter);
  };
  auto gradient = [&](const VectorXd &q) {
    chart.setCoordinates(matter, q);
    return chart.energyGradient(matter);
  };
  const bool moved =
      walkOff(x, hardMode, energy, gradient, spectrum, check, rejectedAt);
  if (!moved) {
    matter.setPositions(savedPos);
    matter.setCell(savedCell);
    return false;
  }
  chart.setCoordinates(matter, x);
  return true;
}

bool leaveHigherIndexAtoms(Matter &matter, VectorXd &hardMode,
                           const IndexSpectrum &spectrum,
                           const IndexCheck &check,
                           std::vector<VectorXd> &rejectedAt) {
  const AtomMatrix saved = matter.getPositions();
  VectorXd x = matter.getPositionsV();
  auto energy = [&](const VectorXd &q) {
    matter.setPositionsV(q);
    return matter.getPotentialEnergy();
  };
  auto gradient = [&](const VectorXd &q) {
    matter.setPositionsV(q);
    const AtomMatrix forces = matter.getForces();
    // Evaluate before `forces` dies. A bare Eigen expression would keep a
    // pointer into that temporary and read freed memory.
    const VectorXd gradient =
        -Eigen::Map<const VectorXd>(forces.data(), forces.size());
    return gradient;
  };
  const bool moved =
      walkOff(x, hardMode, energy, gradient, spectrum, check, rejectedAt);
  if (!moved) {
    matter.setPositions(saved);
    return false;
  }
  matter.setPositionsV(x);
  return true;
}

} // namespace eonc::solid

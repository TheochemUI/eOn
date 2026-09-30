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
#include "eon/Tunneling.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace eonc::tunneling {

double massWeightedDistance(const Matter &a, const Matter &b) {
  if (a.numberOfAtoms() != b.numberOfAtoms()) {
    throw std::invalid_argument("the structures hold different atom counts");
  }
  const AtomMatrix dr = a.pbc(b.getPositions() - a.getPositions());
  const auto mass = a.getMasses();
  double sum = 0.0;
  for (long i = 0; i < a.numberOfAtoms(); ++i) {
    if (mass(i) <= 0.0) {
      throw std::invalid_argument(
          "every atom needs a positive mass for a mass-weighted path");
    }
    sum += mass(i) * dr.row(i).squaredNorm();
  }
  return std::sqrt(sum);
}

std::vector<double>
massWeightedPath(const std::vector<std::shared_ptr<Matter>> &band) {
  std::vector<double> s{0.0};
  s.reserve(band.size());
  for (size_t i = 1; i < band.size(); ++i) {
    s.push_back(s.back() + massWeightedDistance(*band[i - 1], *band[i]));
  }
  return s;
}

Profile::Profile(std::vector<double> s, std::vector<double> v)
    : s_(std::move(s)),
      v_(std::move(v)),
      m_(s_.size(), 0.0) {
  const size_t n = s_.size();
  if (n < 2 || v_.size() != n) {
    throw std::invalid_argument(
        "a profile needs matching s and V with two points or more");
  }
  for (size_t k = 1; k < n; ++k) {
    if (!(s_[k] > s_[k - 1])) {
      throw std::invalid_argument(
          "the path coordinate must increase along the band");
    }
  }
  for (size_t k = 1; k + 1 < n; ++k) {
    const double h0 = s_[k] - s_[k - 1];
    const double h1 = s_[k + 1] - s_[k];
    const double d0 = (v_[k] - v_[k - 1]) / h0;
    const double d1 = (v_[k + 1] - v_[k]) / h1;
    if (d0 * d1 <= 0.0) {
      m_[k] = 0.0;
    } else {
      const double w1 = 2.0 * h1 + h0;
      const double w2 = h1 + 2.0 * h0;
      m_[k] = (w1 + w2) / (w1 / d0 + w2 / d1);
    }
  }
}

double Profile::operator()(double x) const {
  x = std::clamp(x, s_.front(), s_.back());
  auto it = std::upper_bound(s_.begin(), s_.end(), x);
  size_t k = static_cast<size_t>(std::distance(s_.begin(), it));
  k = std::clamp<size_t>(k == 0 ? 0 : k - 1, 0, s_.size() - 2);
  const double h = s_[k + 1] - s_[k];
  const double t = (x - s_[k]) / h;
  const double t2 = t * t;
  const double t3 = t2 * t;
  return (2 * t3 - 3 * t2 + 1) * v_[k] + (t3 - 2 * t2 + t) * h * m_[k] +
         (-2 * t3 + 3 * t2) * v_[k + 1] + (t3 - t2) * h * m_[k + 1];
}

double wellCurvature(const Profile &p, bool leftEnd) {
  const auto &s = p.s();
  const auto &v = p.v();
  const size_t n = s.size();
  const double top = *std::max_element(v.begin(), v.end());
  const double floor = leftEnd ? v.front() : v.back();
  const double half = 0.5 * (top - floor);
  // Sums for the normal equations of y = a x^2 + b x^3.
  double s44 = 0, s45 = 0, s55 = 0, sy2 = 0, sy3 = 0;
  size_t used = 0;
  for (size_t j = 1; j < n; ++j) {
    const size_t i = leftEnd ? j : n - 1 - j;
    const double x = leftEnd ? s[i] - s.front() : s.back() - s[i];
    const double y = v[i] - floor;
    if (y > half) {
      break;
    }
    const double x2 = x * x;
    s44 += x2 * x2;
    s45 += x2 * x2 * x;
    s55 += x2 * x2 * x2;
    sy2 += y * x2;
    sy3 += y * x2 * x;
    ++used;
  }
  if (used == 0) {
    // The next image already stands above half the barrier: the parabola
    // through it is all the band says about this well.
    const size_t i = leftEnd ? 1 : n - 2;
    const double x = leftEnd ? s[i] - s.front() : s.back() - s[i];
    return 2.0 * (v[i] - floor) / (x * x);
  }
  if (used == 1) {
    return 2.0 * sy2 / s44;
  }
  const double det = s44 * s55 - s45 * s45;
  const double a = (sy2 * s55 - sy3 * s45) / det;
  return 2.0 * a;
}

double hbarOmega(double curvature) {
  if (!(curvature > 0.0)) {
    throw std::invalid_argument("a well needs a positive curvature");
  }
  return kHbar * std::sqrt(curvature);
}

double wkbAction(const Profile &p, double energy, int points) {
  const double a = p.s().front();
  const double b = p.s().back();
  const double h = (b - a) / (points - 1);
  double sum = 0.0;
  for (int i = 0; i < points; ++i) {
    const double gap = p(a + i * h) - energy;
    const double f = gap > 0.0 ? std::sqrt(2.0 * gap) : 0.0;
    sum += (i == 0 || i == points - 1) ? 0.5 * f : f;
  }
  return sum * h / kHbar;
}

double Splitting::tlsEnergy() const { return std::hypot(delta, delta0); }

Splitting wkbSplitting(const Profile &p, double hwReactant, double hwProduct) {
  const auto &v = p.v();
  Splitting out;
  const double top = *std::max_element(v.begin(), v.end());
  out.delta = v.back() - v.front();
  out.barrier = top - v.front();
  out.hwReactant = hwReactant;
  out.hwProduct = hwProduct;
  out.referenceEnergy =
      std::max(v.front() + 0.5 * hwReactant, v.back() + 0.5 * hwProduct);
  out.action = wkbAction(p, out.referenceEnergy);
  const double hw = std::sqrt(hwReactant * hwProduct);
  out.delta0 = hw / std::numbers::pi * std::exp(-out.action);
  out.deepWells =
      (top - v.front()) > hwReactant && (top - v.back()) > hwProduct;
  return out;
}

Splitting bandSplitting(const std::vector<std::shared_ptr<Matter>> &band,
                        double referenceEnergy) {
  std::vector<double> v;
  v.reserve(band.size());
  for (const auto &image : band) {
    v.push_back(image->getPotentialEnergy() - referenceEnergy);
  }
  const Profile p(massWeightedPath(band), std::move(v));
  return wkbSplitting(p, hbarOmega(wellCurvature(p, true)),
                      hbarOmega(wellCurvature(p, false)));
}

} // namespace eonc::tunneling

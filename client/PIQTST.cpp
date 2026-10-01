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
#include "eon/PIQTST.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace eonc::piqtst {

namespace {

/// The trapezoid map from mean forces to F: F_j = sum_i W(j, i) F'_i.
MatrixXd trapezoidWeights(const std::vector<Plane> &planes) {
  const long n = static_cast<long>(planes.size());
  MatrixXd w = MatrixXd::Zero(n, n);
  for (long j = 1; j < n; ++j) {
    w.row(j) = w.row(j - 1);
    const double h =
        planes[static_cast<size_t>(j)].s - planes[static_cast<size_t>(j - 1)].s;
    w(j, j - 1) += 0.5 * h;
    w(j, j) += 0.5 * h;
  }
  return w;
}

VectorXd forceErrors(const std::vector<Plane> &planes) {
  VectorXd e(static_cast<long>(planes.size()));
  for (size_t i = 0; i < planes.size(); ++i) {
    e(static_cast<long>(i)) = planes[i].meanForceError;
  }
  return e;
}

double propagated(const VectorXd &gradient, const VectorXd &errors) {
  return std::sqrt(gradient.cwiseProduct(errors).squaredNorm());
}

} // namespace

void integrate(std::vector<Plane> &planes) {
  if (planes.empty()) {
    return;
  }
  const MatrixXd w = trapezoidWeights(planes);
  VectorXd f(static_cast<long>(planes.size()));
  for (size_t i = 0; i < planes.size(); ++i) {
    f(static_cast<long>(i)) = planes[i].meanForce;
  }
  const VectorXd errors = forceErrors(planes);
  const VectorXd values = w * f;
  for (size_t j = 0; j < planes.size(); ++j) {
    planes[j].freeEnergy = values(static_cast<long>(j));
    planes[j].freeEnergyError =
        propagated(w.row(static_cast<long>(j)).transpose(), errors);
  }
}

std::vector<Plane> scan(Potential &pot, const Coordinate &c,
                        const ScanOptions &o) {
  const long dof = 3 * c.atoms;
  if (c.atoms < 1 || static_cast<long>(c.masses.size()) != c.atoms ||
      static_cast<long>(c.free.size()) != dof || c.reference.size() != dof ||
      c.direction.size() != dof) {
    throw std::invalid_argument("piqtst: coordinate sizes do not match");
  }
  if (o.planes.size() < 2) {
    throw std::invalid_argument("piqtst: at least two planes are needed");
  }
  for (size_t j = 1; j < o.planes.size(); ++j) {
    if (!(o.planes[j] > o.planes[j - 1])) {
      throw std::invalid_argument("piqtst: plane positions must ascend");
    }
  }
  if (o.equilibration < 0 || o.production < 2 || o.blocks < 2 ||
      o.production < o.blocks) {
    throw std::invalid_argument(
        "piqtst: sampling needs at least as many steps as blocks, and two "
        "blocks");
  }
  // Cartesian normal a = M^(1/2) n and step b = M^(-1/2) n. With n a unit
  // vector, a . b = 1, so the plane through reference + s b holds s, and
  // dF/ds = -<a . f> / |a|^2 = -<a_hat . f> / |a|.
  VectorXd a = VectorXd::Zero(dof);
  VectorXd b = VectorXd::Zero(dof);
  double nn = 0.0;
  for (long i = 0; i < dof; ++i) {
    if (!c.free[static_cast<size_t>(i)]) {
      continue;
    }
    const double sm = std::sqrt(c.masses[static_cast<size_t>(i / 3)]);
    a(i) = sm * c.direction(i);
    b(i) = c.direction(i) / sm;
    nn += c.direction(i) * c.direction(i);
  }
  if (std::abs(nn - 1.0) > 1e-8) {
    throw std::invalid_argument(
        "piqtst: the direction must be a unit vector on the free coordinates");
  }
  const double aNorm = a.norm();
  auto seedAt = [&](double s) -> VectorXd {
    if (o.seed) {
      VectorXd x = o.seed(s);
      if (x.size() != dof) {
        throw std::invalid_argument("piqtst: a seed has the wrong length");
      }
      return x;
    }
    return c.reference + s * b;
  };

  pathintegral::RingPolymer ring(c.atoms, c.masses, c.numbers, c.free, o.ring);
  const long beads = o.ring.beads;
  const long blockSize = o.production / o.blocks;
  std::vector<Plane> out;
  out.reserve(o.planes.size());
  for (size_t j = 0; j < o.planes.size(); ++j) {
    const double s = o.planes[j];
    const VectorXd origin = c.reference + s * b;
    const VectorXd target = seedAt(s);
    if (j == 0) {
      ring.setAllBeads(target.data());
    } else {
      const VectorXd shift = target - ring.centroid();
      std::vector<VectorXd> moved = ring.beads();
      for (auto &q : moved) {
        q += shift;
      }
      ring.setBeads(moved);
    }
    ring.setHyperplane(a, origin);

    const long batches0 = ring.batches();
    for (long step = 0; step < o.equilibration; ++step) {
      ring.step(pot, c.box, false);
    }
    Plane plane;
    plane.s = s;
    plane.centroid = VectorXd::Zero(dof);
    VectorXd spread2 = VectorXd::Zero(dof);
    double sum = 0.0;
    double blockSum = 0.0;
    std::vector<double> blockMeans;
    for (long step = 0; step < o.production; ++step) {
      ring.resetAverages();
      ring.step(pot, c.box, true);
      const double fn = ring.meanForce();
      sum += fn;
      blockSum += fn;
      if ((step + 1) % blockSize == 0 &&
          static_cast<long>(blockMeans.size()) < o.blocks) {
        blockMeans.push_back(blockSum / static_cast<double>(blockSize));
        blockSum = 0.0;
      }
      const VectorXd centroid = ring.centroid();
      plane.centroid += centroid;
      for (const auto &q : ring.beads()) {
        spread2 += (q - centroid).cwiseAbs2();
      }
    }
    const double steps = static_cast<double>(o.production);
    plane.centroid /= steps;
    spread2 /= steps * static_cast<double>(beads);
    plane.spread.resize(static_cast<size_t>(dof));
    for (long i = 0; i < dof; ++i) {
      plane.spread[static_cast<size_t>(i)] = std::sqrt(spread2(i));
    }
    const double mean = sum / steps;
    double blockMean = 0.0;
    for (const double m : blockMeans) {
      blockMean += m;
    }
    const double nb = static_cast<double>(blockMeans.size());
    blockMean /= nb;
    double var = 0.0;
    for (const double m : blockMeans) {
      var += (m - blockMean) * (m - blockMean);
    }
    var /= nb * (nb - 1.0);
    plane.meanForce = -mean / aNorm;
    plane.meanForceError = std::sqrt(var) / aNorm;
    plane.batches = ring.batches() - batches0;
    out.push_back(std::move(plane));
  }
  integrate(out);
  return out;
}

Rate rate(const std::vector<Plane> &planes, double beta) {
  const long n = static_cast<long>(planes.size());
  if (n < 2) {
    throw std::invalid_argument("piqtst: the rate needs two planes");
  }
  if (!(beta > 0.0)) {
    throw std::invalid_argument("piqtst: the rate needs a positive beta");
  }
  const MatrixXd w = trapezoidWeights(planes);
  const VectorXd errors = forceErrors(planes);
  Rate r;
  double fMin = std::numeric_limits<double>::infinity();
  for (long j = 0; j < n; ++j) {
    const double f = planes[static_cast<size_t>(j)].freeEnergy;
    if (f < fMin) {
      fMin = f;
      r.reactant = j;
    }
  }
  const long top = n - 1;
  const double fTop = planes[static_cast<size_t>(top)].freeEnergy;
  r.barrier = fTop - fMin;
  r.barrierError =
      propagated((w.row(top) - w.row(r.reactant)).transpose(), errors);
  r.firstPlaneHeight = beta * (planes.front().freeEnergy - fMin);

  // Z = int exp(-beta F) ds by the trapezoid rule, referenced to fMin.
  VectorXd weight(n);
  double z = 0.0;
  for (long j = 0; j < n; ++j) {
    const double left = j > 0 ? planes[static_cast<size_t>(j)].s -
                                    planes[static_cast<size_t>(j - 1)].s
                              : 0.0;
    const double right = j + 1 < n ? planes[static_cast<size_t>(j + 1)].s -
                                         planes[static_cast<size_t>(j)].s
                                   : 0.0;
    weight(j) =
        0.5 * (left + right) *
        std::exp(-beta * (planes[static_cast<size_t>(j)].freeEnergy - fMin));
    z += weight(j);
  }
  weight /= z;
  r.logRate = std::log(0.5 * std::sqrt(2.0 / (std::numbers::pi * beta))) -
              beta * (fTop - fMin) - std::log(z);
  const VectorXd gradient =
      -beta * w.row(top).transpose() + beta * (w.transpose() * weight);
  r.logRateError = propagated(gradient, errors);
  return r;
}

} // namespace eonc::piqtst

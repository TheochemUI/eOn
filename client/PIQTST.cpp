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

/// Cartesian normal a = M^(1/2) n and step b = M^(-1/2) n. With n a unit
/// vector, a . b = 1, so the plane through reference + s b holds s, and
/// dF/ds = -<a . f> / |a|^2 = -<a_hat . f> / |a|.
struct Axes {
  VectorXd a;
  VectorXd b;
};

Axes axes(const Coordinate &c) {
  const long dof = 3 * c.atoms;
  if (c.atoms < 1 || static_cast<long>(c.masses.size()) != c.atoms ||
      static_cast<long>(c.free.size()) != dof || c.reference.size() != dof ||
      c.direction.size() != dof) {
    throw std::invalid_argument("piqtst: coordinate sizes do not match");
  }
  Axes out{VectorXd::Zero(dof), VectorXd::Zero(dof)};
  double nn = 0.0;
  for (long i = 0; i < dof; ++i) {
    if (!c.free[static_cast<size_t>(i)]) {
      continue;
    }
    const double sm = std::sqrt(c.masses[static_cast<size_t>(i / 3)]);
    out.a(i) = sm * c.direction(i);
    out.b(i) = c.direction(i) / sm;
    nn += c.direction(i) * c.direction(i);
  }
  if (std::abs(nn - 1.0) > 1e-8) {
    throw std::invalid_argument(
        "piqtst: the direction must be a unit vector on the free coordinates");
  }
  return out;
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
  const Axes ax = axes(c);
  const VectorXd &a = ax.a;
  const VectorXd &b = ax.b;
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
  // The reactant is the lowest F before the dividing plane; s* itself is
  // the top, and a minimum there is a plane that is no barrier, which the
  // barrier's sign then shows.
  const long top = n - 1;
  double fMin = std::numeric_limits<double>::infinity();
  for (long j = 0; j < top; ++j) {
    const double f = planes[static_cast<size_t>(j)].freeEnergy;
    if (f < fMin) {
      fMin = f;
      r.reactant = j;
    }
  }
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

Recrossing recrossing(Potential &pot, const Coordinate &c,
                      const RecrossingOptions &o) {
  const long dof = 3 * c.atoms;
  const Axes ax = axes(c);
  if (o.parents < 2 || o.children < 1 || o.spacing < 1 || o.steps < 4 ||
      o.equilibration < 0) {
    throw std::invalid_argument(
        "piqtst: recrossing needs two parents, one child, a positive "
        "spacing and four steps");
  }
  const VectorXd origin = c.reference + o.s * ax.b;
  VectorXd start = origin;
  if (o.seed) {
    start = o.seed(o.s);
    if (start.size() != dof) {
      throw std::invalid_argument("piqtst: a seed has the wrong length");
    }
  }
  // Bennett-Chandler needs parents from the ring polymer's own constrained
  // Boltzmann distribution. PIGLET's coloured noise samples another one,
  // tuned for bead averages, so the parents always take PILE.
  pathintegral::Options parentOptions = o.ring;
  parentOptions.thermostat = pathintegral::Thermostat::Pile;
  parentOptions.gleFile.clear();
  pathintegral::RingPolymer parent(c.atoms, c.masses, c.numbers, c.free,
                                   parentOptions);
  parent.setAllBeads(start.data());
  parent.setHyperplane(ax.a, origin);
  for (long step = 0; step < o.equilibration; ++step) {
    parent.step(pot, c.box, false);
  }

  // The children draw every normal mode from the free-ring Boltzmann
  // distribution at beta / P, which PILE's initial momenta are, and carry
  // their own stream so parent and child noise stay independent.
  pathintegral::Options childOptions = o.ring;
  childOptions.thermostat = pathintegral::Thermostat::Pile;
  childOptions.gleFile.clear();
  childOptions.seed = o.ring.seed + 0x9E3779B97F4A7C15ULL;
  pathintegral::RingPolymer child(c.atoms, c.masses, c.numbers, c.free,
                                  childOptions);

  const long n = o.steps + 1;
  // Per parent: sum over children of sdot(0) h(s(t) - s*) at each time,
  // and of sdot(0) h(sdot(0)).
  std::vector<VectorXd> numerator(static_cast<size_t>(o.parents),
                                  VectorXd::Zero(n));
  std::vector<double> denominator(static_cast<size_t>(o.parents), 0.0);
  Recrossing out;
  for (long p = 0; p < o.parents; ++p) {
    for (long step = 0; step < o.spacing; ++step) {
      parent.step(pot, c.box, false);
    }
    const std::vector<VectorXd> beads = parent.beads();
    VectorXd &num = numerator[static_cast<size_t>(p)];
    double &den = denominator[static_cast<size_t>(p)];
    for (long k = 0; k < o.children; ++k) {
      child.setBeads(beads);
      child.thermalMomenta();
      std::vector<VectorXd> reversed = child.momenta();
      const double forward = ax.a.dot(child.centroidVelocity());
      for (int sign = 0; sign < 2; ++sign) {
        if (sign == 1) {
          for (auto &v : reversed) {
            v = -v;
          }
          child.setBeads(beads);
          child.setMomenta(reversed);
        }
        const double sdot = sign == 0 ? forward : -forward;
        const double flux = sdot > 0.0 ? sdot : 0.0;
        den += flux;
        num(0) += flux;
        for (long i = 1; i < n; ++i) {
          child.nveStep(pot, c.box);
          const double s = ax.a.dot(child.centroid() - c.reference);
          if (s > o.s) {
            num(i) += sdot;
          }
        }
        ++out.trajectories;
      }
    }
  }

  VectorXd total = VectorXd::Zero(n);
  double totalDen = 0.0;
  for (long p = 0; p < o.parents; ++p) {
    total += numerator[static_cast<size_t>(p)];
    totalDen += denominator[static_cast<size_t>(p)];
  }
  if (!(totalDen > 0.0)) {
    throw std::runtime_error("piqtst: no child left the plane forward");
  }
  out.time.resize(static_cast<size_t>(n));
  out.kappa.resize(static_cast<size_t>(n));
  for (long i = 0; i < n; ++i) {
    out.time[static_cast<size_t>(i)] = static_cast<double>(i) * o.ring.dt;
    out.kappa[static_cast<size_t>(i)] = total(i) / totalDen;
  }

  // Plateau numerators per parent, averaged over the last quarter, and the
  // leave-one-parent-out ratios.
  const long first = n - 1 - o.steps / 4;
  const double span = static_cast<double>(n - first);
  std::vector<double> plateauNum(static_cast<size_t>(o.parents), 0.0);
  double sumNum = 0.0;
  for (long p = 0; p < o.parents; ++p) {
    plateauNum[static_cast<size_t>(p)] =
        numerator[static_cast<size_t>(p)].tail(n - first).sum() / span;
    sumNum += plateauNum[static_cast<size_t>(p)];
  }
  out.plateau = sumNum / totalDen;
  const double np = static_cast<double>(o.parents);
  std::vector<double> leaveOut(static_cast<size_t>(o.parents), 0.0);
  double meanLeave = 0.0;
  for (long p = 0; p < o.parents; ++p) {
    leaveOut[static_cast<size_t>(p)] =
        (sumNum - plateauNum[static_cast<size_t>(p)]) /
        (totalDen - denominator[static_cast<size_t>(p)]);
    meanLeave += leaveOut[static_cast<size_t>(p)];
  }
  meanLeave /= np;
  double var = 0.0;
  for (const double v : leaveOut) {
    var += (v - meanLeave) * (v - meanLeave);
  }
  out.plateauError = std::sqrt((np - 1.0) / np * var);
  out.batches = parent.batches() + child.batches();
  return out;
}

} // namespace eonc::piqtst

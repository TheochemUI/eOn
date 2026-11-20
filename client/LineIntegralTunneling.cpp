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
#include "eon/LineIntegralTunneling.h"

#include "eon/NEBForceProjection.h"
#include "eon/ObjectiveFunction.h"
#include "eon/Optimizer.h"
#include "magic_enum/magic_enum.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>

namespace eonc::tunneling {

void validateOptions(const line_integral_tunneling_options_t &o) {
  auto require = [](bool valid, const char *key, const char *requirement) {
    if (!valid) {
      throw std::invalid_argument(std::string("[Line Integral Tunneling] ") +
                                  key + " " + requirement);
    }
  };
  require(o.images >= 3, "images", "must be at least 3");
  require(o.replicas >= 4, "replicas", "must be at least 4");
  require(std::isfinite(o.energy_offset) && o.energy_offset > 0.0,
          "energy_offset", "must be finite and positive");
  require(std::isfinite(o.lowest_temperature) && o.lowest_temperature > 0.0,
          "lowest_temperature",
          "must be finite and positive");
  require(std::isfinite(o.energy_step) && o.energy_step > 0.0, "energy_step",
          "must be finite and positive");
  require(std::isfinite(o.log_rate_tolerance) && o.log_rate_tolerance > 0.0,
          "log_rate_tolerance", "must be finite and positive");
  require(std::isfinite(o.spring) && o.spring > 0.0, "spring",
          "must be finite and positive");
  require(std::isfinite(o.contour) && o.contour > 0.0, "contour",
          "must be finite and positive");
  require(std::isfinite(o.force_tolerance) && o.force_tolerance > 0.0,
          "force_tolerance", "must be finite and positive");
  require(o.max_iterations >= 1, "max_iterations", "must be at least 1");
  require(std::isfinite(o.li_neb_time_step) && o.li_neb_time_step > 0.0,
          "li_neb_time_step", "must be finite and positive");
  require(std::isfinite(o.period_tolerance) && o.period_tolerance > 0.0,
          "period_tolerance", "must be finite and positive");
  require(std::isfinite(o.determinant_estimated_error) &&
              o.determinant_estimated_error > 0.0,
          "determinant_estimated_error", "must be finite and positive");
  if (o.optimizer != "velocity_verlet") {
    const auto method =
        magic_enum::enum_cast<OptType>(o.optimizer,
                                      magic_enum::case_insensitive);
    require(method && *method != OptType::None && *method != OptType::Unknown,
            "optimizer",
            "must be velocity_verlet or a recognized eOn optimizer token");
  }
}

namespace {

void requireImageCount(long images) {
  if (images < 3) {
    throw std::invalid_argument(
        "line-integral tunneling needs at least three images");
  }
}

void requirePath(const std::vector<VectorXd> &path) {
  requireImageCount(static_cast<long>(path.size()));
  const long freedom = path.front().size();
  if (freedom < 1) {
    throw std::invalid_argument("line-integral path has no coordinates");
  }
  for (const auto &q : path) {
    if (q.size() != freedom || !q.allFinite()) {
      throw std::invalid_argument(
          "line-integral path coordinates are inconsistent");
    }
  }
}

std::vector<double> chordCoordinate(const std::vector<VectorXd> &path) {
  requirePath(path);
  std::vector<double> s(path.size(), 0.0);
  for (size_t i = 1; i < path.size(); ++i) {
    const double length = (path[i] - path[i - 1]).norm();
    if (!(length > 0.0)) {
      throw std::invalid_argument(
          "line-integral path contains adjacent duplicate images");
    }
    s[i] = s[i - 1] + length;
  }
  return s;
}

class NaturalCubic {
public:
  NaturalCubic(std::vector<double> x, std::vector<double> y)
      : x_(std::move(x)),
        y_(std::move(y)),
        second_(y_.size(), 0.0) {
    if (x_.size() < 2 || x_.size() != y_.size()) {
      throw std::invalid_argument("cubic spline needs matching samples");
    }
    const size_t n = x_.size();
    for (size_t i = 1; i < n; ++i) {
      if (!(x_[i] > x_[i - 1])) {
        throw std::invalid_argument(
            "cubic spline coordinates must be strictly increasing");
      }
    }
    if (n == 2) {
      return;
    }
    MatrixXd system =
        MatrixXd::Zero(static_cast<long>(n - 2), static_cast<long>(n - 2));
    VectorXd rhs(static_cast<long>(n - 2));
    for (size_t i = 1; i + 1 < n; ++i) {
      const double left = x_[i] - x_[i - 1];
      const double right = x_[i + 1] - x_[i];
      const long row = static_cast<long>(i - 1);
      system(row, row) = 2.0 * (left + right);
      if (i > 1) {
        system(row, row - 1) = left;
      }
      if (i + 2 < n) {
        system(row, row + 1) = right;
      }
      rhs(row) =
          6.0 * ((y_[i + 1] - y_[i]) / right - (y_[i] - y_[i - 1]) / left);
    }
    const VectorXd interior = system.ldlt().solve(rhs);
    if (!interior.allFinite()) {
      throw std::runtime_error("cubic spline system is singular");
    }
    for (size_t i = 1; i + 1 < n; ++i) {
      second_[i] = interior(static_cast<long>(i - 1));
    }
  }

  double value(double x) const {
    const size_t i = interval(x);
    const double h = x_[i + 1] - x_[i];
    const double a = (x_[i + 1] - x) / h;
    const double b = (x - x_[i]) / h;
    return a * y_[i] + b * y_[i + 1] +
           ((a * a * a - a) * second_[i] + (b * b * b - b) * second_[i + 1]) *
               h * h / 6.0;
  }

  double derivative(double x) const {
    const size_t i = interval(x);
    const double h = x_[i + 1] - x_[i];
    const double a = (x_[i + 1] - x) / h;
    const double b = (x - x_[i]) / h;
    return (y_[i + 1] - y_[i]) / h + ((-3.0 * a * a + 1.0) * second_[i] +
                                      (3.0 * b * b - 1.0) * second_[i + 1]) *
                                         h / 6.0;
  }

private:
  size_t interval(double x) const {
    if (!(x >= x_.front() && x <= x_.back())) {
      throw std::out_of_range("cubic spline query is outside its path");
    }
    if (x >= x_.back()) {
      return x_.size() - 2;
    }
    const auto upper = std::upper_bound(x_.begin(), x_.end(), x);
    return static_cast<size_t>(upper - x_.begin() - 1);
  }

  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> second_;
};

class CubicVectorPath {
public:
  explicit CubicVectorPath(const std::vector<VectorXd> &path)
      : coordinate_(chordCoordinate(path)) {
    const long freedom = path.front().size();
    splines_.reserve(static_cast<size_t>(freedom));
    for (long a = 0; a < freedom; ++a) {
      std::vector<double> values(path.size());
      for (size_t i = 0; i < path.size(); ++i) {
        values[i] = path[i](a);
      }
      splines_.emplace_back(coordinate_, std::move(values));
    }
  }

  VectorXd value(double s) const {
    VectorXd q(static_cast<long>(splines_.size()));
    for (long a = 0; a < q.size(); ++a) {
      q(a) = splines_[static_cast<size_t>(a)].value(s);
    }
    return q;
  }

  VectorXd derivative(double s) const {
    VectorXd q(static_cast<long>(splines_.size()));
    for (long a = 0; a < q.size(); ++a) {
      q(a) = splines_[static_cast<size_t>(a)].derivative(s);
    }
    return q;
  }

  double length() const { return coordinate_.back(); }
  const std::vector<double> &coordinate() const { return coordinate_; }

private:
  std::vector<double> coordinate_;
  std::vector<NaturalCubic> splines_;
};

class HermiteEnergy {
public:
  HermiteEnergy(std::vector<double> x, std::vector<double> value,
                std::vector<double> slope)
      : x_(std::move(x)),
        value_(std::move(value)),
        slope_(std::move(slope)) {
    if (x_.size() < 2 || x_.size() != value_.size() ||
        x_.size() != slope_.size()) {
      throw std::invalid_argument("Hermite path needs matching samples");
    }
  }

  double value(double x) const {
    const auto [i, u] = interval(x);
    const double h = x_[i + 1] - x_[i];
    const double u2 = u * u;
    const double u3 = u2 * u;
    return (2.0 * u3 - 3.0 * u2 + 1.0) * value_[i] +
           (u3 - 2.0 * u2 + u) * h * slope_[i] +
           (-2.0 * u3 + 3.0 * u2) * value_[i + 1] +
           (u3 - u2) * h * slope_[i + 1];
  }

  double derivative(double x) const {
    const auto [i, u] = interval(x);
    const double h = x_[i + 1] - x_[i];
    const double u2 = u * u;
    return ((6.0 * u2 - 6.0 * u) * value_[i] +
            (3.0 * u2 - 4.0 * u + 1.0) * h * slope_[i] +
            (-6.0 * u2 + 6.0 * u) * value_[i + 1] +
            (3.0 * u2 - 2.0 * u) * h * slope_[i + 1]) /
           h;
  }

private:
  std::pair<size_t, double> interval(double x) const {
    if (!(x >= x_.front() && x <= x_.back())) {
      throw std::out_of_range("Hermite query is outside its path");
    }
    size_t i = x_.size() - 2;
    if (x < x_.back()) {
      i = static_cast<size_t>(std::upper_bound(x_.begin(), x_.end(), x) -
                              x_.begin() - 1);
    }
    return {i, (x - x_[i]) / (x_[i + 1] - x_[i])};
  }

  std::vector<double> x_;
  std::vector<double> value_;
  std::vector<double> slope_;
};

HermiteEnergy energySpline(const LineIntegralPath &path,
                           const CubicVectorPath &coordinates) {
  if (path.energies.size() != path.images.size() ||
      path.gradients.size() != path.images.size()) {
    throw std::invalid_argument(
        "line-integral path is missing energies or gradients");
  }
  std::vector<double> values = path.energies;
  values.front() = path.energy;
  values.back() = path.energy;
  std::vector<double> slopes(path.images.size());
  const auto &s = coordinates.coordinate();
  for (size_t i = 0; i < path.images.size(); ++i) {
    slopes[i] = path.gradients[i].dot(coordinates.derivative(s[i]));
  }
  return HermiteEnergy(s, std::move(values), std::move(slopes));
}

std::vector<VectorXd> resample(const std::vector<VectorXd> &path, long images) {
  requireImageCount(images);
  const CubicVectorPath spline(path);
  std::vector<VectorXd> result;
  result.reserve(static_cast<size_t>(images));
  for (long i = 0; i < images; ++i) {
    const double fraction =
        static_cast<double>(i) / static_cast<double>(images - 1);
    result.push_back(spline.value(fraction * spline.length()));
  }
  return result;
}

VectorXd crossing(const VectorXd &lower, double lowerEnergy,
                  const VectorXd &upper, double upperEnergy, double energy) {
  if (!(lowerEnergy < energy) || !(upperEnergy >= energy)) {
    throw std::invalid_argument("MEP points do not bracket the path energy");
  }
  const double fraction = (energy - lowerEnergy) / (upperEnergy - lowerEnergy);
  return lower + fraction * (upper - lower);
}

struct BracketedPath {
  VectorXd left;
  VectorXd right;
  size_t leftLower;
  size_t rightLower;
};

BracketedPath brackets(const std::vector<VectorXd> &mep,
                       const std::vector<double> &energies, double energy) {
  if (mep.size() < 3 || mep.size() != energies.size()) {
    throw std::invalid_argument("MEP coordinates and energies do not match");
  }
  requirePath(mep);
  const size_t saddle = static_cast<size_t>(
      std::max_element(energies.begin(), energies.end()) - energies.begin());
  if (saddle == 0 || saddle + 1 == mep.size() || !(energies[saddle] > energy)) {
    throw std::invalid_argument(
        "MEP has no interior saddle above the requested energy");
  }
  size_t left = saddle;
  while (left > 0 && !(energies[left - 1] < energy)) {
    --left;
  }
  if (left == 0) {
    throw std::invalid_argument(
        "MEP reactant side does not cross the requested energy");
  }
  const size_t leftLower = left - 1;
  size_t right = saddle;
  while (right + 1 < mep.size() && !(energies[right + 1] < energy)) {
    ++right;
  }
  if (right + 1 == mep.size()) {
    throw std::invalid_argument(
        "MEP product side does not cross the requested energy");
  }
  for (size_t i = 0; i < left; ++i) {
    if (energies[i] >= energy) {
      throw std::invalid_argument(
          "MEP has more than one disjoint barrier above the requested energy");
    }
  }
  for (size_t i = right + 1; i < mep.size(); ++i) {
    if (energies[i] >= energy) {
      throw std::invalid_argument(
          "MEP has more than one disjoint barrier above the requested energy");
    }
  }
  const size_t rightLower = right + 1;
  return {crossing(mep[leftLower], energies[leftLower], mep[left],
                   energies[left], energy),
          crossing(mep[rightLower], energies[rightLower], mep[right],
                   energies[right], energy),
          leftLower, rightLower};
}

void appendDistinct(std::vector<VectorXd> &path, const VectorXd &point) {
  if (path.empty() || (point - path.back()).norm() > 0.0) {
    path.push_back(point);
  }
}

class LineIntegralEvaluation {
public:
  LineIntegralEvaluation(std::vector<VectorXd> images, double energy,
                         const BatchPotential &potential,
                         const LineIntegralOptions &options,
                         PathAlignment align)
      : images_(std::move(images)),
        energy_(energy),
        potential_(potential),
        options_(options),
        align_(std::move(align)),
        gradient_(static_cast<long>(images_.size()) * images_.front().size()) {
    requirePath(images_);
    if (!std::isfinite(energy_) || !std::isfinite(options_.spring) ||
        !std::isfinite(options_.contour) ||
        !std::isfinite(options_.forceTolerance) ||
        !std::isfinite(options_.timeStep) || !(options_.spring > 0.0) ||
        !(options_.contour > 0.0) || !(options_.forceTolerance > 0.0) ||
        !(options_.timeStep > 0.0) ||
        options_.maxIterations < 1) {
      throw std::invalid_argument(
          "line-integral energy and options must be finite and positive");
    }
    evaluate();
  }

  long freedom() const { return images_.front().size(); }
  long degrees() const { return static_cast<long>(images_.size()) * freedom(); }
  const VectorXd &gradient() const { return gradient_; }
  const std::vector<VectorXd> &images() const { return images_; }
  double objective() const { return objective_; }
  double interiorGradient() const { return interiorGradient_; }
  double endpointGradient() const { return endpointGradient_; }
  bool converged() const {
    return interiorGradient_ < options_.forceTolerance &&
           endpointGradient_ < options_.forceTolerance;
  }

  VectorXd positions() const {
    VectorXd flat(degrees());
    for (long i = 0; i < static_cast<long>(images_.size()); ++i) {
      flat.segment(i * freedom(), freedom()) = images_[static_cast<size_t>(i)];
    }
    return flat;
  }

  void positions(const VectorXd &flat) {
    if (flat.size() != degrees()) {
      throw std::invalid_argument("LI-NEB optimizer changed its dimension");
    }
    for (long i = 0; i < static_cast<long>(images_.size()); ++i) {
      images_[static_cast<size_t>(i)] = flat.segment(i * freedom(), freedom());
    }
    if (align_) {
      align_(images_);
    }
    evaluate();
  }

  LineIntegralPath result(long iterations) const {
    LineIntegralPath out;
    out.images = images_;
    out.energies = energies_;
    out.gradients = physicalGradients_;
    out.energy = energy_;
    out.action = action();
    out.interiorGradient = interiorGradient_;
    out.endpointGradient = endpointGradient_;
    out.iterations = iterations;
    out.converged = converged();
    return out;
  }

private:
  void evaluate() {
    potential_(images_, energies_, physicalGradients_);
    if (energies_.size() != images_.size() ||
        physicalGradients_.size() != images_.size()) {
      throw std::runtime_error(
          "LI-NEB potential returned the wrong number of images");
    }
    for (size_t i = 0; i < images_.size(); ++i) {
      if (physicalGradients_[i].size() != freedom() ||
          !physicalGradients_[i].allFinite() || !std::isfinite(energies_[i])) {
        throw std::runtime_error("LI-NEB potential returned invalid data");
      }
    }

    const size_t n = images_.size();
    std::vector<double> momentum(n, 0.0);
    for (size_t i = 1; i + 1 < n; ++i) {
      if (!(energies_[i] > energy_)) {
        throw std::runtime_error(
            "an interior LI-NEB image is not above the path energy");
      }
      momentum[i] = std::sqrt(2.0 * (energies_[i] - energy_));
    }
    std::vector<VectorXd> opt(n, VectorXd::Zero(freedom()));
    objective_ = 0.0;
    for (size_t i = 0; i + 1 < n; ++i) {
      const double length = (images_[i + 1] - images_[i]).norm();
      if (!(length > 0.0)) {
        throw std::runtime_error("LI-NEB images have collided");
      }
      objective_ += 0.5 * (momentum[i] + momentum[i + 1]) * length;
    }
    interiorGradient_ = 0.0;
    for (size_t i = 1; i + 1 < n; ++i) {
      const VectorXd previous = images_[i] - images_[i - 1];
      const VectorXd next = images_[i + 1] - images_[i];
      const double previousLength = previous.norm();
      const double nextLength = next.norm();
      if (!(previousLength > 0.0) || !(nextLength > 0.0)) {
        throw std::runtime_error("LI-NEB images have collided");
      }
      const VectorXd previousUnit = previous / previousLength;
      const VectorXd nextUnit = next / nextLength;
      const VectorXd tangent =
          neb::computeTangent(next, previous, energies_[i], energies_[i - 1],
                              energies_[i + 1], false);
      const VectorXd actionGradient =
          0.5 *
          ((momentum[i] + momentum[i - 1]) * previousUnit -
           (momentum[i] + momentum[i + 1]) * nextUnit +
           physicalGradients_[i] * (previousLength + nextLength) / momentum[i]);
      const VectorXd perpendicular = neb::forcePerp(actionGradient, tangent);
      const VectorXd spring =
          options_.spring * (nextLength - previousLength) * tangent;
      opt[i] = perpendicular - spring;
      const double norm = perpendicular.norm();
      if (!std::isfinite(norm)) {
        throw std::runtime_error("LI-NEB interior gradient is not finite");
      }
      interiorGradient_ = std::max(interiorGradient_, norm);
    }

    endpointGradient_ = 0.0;
    for (const size_t i : {size_t{0}, n - 1}) {
      const size_t neighbour = i == 0 ? 1 : n - 2;
      const VectorXd spring =
          options_.spring * (images_[neighbour] - images_[i]);
      const VectorXd force = -physicalGradients_[i];
      const double forceNorm = force.norm();
      if (!(forceNorm > 0.0)) {
        throw std::runtime_error(
            "LI-NEB endpoint has no energy-contour normal");
      }
      const VectorXd normal = force / forceNorm;
      const VectorXd endpointForce =
          spring -
          (spring.dot(normal) - options_.contour * (energies_[i] - energy_)) *
              normal;
      opt[i] = -endpointForce;
      const double norm = endpointForce.norm();
      if (!std::isfinite(norm)) {
        throw std::runtime_error("LI-NEB endpoint gradient is not finite");
      }
      endpointGradient_ = std::max(endpointGradient_, norm);
    }

    for (long i = 0; i < static_cast<long>(n); ++i) {
      gradient_.segment(i * freedom(), freedom()) = opt[static_cast<size_t>(i)];
    }
  }

  double action() const {
    const CubicVectorPath coordinates(images_);
    LineIntegralPath current;
    current.images = images_;
    current.energies = energies_;
    current.gradients = physicalGradients_;
    current.energy = energy_;
    const HermiteEnergy energyPath = energySpline(current, coordinates);
    static constexpr double nodes[] = {-0.9602898564975363, -0.7966664774136267,
                                       -0.5255324099163290, -0.1834346424956498,
                                       0.1834346424956498,  0.5255324099163290,
                                       0.7966664774136267,  0.9602898564975363};
    static constexpr double weights[] = {
        0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
        0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763};
    const auto &s = coordinates.coordinate();
    double result = 0.0;
    for (size_t i = 0; i + 1 < s.size(); ++i) {
      const double midpoint = 0.5 * (s[i] + s[i + 1]);
      const double half = 0.5 * (s[i + 1] - s[i]);
      for (size_t q = 0; q < std::size(nodes); ++q) {
        const double sample = midpoint + half * nodes[q];
        const double excess = energyPath.value(sample) - energy_;
        if (!(excess >= 0.0)) {
          throw std::runtime_error(
              "cubic LI-NEB path falls below its turning-point energy");
        }
        result += half * weights[q] * std::sqrt(2.0 * excess) *
                  coordinates.derivative(sample).norm();
      }
    }
    return result;
  }

  std::vector<VectorXd> images_;
  double energy_;
  const BatchPotential &potential_;
  LineIntegralOptions options_;
  PathAlignment align_;
  std::vector<double> energies_;
  std::vector<VectorXd> physicalGradients_;
  VectorXd gradient_;
  double objective_{0.0};
  double interiorGradient_{0.0};
  double endpointGradient_{0.0};
};

class LineIntegralObjective final : public ObjectiveFunction {
public:
  LineIntegralObjective(const Parameters &parameters,
                        std::shared_ptr<LineIntegralEvaluation> evaluation)
      : ObjectiveFunction(parameters),
        evaluation_(std::move(evaluation)) {}

  double getEnergy() override { return evaluation_->objective(); }
  VectorXd getGradient(bool) override { return evaluation_->gradient(); }
  void setPositions(const VectorXd &x) override { evaluation_->positions(x); }
  VectorXd getPositions() override { return evaluation_->positions(); }
  int degreesOfFreedom() override {
    return static_cast<int>(evaluation_->degrees());
  }
  bool isConverged() override { return evaluation_->converged(); }
  double getConvergence() override {
    return std::max(evaluation_->interiorGradient(),
                    evaluation_->endpointGradient());
  }
  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    return a - b;
  }
  bool supportsFiniteDifferenceCurvature() const override { return false; }

private:
  std::shared_ptr<LineIntegralEvaluation> evaluation_;
};

struct TimeInterval {
  double left;
  double right;
  double integral;
  double error;
};

template <class Function>
double gaussLegendreEight(const Function &function, double left, double right) {
  constexpr double nodes[] = {0.18343464249564980494, 0.52553240991632898582,
                              0.79666647741362673959, 0.96028985649753623168};
  constexpr double weights[] = {0.36268378337836198297, 0.31370664587788728734,
                                0.22238103445337447054, 0.10122853629037625915};
  const double midpoint = 0.5 * (left + right);
  const double radius = 0.5 * (right - left);
  double sum = 0.0;
  for (size_t i = 0; i < std::size(nodes); ++i) {
    const double offset = radius * nodes[i];
    sum += weights[i] *
           (function(midpoint - offset) + function(midpoint + offset));
  }
  return radius * sum;
}

template <class Function>
TimeInterval timeInterval(const Function &function, double left, double right) {
  const double midpoint = 0.5 * (left + right);
  if (!(midpoint > left && midpoint < right)) {
    throw std::runtime_error(
        "imaginary-time quadrature exhausted coordinate precision");
  }
  const double whole = gaussLegendreEight(function, left, right);
  const double split = gaussLegendreEight(function, left, midpoint) +
                       gaussLegendreEight(function, midpoint, right);
  constexpr double refinementRatio = static_cast<double>((size_t{1} << 16) - 1);
  return {left, right, split, std::abs(split - whole) / refinementRatio};
}

template <class Function>
std::vector<TimeInterval> timeQuadrature(const Function &function,
                                         const std::vector<double> &boundaries,
                                         double tolerance) {
  if (boundaries.size() < 2) {
    throw std::invalid_argument(
        "imaginary-time quadrature needs interval boundaries");
  }
  std::vector<TimeInterval> intervals;
  intervals.reserve(boundaries.size() - 1);
  for (size_t i = 1; i < boundaries.size(); ++i) {
    intervals.push_back(
        timeInterval(function, boundaries[i - 1], boundaries[i]));
  }
  auto totalError = [&]() {
    double error = 0.0;
    for (const TimeInterval &interval : intervals) {
      error += interval.error;
    }
    return error;
  };
  while (totalError() > tolerance) {
    const auto largest =
        std::max_element(intervals.begin(), intervals.end(),
                         [](const TimeInterval &a, const TimeInterval &b) {
                           return a.error < b.error;
                         });
    const double midpoint = 0.5 * (largest->left + largest->right);
    const TimeInterval left = timeInterval(function, largest->left, midpoint);
    const TimeInterval right = timeInterval(function, midpoint, largest->right);
    *largest = left;
    intervals.push_back(right);
  }
  std::sort(intervals.begin(), intervals.end(),
            [](const TimeInterval &a, const TimeInterval &b) {
              return a.left < b.left;
            });
  return intervals;
}

template <class Function>
double invertTime(const Function &function, const TimeInterval &interval,
                  double timeFromLeft, double tolerance) {
  if (std::abs(timeFromLeft) <= tolerance) {
    return interval.left;
  }
  if (std::abs(interval.integral - timeFromLeft) <= tolerance) {
    return interval.right;
  }
  double left = interval.left;
  double right = interval.right;
  while (true) {
    const double midpoint = 0.5 * (left + right);
    const double time = gaussLegendreEight(function, interval.left, midpoint);
    if (std::abs(time - timeFromLeft) <= tolerance) {
      return midpoint;
    }
    if (!(midpoint > left && midpoint < right)) {
      throw std::runtime_error(
          "imaginary-time inversion did not meet tolerance");
    }
    if (time < timeFromLeft) {
      left = midpoint;
    } else {
      right = midpoint;
    }
  }
}

struct PeriodEstimate {
  double period;
  double quadratureError;
};

PeriodEstimate periodEstimate(const LineIntegralPath &path, double tolerance) {
  const CubicVectorPath coordinates(path.images);
  const HermiteEnergy energy = energySpline(path, coordinates);
  const auto integrand = [&](double theta) {
    const double sine = std::sin(theta);
    const double s = coordinates.length() * sine * sine;
    const double excess = energy.value(s) - path.energy;
    if (!(excess > 0.0)) {
      throw std::runtime_error(
          "cubic path energy does not remain above the turning contour");
    }
    const double jacobian = coordinates.length() * std::sin(2.0 * theta);
    return jacobian * coordinates.derivative(s).norm() /
           std::sqrt(2.0 * excess);
  };
  std::vector<double> thetaBoundaries;
  thetaBoundaries.reserve(coordinates.coordinate().size());
  for (double s : coordinates.coordinate()) {
    thetaBoundaries.push_back(std::asin(std::sqrt(s / coordinates.length())));
  }
  const std::vector<TimeInterval> quadrature =
      timeQuadrature(integrand, thetaBoundaries, 0.5 * tolerance);
  double halfPeriod = 0.0;
  double halfError = 0.0;
  for (const TimeInterval &interval : quadrature) {
    halfPeriod += interval.integral;
    halfError += interval.error;
  }
  return {2.0 * halfPeriod, 2.0 * halfError};
}

LineIntegralPath coarsenedPeriodPath(const LineIntegralPath &path) {
  if (path.images.size() < 5) {
    throw std::invalid_argument(
        "period discretization error needs at least five LI images");
  }
  LineIntegralPath coarse;
  coarse.energy = path.energy;
  for (size_t i = 0; i < path.images.size(); i += 2) {
    coarse.images.push_back(path.images[i]);
    coarse.energies.push_back(path.energies[i]);
    coarse.gradients.push_back(path.gradients[i]);
  }
  if ((path.images.size() - 1) % 2 != 0) {
    coarse.images.push_back(path.images.back());
    coarse.energies.push_back(path.energies.back());
    coarse.gradients.push_back(path.gradients.back());
  }
  return coarse;
}

LineIntegralPath secantPeriodPath(const LineIntegralPath &path) {
  LineIntegralPath secant = path;
  const CubicVectorPath coordinates(path.images);
  const std::vector<double> &s = coordinates.coordinate();
  for (size_t i = 0; i < path.images.size(); ++i) {
    double slope = 0.0;
    if (i == 0) {
      slope = (path.energies[1] - path.energies[0]) / (s[1] - s[0]);
    } else if (i + 1 == path.images.size()) {
      slope = (path.energies[i] - path.energies[i - 1]) / (s[i] - s[i - 1]);
    } else {
      slope =
          (path.energies[i + 1] - path.energies[i - 1]) / (s[i + 1] - s[i - 1]);
    }
    const VectorXd tangent = coordinates.derivative(s[i]);
    if (!(tangent.squaredNorm() > 0.0)) {
      throw std::runtime_error(
          "period discretization estimate has no path tangent");
    }
    secant.gradients[i] = slope * tangent / tangent.squaredNorm();
  }
  return secant;
}

double endpointPeriodError(const LineIntegralPath &path,
                           const CubicVectorPath &coordinates) {
  double halfPeriodError = 0.0;
  const size_t endpoints[] = {0, path.images.size() - 1};
  for (size_t endpoint : endpoints) {
    const double s = endpoint == 0 ? 0.0 : coordinates.coordinate().back();
    const VectorXd derivative = coordinates.derivative(s);
    const double energyError = std::abs(path.energies[endpoint] - path.energy);
    if (energyError > 0.0) {
      const double slope = std::abs(path.gradients[endpoint].dot(derivative));
      if (!(slope > 0.0)) {
        throw std::runtime_error(
            "turning-contour error has no endpoint energy slope");
      }
      halfPeriodError +=
          derivative.norm() * std::sqrt(2.0 * energyError) / slope;
    }
  }
  return 2.0 * halfPeriodError;
}

} // namespace

std::vector<VectorXd> lineIntegralSeed(const std::vector<VectorXd> &mep,
                                       const std::vector<double> &energies,
                                       double energy, long images) {
  const BracketedPath b = brackets(mep, energies, energy);
  std::vector<VectorXd> exposed;
  appendDistinct(exposed, b.left);
  for (size_t i = b.leftLower + 1; i < b.rightLower; ++i) {
    appendDistinct(exposed, mep[i]);
  }
  appendDistinct(exposed, b.right);
  return resample(exposed, images);
}

std::vector<VectorXd> lineIntegralContinuation(
    const std::vector<VectorXd> &mep, const std::vector<double> &energies,
    const LineIntegralPath &previous, double energy, long images) {
  requirePath(previous.images);
  if (!(energy < previous.energy)) {
    throw std::invalid_argument("LI-NEB continuation energy must decrease");
  }
  const BracketedPath b = brackets(mep, energies, energy);
  const BracketedPath old = brackets(mep, energies, previous.energy);
  if (b.leftLower > old.leftLower || b.rightLower < old.rightLower) {
    throw std::runtime_error(
        "lower-energy MEP crossings do not enclose the previous path");
  }
  std::vector<VectorXd> exposed;
  appendDistinct(exposed, b.left);
  for (size_t i = b.leftLower + 1; i <= old.leftLower; ++i) {
    appendDistinct(exposed, mep[i]);
  }
  for (const auto &q : previous.images) {
    appendDistinct(exposed, q);
  }
  for (size_t i = old.rightLower; i < b.rightLower; ++i) {
    appendDistinct(exposed, mep[i]);
  }
  appendDistinct(exposed, b.right);
  return resample(exposed, images);
}

std::vector<VectorXd>
resampleLineIntegralPath(const std::vector<VectorXd> &path, long images) {
  return resample(path, images);
}

double
nextLineIntegralEnergy(const std::vector<LineIntegralContinuationPoint> &points,
                       double maximumEnergyStep, double logRateTolerance,
                       double lowestTemperature) {
  if (points.empty() || !(maximumEnergyStep > 0.0) ||
      !(logRateTolerance > 0.0) || !(lowestTemperature > 0.0)) {
    throw std::invalid_argument(
        "energy continuation needs a point and positive controls");
  }
  const LineIntegralContinuationPoint &current = points.back();
  if (!(current.temperature > lowestTemperature)) {
    throw std::invalid_argument(
        "energy continuation has already reached the lowest temperature");
  }
  double decrement = maximumEnergyStep;
  if (points.size() > 1) {
    const LineIntegralContinuationPoint &previous = points[points.size() - 2];
    const double energyDifference = previous.energy - current.energy;
    const double temperatureDifference =
        previous.temperature - current.temperature;
    const double logRateDifference =
        std::abs(current.logRate - previous.logRate);
    if (!(energyDifference > 0.0) || !(temperatureDifference > 0.0) ||
        !(logRateDifference > 0.0)) {
      throw std::runtime_error(
          "energy continuation requires decreasing energy and temperature "
          "and a changing log rate");
    }
    const double rateControlled =
        logRateTolerance * energyDifference / logRateDifference;
    decrement = std::min(maximumEnergyStep, rateControlled);
  }
  if (!(decrement > 0.0) || !std::isfinite(decrement)) {
    throw std::runtime_error(
        "energy continuation produced an invalid energy step");
  }
  return current.energy - decrement;
}

LineIntegralPath optimizeLineIntegralPath(std::vector<VectorXd> guess,
                                          double energy,
                                          const BatchPotential &potential,
                                          const LineIntegralOptions &options,
                                          const PathAlignment &align) {
  auto evaluation = std::make_shared<LineIntegralEvaluation>(
      std::move(guess), energy, potential, options, align);
  VectorXd velocity = VectorXd::Zero(evaluation->degrees());
  long iteration = 0;
  while (!evaluation->converged() && iteration < options.maxIterations) {
    const VectorXd acceleration = -evaluation->gradient();
    VectorXd positions =
        evaluation->positions() + options.timeStep * velocity +
        0.5 * options.timeStep * options.timeStep * acceleration;
    evaluation->positions(positions);
    const VectorXd nextAcceleration = -evaluation->gradient();
    velocity += 0.5 * options.timeStep * (acceleration + nextAcceleration);
    const double along = velocity.dot(nextAcceleration);
    if (along > 0.0) {
      velocity = along / nextAcceleration.squaredNorm() * nextAcceleration;
    } else {
      velocity.setZero();
    }
    ++iteration;
  }
  return evaluation->result(iteration);
}

LineIntegralPath optimizeLineIntegralPath(
    std::vector<VectorXd> guess, double energy, const BatchPotential &potential,
    const LineIntegralOptions &options, const Parameters &parameters,
    OptType optimizer, const PathAlignment &align) {
  auto evaluation = std::make_shared<LineIntegralEvaluation>(
      std::move(guess), energy, potential, options, align);
  auto objective =
      std::make_shared<LineIntegralObjective>(parameters, evaluation);
  auto algorithm = helpers::create::mkOptim(objective, optimizer, parameters);
  long iteration = 0;
  while (!evaluation->converged() && iteration < options.maxIterations) {
    algorithm->step(parameters.optimizer_options().max_move);
    ++iteration;
  }
  return evaluation->result(iteration);
}

LineIntegralOrbit lineIntegralOrbit(const LineIntegralPath &path, long replicas,
                                    double periodTolerance,
                                    double crossoverPeriod) {
  requirePath(path.images);
  if (!path.converged) {
    throw std::invalid_argument(
        "cannot propagate an unconverged line-integral path");
  }
  if (!std::isfinite(path.energy) || !std::isfinite(path.action) ||
      path.action < 0.0) {
    throw std::invalid_argument(
        "line-integral path energy and action must be finite and action "
        "non-negative");
  }
  if (replicas < 4 || !std::isfinite(periodTolerance) ||
      !std::isfinite(crossoverPeriod) || !(periodTolerance > 0.0) ||
      !(crossoverPeriod > 0.0)) {
    throw std::invalid_argument(
        "line-integral orbit needs at least four replicas and positive "
        "period controls");
  }
  const CubicVectorPath coordinates(path.images);
  const HermiteEnergy energy = energySpline(path, coordinates);
  const auto integrand = [&](double theta) {
    const double sine = std::sin(theta);
    const double s = coordinates.length() * sine * sine;
    const double excess = energy.value(s) - path.energy;
    if (!(excess > 0.0)) {
      throw std::runtime_error(
          "cubic path energy does not remain above the turning contour");
    }
    const double jacobian = coordinates.length() * std::sin(2.0 * theta);
    return jacobian * coordinates.derivative(s).norm() /
           std::sqrt(2.0 * excess);
  };
  std::vector<double> thetaBoundaries;
  thetaBoundaries.reserve(coordinates.coordinate().size());
  for (double s : coordinates.coordinate()) {
    thetaBoundaries.push_back(std::asin(std::sqrt(s / coordinates.length())));
  }
  const std::vector<TimeInterval> quadrature =
      timeQuadrature(integrand, thetaBoundaries, 0.5 * periodTolerance);
  double halfPeriod = 0.0;
  double halfPeriodError = 0.0;
  for (const TimeInterval &interval : quadrature) {
    halfPeriod += interval.integral;
    halfPeriodError += interval.error;
  }
  const double period = 2.0 * halfPeriod;
  const PeriodEstimate coarse =
      periodEstimate(coarsenedPeriodPath(path), periodTolerance);
  const PeriodEstimate secant =
      periodEstimate(secantPeriodPath(path), periodTolerance);
  const double discretizationError = std::max(std::abs(period - coarse.period),
                                              std::abs(period - secant.period));
  const double periodError = 2.0 * halfPeriodError + discretizationError +
                             endpointPeriodError(path, coordinates);
  if (!(period - periodError > crossoverPeriod)) {
    throw CrossoverTemperatureError(
        "line-integral orbit is not resolved below the crossover "
        "temperature");
  }
  LineIntegralOrbit orbit;
  orbit.period = period;
  orbit.periodError = periodError;
  orbit.temperature = kHbar / (kBoltzmann * period);
  orbit.energy = path.energy;
  orbit.action = path.action;
  orbit.replicas.reserve(static_cast<size_t>(replicas));
  orbit.energies.reserve(static_cast<size_t>(replicas));
  orbit.coordinate.reserve(static_cast<size_t>(replicas));
  for (long j = 0; j < replicas; ++j) {
    const double sample =
        period * static_cast<double>(j) / static_cast<double>(replicas);
    const double halfTime = sample <= halfPeriod ? sample : period - sample;
    double theta = 0.0;
    if (halfTime >= halfPeriod) {
      theta = std::numbers::pi / 2.0;
    } else if (halfTime > 0.0) {
      double elapsed = 0.0;
      for (const TimeInterval &interval : quadrature) {
        if (elapsed + interval.integral >= halfTime) {
          const double inversionTolerance =
              interval.error +
              periodTolerance / (2.0 * static_cast<double>(replicas));
          theta = invertTime(integrand, interval, halfTime - elapsed,
                             inversionTolerance);
          break;
        }
        elapsed += interval.integral;
      }
    }
    const double sine = std::sin(theta);
    const double s = coordinates.length() * sine * sine;
    orbit.coordinate.push_back(s);
    orbit.replicas.push_back(coordinates.value(s));
    orbit.energies.push_back(energy.value(s));
  }
  return orbit;
}

RateInstanton hqtstRate(const LineIntegralOrbit &orbit,
                        const HqtstProducts &products, double reactantEnergy,
                        const MatrixXd &reactantRigid,
                        const HqtstOptions &options,
                        const RingRigidBodies &rigidBodies) {
  const long replicas = static_cast<long>(orbit.replicas.size());
  if (replicas < 4 || orbit.energies.size() != orbit.replicas.size() ||
      !std::isfinite(orbit.energy) || !std::isfinite(orbit.action) ||
      orbit.action < 0.0 || !std::isfinite(orbit.period) ||
      !(orbit.period > 0.0) || !std::isfinite(reactantEnergy) ||
      !std::isfinite(options.estimatedLogRateErrorTarget) ||
      !(options.estimatedLogRateErrorTarget > 0.0) ||
      !products.ring || !products.reactant) {
    throw std::invalid_argument("HQTST orbit data are inconsistent");
  }
  for (size_t j = 0; j < orbit.replicas.size(); ++j) {
    if (!orbit.replicas[j].allFinite() ||
        !std::isfinite(orbit.energies[j])) {
      throw std::invalid_argument("HQTST orbit data are not finite");
    }
  }
  const long freedom = orbit.replicas.front().size();
  if (reactantRigid.rows() != freedom || reactantRigid.cols() >= freedom) {
    throw std::invalid_argument(
        "HQTST reactant rigid directions do not match the coordinates");
  }
  RateInstanton result;
  result.beads = orbit.replicas;
  result.energies = orbit.energies;
  result.beta = orbit.period / kHbar;
  result.betaN = result.beta / static_cast<double>(replicas);
  result.temperature = orbit.temperature;
  for (long j = 0; j < replicas; ++j) {
    const VectorXd displacement =
        result.beads[static_cast<size_t>((j + 1) % replicas)] -
        result.beads[static_cast<size_t>(j)];
    result.bN += displacement.squaredNorm();
  }
  if (!std::isfinite(result.bN) || !(result.bN > 0.0)) {
    throw std::runtime_error("HQTST orbit norm is not finite and positive");
  }
  // The closed orbit traverses the LI-NEB path twice, so its Euclidean
  // action is S(E) = 2 W(E) + E tau, and beta_N U_N = S(E) / hbar.
  result.ringPotential = static_cast<double>(replicas) *
                         (orbit.energy + 2.0 * orbit.action / orbit.period);

  const double bnh = result.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);
  std::vector<double> frequencies(static_cast<size_t>(replicas));
  for (long k = 0; k < replicas; ++k) {
    const double s = std::sin(std::numbers::pi * static_cast<double>(k) /
                              static_cast<double>(replicas));
    frequencies[static_cast<size_t>(k)] = 4.0 * c * s * s;
  }
  const long rigidModes = reactantRigid.cols();
  const StochasticLogDet well =
      wellLogDet(products.reactant, freedom, reactantRigid, frequencies,
                 options.estimatedLogRateErrorTarget, options.seed);

  std::vector<std::vector<VectorXd>> rigid;
  if (!rigidBodies.sqrtMasses.empty()) {
    rigid = ringRigidMotions(result.beads, rigidBodies);
    if (static_cast<long>(rigid.size()) != rigidModes) {
      throw std::runtime_error(
          "HQTST ring has " + std::to_string(rigid.size()) +
          " rigid motions and the reactant " + std::to_string(rigidModes));
    }
  } else {
    for (long r = 0; r < rigidModes; ++r) {
      rigid.emplace_back(
          static_cast<size_t>(replicas),
          (reactantRigid.col(r) / std::sqrt(static_cast<double>(replicas)))
              .eval());
    }
  }
  std::vector<VectorXd> cycle(static_cast<size_t>(replicas));
  for (long j = 0; j < replicas; ++j) {
    cycle[static_cast<size_t>(j)] =
        0.5 *
        (result.beads[static_cast<size_t>((j + 1) % replicas)] -
         result.beads[static_cast<size_t>((j + replicas - 1) % replicas)]);
  }
  // The kinetic preconditioner's shift is the reactant's mean curvature,
  // the scale the bead Hessians share away from the barrier.
  if (!(well.meanCurvature > 0.0)) {
    throw std::runtime_error("HQTST reactant has no positive mean curvature");
  }
  const StochasticLogDet ring =
      ringLogDetPrime(products.ring, c, well.meanCurvature, cycle, rigid,
                      options.estimatedLogRateErrorTarget, options.seed + 1);
  result.zeroEigenvalue = ring.zeroEigenvalue;
  result.negativeEigenvalue = std::min(0.0, ring.lowestEigenvalue);
  result.negativeModes = ring.negativeModes;

  const long nDrop = 1 + rigidModes;
  const double logProd =
      static_cast<double>(replicas * freedom - nDrop) * std::log(bnh) +
      0.5 * ring.value;
  result.logRateTimesZr =
      -std::log(bnh) +
      0.5 * std::log(result.bN /
                     (2.0 * std::numbers::pi * result.betaN * kHbar * kHbar)) -
      logProd - result.betaN * result.ringPotential;
  // The rigid modes leave the centroid (k = 0) factor, as the ring's own
  // rigid modes leave its product; both carry them as free particles for
  // k > 0.
  double rigidFree = 0.0;
  for (long k = 1; k < replicas; ++k) {
    rigidFree += std::log(frequencies[static_cast<size_t>(k)]);
  }
  result.logZr =
      -result.beta * reactantEnergy -
      static_cast<double>(replicas * freedom - rigidModes) * std::log(bnh) -
      0.5 * (well.value + static_cast<double>(rigidModes) * rigidFree);
  result.logRate = result.logRateTimesZr - result.logZr;
  if (!std::isfinite(result.logRate)) {
    throw std::runtime_error("HQTST log rate is not finite");
  }
  result.rate = std::exp(result.logRate) / kTimeUnitSeconds;
  if (!std::isfinite(result.rate)) {
    throw std::runtime_error("HQTST rate is not finite");
  }
  result.effectiveBarrier =
      -std::log(2.0 * std::numbers::pi * kHbar * result.beta) / result.beta -
      result.logRate / result.beta;
  result.logRateError = 0.5 * (ring.error() + well.error());
  result.converged =
      result.negativeModes == 1 &&
      result.logRateError <= options.estimatedLogRateErrorTarget;
  result.hessianProducts = ring.products * replicas + well.products;
  return result;
}

} // namespace eonc::tunneling

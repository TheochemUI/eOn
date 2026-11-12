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

#include "eon/HelperFunctions.h"
#include "eon/JobResult.h"
#include "eon/LowestEigenmode.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/Optimizer.h"
#include "eon/SafeMath.h"
#include "eon/SolidStateNEB.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <sstream>
#include <stdexcept>

namespace eonc::solid {
namespace {

void zeroUpper(Matrix3d &cell) {
  cell(0, 1) = 0.0;
  cell(0, 2) = 0.0;
  cell(1, 2) = 0.0;
}

double wrapFrac(double value) {
  double wrapped = value - std::floor(value);
  if (wrapped >= 0.5) {
    wrapped -= 1.0;
  }
  return wrapped;
}

Matrix3d symmetricPart(const Matrix3d &tensor) {
  return 0.5 * (tensor + tensor.transpose());
}

int activeCount(const Options &options) {
  return static_cast<int>(
      std::count(options.mask.begin(), options.mask.end(), true));
}

void requireChart(const Chart &chart) {
  if (chart.degreesOfFreedom() < 1) {
    throw std::invalid_argument(
        "solid-state chart has no free atoms and no free strain");
  }
}

double atanRatioOrHalfPi(double num, double denom) {
  if (!std::isfinite(num) || !std::isfinite(denom)) {
    return 0.0;
  }
  if (std::abs(denom) < eonc::safemath::eps) {
    if (std::abs(num) < eonc::safemath::eps) {
      return 0.0;
    }
    return std::copysign(eonc::helpers::pi * 0.5, num);
  }
  return std::atan(num / denom);
}

class VariableCellObjective final : public ObjectiveFunction {
public:
  VariableCellObjective(Matter &matter, Chart chart, const Parameters &params)
      : ObjectiveFunction(params), matter_{matter}, chart_{std::move(chart)} {
    eonc::helpers::requireKnownConvergenceMetric(
        params.optimizer_options().convergence_metric, "[Solid State]");
  }

  double getEnergy() override { return matter_.getPotentialEnergy(); }

  VectorXd getGradient(bool) override { return chart_.energyGradient(matter_); }

  void setPositions(const VectorXd &x) override {
    chart_.setCoordinates(matter_, x);
  }

  VectorXd getPositions() override { return chart_.coordinates(matter_); }

  int degreesOfFreedom() override { return chart_.degreesOfFreedom(); }

  bool isConverged() override {
    return getConvergence() < params.optimizer_options().converged_force;
  }

  double getConvergence() override { return convergence(); }

  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    return chart_.difference(a, b);
  }

  bool getPeriodic() const override { return true; }

private:
  Matter &matter_;
  Chart chart_;

  double convergence() {
    const VectorXd force = chart_.generalizedForce(matter_);
    const std::string &metric = params.optimizer_options().convergence_metric;
    if (metric == "norm") {
      return force.norm();
    }
    if (metric == "rms") {
      const auto n = force.size();
      return n > 0 ? force.norm() / std::sqrt(static_cast<double>(n)) : 0.0;
    }
    if (metric == "max_atom") {
      double peak = 0.0;
      const AtomMatrix cartesian = matter_.getForces();
      for (long atom : chart_.freeAtoms()) {
        peak = std::max(peak, cartesian.row(atom).norm());
      }
      return peak;
    }
    if (metric == "max_component") {
      return force.size() > 0 ? force.cwiseAbs().maxCoeff() : 0.0;
    }
    throw std::invalid_argument(
        std::format("[Solid State] unknown convergence_metric: {}", metric));
  }
};

class ClimbObjective final : public ObjectiveFunction {
public:
  ClimbObjective(Matter &matter, const Chart &chart, const Parameters &params)
      : ObjectiveFunction(params), matter_{matter}, chart_{chart} {}

  void setMode(const VectorXd &mode, double curvature) {
    mode_ = mode;
    curvature_ = curvature;
  }

  bool gradientMatchesEnergy() const override { return false; }

  double getEnergy() override { return matter_.getPotentialEnergy(); }

  VectorXd getGradient(bool) override {
    VectorXd force = chart_.generalizedForce(matter_);
    const double proj = force.dot(mode_);
    if (curvature_ > 0.0) {
      force = -proj * mode_;
    } else {
      force -= 2.0 * proj * mode_;
    }
    return -force;
  }

  void setPositions(const VectorXd &x) override {
    chart_.setCoordinates(matter_, x);
  }

  VectorXd getPositions() override { return chart_.coordinates(matter_); }

  int degreesOfFreedom() override { return chart_.degreesOfFreedom(); }

  bool isConverged() override {
    const double tol = params.saddle_search_options().converged_force > 0.0
                           ? params.saddle_search_options().converged_force
                           : params.optimizer_options().converged_force;
    return chart_.generalizedForce(matter_).norm() < tol;
  }

  double getConvergence() override {
    return chart_.generalizedForce(matter_).norm();
  }

  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    return chart_.difference(a, b);
  }

  bool getPeriodic() const override { return true; }

private:
  Matter &matter_;
  const Chart &chart_;
  VectorXd mode_;
  double curvature_{0.0};
};

ModeEstimate dimerMode(Matter &matter, const Chart &chart, VectorXd direction,
                       const Parameters &params) {
  const double fd = params.main_options().finiteDifference;
  if (!(fd > 0.0)) {
    throw std::invalid_argument(
        "solid-state dimer needs a positive finiteDifference");
  }
  eonc::safemath::safe_normalize_inplace(direction);
  const double angle = params.dimer_options().rotation_angle;
  VectorXd plane = VectorXd::Zero(direction.size());
  VectorXd rotationalOld = VectorXd::Zero(direction.size());
  VectorXd planeOld = plane;
  double lengthOld = 0.0;
  double curvature = 0.0;
  long rotations = 0;
  bool done = false;
  while (!done) {
    const VectorXd x = chart.coordinates(matter);
    const AtomMatrix savedPos = matter.getPositions();
    const Matrix3d savedCell = matter.getCell();
    const VectorXd forceCenter = chart.generalizedForce(matter);
    chart.setCoordinates(matter, x + fd * direction);
    const VectorXd forcePlus = chart.generalizedForce(matter);
    matter.setPositions(savedPos);
    matter.setCell(savedCell);

    const VectorXd forceMinus = 2.0 * forceCenter - forcePlus;
    const double projPlus = direction.dot(forcePlus);
    const double projMinus = direction.dot(forceMinus);
    VectorXd rotational =
        ((forcePlus - projPlus * direction) -
         (forceMinus - projMinus * direction)) /
        (2.0 * fd);
    curvature = (projMinus - projPlus) / (2.0 * fd);

    double gamma = 0.0;
    const double overlap = std::abs(rotational.dot(rotationalOld));
    const double oldNorm2 = rotationalOld.squaredNorm();
    if (overlap < 0.5 * oldNorm2 && oldNorm2 > 0.0) {
      gamma = rotational.dot(rotational - rotationalOld) / oldNorm2;
    }
    plane = rotational + planeOld * lengthOld * gamma;
    lengthOld = plane.norm();
    plane -= direction * plane.dot(direction);
    eonc::safemath::safe_normalize_inplace(plane);
    rotationalOld = rotational;

    const double torque = rotational.norm();
    const auto &dimer = params.dimer_options();
    if ((torque > dimer.torque_max && rotations >= dimer.rotations_max) ||
        (torque < dimer.torque_max && torque >= dimer.torque_min &&
         rotations >= dimer.rotations_min) ||
        (torque < dimer.torque_min)) {
      done = true;
    }
    if (!done) {
      const double rot1 = rotational.dot(plane);
      auto rotate = [&](double theta) {
        const double c = std::cos(theta);
        const double s = std::sin(theta);
        const VectorXd next = direction * c + plane * s;
        plane = plane * c - direction * s;
        direction = next;
        eonc::safemath::safe_normalize_inplace(direction);
      };
      rotate(angle);
      const VectorXd x2 = chart.coordinates(matter);
      const AtomMatrix pos2 = matter.getPositions();
      const Matrix3d cell2 = matter.getCell();
      const VectorXd center2 = chart.generalizedForce(matter);
      chart.setCoordinates(matter, x2 + fd * direction);
      const VectorXd plus2 = chart.generalizedForce(matter);
      matter.setPositions(pos2);
      matter.setCell(cell2);
      const VectorXd minus2 = 2.0 * center2 - plus2;
      const double pPlus = direction.dot(plus2);
      const double pMinus = direction.dot(minus2);
      rotational = ((plus2 - pPlus * direction) - (minus2 - pMinus * direction)) /
                   (2.0 * fd);
      curvature = (pMinus - pPlus) / (2.0 * fd);
      const double rot2 = rotational.dot(plane);
      const double rotChange = (rot1 - rot2) / angle;
      const double forceDimer = 0.5 * (rot1 + rot2);
      double rotationAngle =
          atanRatioOrHalfPi(2.0 * forceDimer, rotChange) / 2.0 - angle / 2.0;
      if (rotChange < 0.0) {
        rotationAngle += eonc::helpers::pi / 2.0;
      }
      rotate(rotationAngle);
      planeOld = plane;
      ++rotations;
    }
  }
  eonc::safemath::safe_normalize_inplace(direction);
  return {direction, curvature};
}

ModeEstimate lanczosMode(Matter &matter, const Chart &chart, VectorXd direction,
                         const Parameters &params) {
  const int size = static_cast<int>(direction.size());
  const long maxIters =
      std::max(1L, std::min(params.lanczos_options().max_iterations,
                            static_cast<long>(size)));
  const double dr = params.main_options().finiteDifference;
  if (!(dr > 0.0)) {
    throw std::invalid_argument(
        "solid-state Lanczos needs a positive finiteDifference");
  }
  const VectorXd origin = chart.coordinates(matter);
  const AtomMatrix savedPos = matter.getPositions();
  const Matrix3d savedCell = matter.getCell();
  const VectorXd force0 = chart.generalizedForce(matter);
  auto product = [&](const VectorXd &v) {
    chart.setCoordinates(matter, origin + dr * v);
    const VectorXd force = chart.generalizedForce(matter);
    matter.setPositions(savedPos);
    matter.setCell(savedCell);
    return -(force - force0) / dr;
  };

  MatrixXd basis(size, maxIters);
  MatrixXd tridiagonal = MatrixXd::Zero(maxIters, maxIters);
  double beta = direction.norm();
  if (beta < eonc::safemath::eps) {
    direction.setZero();
    if (size > 0) {
      direction[0] = 1.0;
    }
    beta = 1.0;
  }
  double previous = 0.0;
  VectorXd estimate = direction / beta;
  double curvature = 0.0;
  for (long i = 0; i < maxIters; ++i) {
    basis.col(i) = direction / beta;
    VectorXd u = product(basis.col(i));
    if (i > 0) {
      u -= beta * basis.col(i - 1);
    }
    const double alpha = basis.col(i).dot(u);
    direction = u - alpha * basis.col(i);
    tridiagonal(i, i) = alpha;
    if (i > 0) {
      tridiagonal(i - 1, i) = beta;
      tridiagonal(i, i - 1) = beta;
    }
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(
        tridiagonal.block(0, 0, i + 1, i + 1));
    curvature = solver.eigenvalues()(0);
    estimate = basis.leftCols(i + 1) * solver.eigenvectors().col(0);
    eonc::safemath::safe_normalize_inplace(estimate);
    beta = direction.norm();
    const double change = eonc::safemath::safe_div(std::abs(curvature - previous),
                                                   std::abs(previous), 1.0);
    previous = curvature;
    if (beta <= 1e-10 * std::abs(alpha)) {
      break;
    }
    if (i >= 1 && change < params.lanczos_options().tolerance) {
      break;
    }
  }
  matter.setPositions(savedPos);
  matter.setCell(savedCell);
  return {estimate, curvature};
}

void appendTensor(JobResultEnvelope &env, std::string_view prefix,
                  std::string_view name, const Matrix3d &tensor) {
  const char *axis = "xyz";
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      env.extras.emplace_back(
          std::format("{}{}_{}{}", prefix, name, axis[row], axis[col]),
          tensor(row, col));
    }
  }
}

void writeTensor(std::ostream &out, std::string_view prefix,
                 std::string_view name, const Matrix3d &tensor) {
  const char *axis = "xyz";
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      out << std::format("{:.12e} {}{}_{}{}\n", tensor(row, col), prefix, name,
                         axis[row], axis[col]);
    }
  }
}

bool readStress(Matter &matter, Matrix3d &sigma, double &volume) {
  if (!matter.getPotential() || !matter.getPotential()->computesStress()) {
    return false;
  }
  try {
    sigma = matter.cauchyStress();
  } catch (const std::exception &) {
    return false;
  }
  if (!sigma.allFinite()) {
    return false;
  }
  volume = std::abs(matter.getCell().determinant());
  return volume > 0.0 && std::isfinite(volume);
}

} // namespace

Options optionsFrom(const Parameters &params) {
  const auto &stored = params.solid_state_options();
  Options options;
  options.enabled = stored.enabled;
  options.weight = stored.weight;
  options.external = Matrix3d::Identity() * stored.pressure;
  for (int i = 0; i < 9; ++i) {
    options.external(i / 3, i % 3) += stored.stress[static_cast<size_t>(i)];
  }
  options.external = symmetricPart(options.external);
  for (int i = 0; i < 6; ++i) {
    options.mask[static_cast<size_t>(i)] = stored.mask[static_cast<size_t>(i)] != 0;
  }
  return options;
}

bool movesCell(const Parameters &params) {
  const Options options = optionsFrom(params);
  if (!options.enabled) {
    return false;
  }
  if (!(options.weight > 0.0)) {
    throw std::invalid_argument("solid_state weight must be positive");
  }
  return activeCount(options) > 0;
}

Chart::Chart(Matter &matter, const Options &options, AtomMatrix *mode)
    : options_{options} {
  Matrix3d rotation = Matrix3d::Identity();
  Matrix3d cell = matter.getCell();
  AtomMatrix positions = matter.getPositions();
  if (!eonc::neb::orientCellLowerTriangular(cell, positions, &rotation)) {
    throw std::invalid_argument(
        "solid-state search could not put the cell in lower-triangular form");
  }
  matter.setCell(cell);
  if (positions.rows() > 0) {
    matter.setPositions(positions);
  }
  if (mode != nullptr && mode->rows() == positions.rows()) {
    *mode = (*mode) * rotation;
    for (long atom = 0; atom < matter.numberOfAtoms(); ++atom) {
      if (matter.getFixed(atom)) {
        mode->row(atom).setZero();
      }
    }
  }
  h0_ = cell;
  r0_ = matter.getPositions();
  referenceVolume_ = std::abs(h0_.determinant());
  if (!(referenceVolume_ > 0.0)) {
    throw std::invalid_argument("solid-state search needs a nonzero cell");
  }
  jacobian_ = eonc::neb::solidStateJacobian(
      referenceVolume_, static_cast<int>(matter.numberOfAtoms()),
      options_.weight);
  for (long atom = 0; atom < matter.numberOfAtoms(); ++atom) {
    if (!matter.getFixed(atom)) {
      free_.push_back(atom);
    }
  }
  for (int comp = 0; comp < 6; ++comp) {
    if (options_.mask[static_cast<size_t>(comp)]) {
      active_.push_back(comp);
    }
  }
}

Matrix3d Chart::deformation(const Matter &matter) const {
  return h0_.inverse() * matter.getCell();
}

void Chart::readStrain(const Matrix3d &deformation, Matrix3d &eps) const {
  eps = deformation - Matrix3d::Identity();
  zeroUpper(eps);
  for (int comp = 0; comp < 6; ++comp) {
    if (!options_.mask[static_cast<size_t>(comp)]) {
      eps(kStrainRows[comp], kStrainCols[comp]) = 0.0;
    }
  }
}

VectorXd Chart::coordinates(const Matter &matter) const {
  const Matrix3d deform = deformation(matter);
  Matrix3d eps;
  readStrain(deform, eps);
  const Matrix3d inverse = deform.inverse();
  VectorXd x(degreesOfFreedom());
  int slot = 0;
  for (long atom : free_) {
    const Eigen::RowVector3d reference =
        matter.getPositions().row(atom) * inverse;
    x.segment<3>(slot) = reference.transpose();
    slot += 3;
  }
  for (int comp : active_) {
    x[slot++] = jacobian_ * eps(kStrainRows[comp], kStrainCols[comp]);
  }
  return x;
}

void Chart::setCoordinates(Matter &matter, const VectorXd &x) const {
  if (x.size() != degreesOfFreedom()) {
    throw std::invalid_argument("solid-state coordinate size mismatch");
  }
  Matrix3d eps = Matrix3d::Zero();
  int slot = 3 * static_cast<int>(free_.size());
  for (int comp : active_) {
    eps(kStrainRows[comp], kStrainCols[comp]) = x[slot++] / jacobian_;
  }
  const Matrix3d deform = Matrix3d::Identity() + eps;
  Matrix3d cell = h0_ * deform;
  zeroUpper(cell);
  AtomMatrix positions = r0_;
  slot = 0;
  for (long atom = 0; atom < positions.rows(); ++atom) {
    Eigen::RowVector3d reference = r0_.row(atom);
    if (!matter.getFixed(atom)) {
      reference = x.segment<3>(slot).transpose();
      slot += 3;
    }
    positions.row(atom) = reference * deform;
  }
  matter.setCell(cell);
  matter.setPositions(positions);
}

VectorXd Chart::energyGradient(Matter &matter) const {
  const bool dirty = matter.needsForceUpdate();
  const AtomMatrix forces = matter.getForces();
  Matrix3d sigma;
  if (dirty) {
    if (!matter.getPotential() || !matter.getPotential()->computesStress()) {
      throw std::logic_error(
          "solid-state search requires a potential that reports the Cauchy stress");
    }
    sigma = matter.getPotential()->cauchyStress();
  } else {
    sigma = matter.cauchyStress();
  }
  if (!sigma.allFinite()) {
    throw std::runtime_error("potential returned a non-finite stress tensor");
  }
  const Matrix3d deform = deformation(matter);
  const Matrix3d deformT = deform.transpose();
  const double volume = std::abs(matter.getCell().determinant());
  if (!(volume > 0.0)) {
    throw std::runtime_error("solid-state search cell collapsed");
  }
  // dE/dε at fixed reference coordinates. The current right strain of
  // h = h0 (I+ε) is deform^{-1} dε, so the chain rule brings deform^{-T}.
  const Matrix3d dEnergy =
      volume * (deformT.inverse() * sigma) + referenceVolume_ * options_.external;

  VectorXd gradient(degreesOfFreedom());
  int slot = 0;
  for (long atom : free_) {
    const Eigen::RowVector3d row = -forces.row(atom) * deformT;
    gradient.segment<3>(slot) = row.transpose();
    slot += 3;
  }
  for (int comp : active_) {
    gradient[slot++] =
        dEnergy(kStrainRows[comp], kStrainCols[comp]) / jacobian_;
  }
  return gradient;
}

VectorXd Chart::difference(const VectorXd &a, const VectorXd &b) const {
  VectorXd delta = a - b;
  const Matrix3d inverse = h0_.inverse();
  for (int i = 0; i < static_cast<int>(free_.size()); ++i) {
    Eigen::RowVector3d frac = delta.segment<3>(3 * i).transpose() * inverse;
    for (int axis = 0; axis < 3; ++axis) {
      frac[axis] = wrapFrac(frac[axis]);
    }
    delta.segment<3>(3 * i) = (frac * h0_).transpose();
  }
  return delta;
}

AtomMatrix Chart::unpackAtoms(const VectorXd &joint) const {
  AtomMatrix mode = AtomMatrix::Zero(r0_.rows(), 3);
  int slot = 0;
  for (long atom : free_) {
    mode.row(atom) = joint.segment<3>(slot).transpose();
    slot += 3;
  }
  return mode;
}

Matrix3d Chart::unpackStrain(const VectorXd &joint) const {
  Matrix3d eps = Matrix3d::Zero();
  int slot = 3 * static_cast<int>(free_.size());
  for (int comp : active_) {
    eps(kStrainRows[comp], kStrainCols[comp]) = joint[slot++] / jacobian_;
  }
  return eps;
}

VectorXd Chart::packAtoms(const AtomMatrix &mode) const {
  VectorXd joint = VectorXd::Zero(degreesOfFreedom());
  int slot = 0;
  for (long atom : free_) {
    if (atom < mode.rows()) {
      joint.segment<3>(slot) = mode.row(atom).transpose();
    }
    slot += 3;
  }
  return joint;
}

std::shared_ptr<ObjectiveFunction>
makeVariableCellObjective(Matter &matter, const Parameters &params) {
  if (!matter.getPeriodic()) {
    throw std::invalid_argument(
        "variable-cell minimization requires periodic boundaries");
  }
  if (!matter.getPotential() || !matter.getPotential()->computesStress()) {
    throw std::logic_error(
        "variable-cell minimization requires a potential that reports the Cauchy stress");
  }
  const auto &opt = params.optimizer_options();
  if (opt.method == OptType::LBFGS && opt.lbfgs.precon != "none") {
    throw std::invalid_argument(
        "variable-cell minimization does not use an LBFGS pair preconditioner");
  }
  if (opt.lbfgs.project_rigid) {
    throw std::invalid_argument(
        "variable-cell minimization does not project rigid motion");
  }
  Chart chart(matter, optionsFrom(params));
  requireChart(chart);
  return std::make_shared<VariableCellObjective>(matter, std::move(chart),
                                                 params);
}

ModeEstimate lowestMode(Matter &matter, const Chart &chart, VectorXd direction,
                        const Parameters &params) {
  requireChart(chart);
  const std::string &method = params.saddle_search_options().minmode_method;
  if (method == LowestEigenmode::MINMODE_DIMER) {
    return dimerMode(matter, chart, std::move(direction), params);
  }
  if (method == LowestEigenmode::MINMODE_LANCZOS) {
    return lanczosMode(matter, chart, std::move(direction), params);
  }
  throw std::invalid_argument(
      "solid-state min-mode supports min_mode_method dimer and lanczos");
}

std::unique_ptr<SearchState>
runMinMode(Matter &matter, AtomMatrix mode, double reactantEnergy,
           const Parameters &params, long maxIterations) {
  if (params.saddle_search_options().prfo_refine.enabled) {
    throw std::invalid_argument(
        "solid-state min-mode does not run prfo_refine");
  }
  if (!matter.getPeriodic()) {
    throw std::invalid_argument(
        "solid-state dimer requires periodic boundaries");
  }
  if (!matter.getPotential() || !matter.getPotential()->computesStress()) {
    throw std::logic_error(
        "solid-state dimer requires a potential that reports the Cauchy stress");
  }
  const std::string &method = params.saddle_search_options().minmode_method;
  if (method != LowestEigenmode::MINMODE_DIMER &&
      method != LowestEigenmode::MINMODE_LANCZOS) {
    throw std::invalid_argument(
        "solid-state min-mode supports min_mode_method dimer and lanczos");
  }
  auto state = std::make_unique<SearchState>(
      Chart(matter, optionsFrom(params), &mode));
  requireChart(state->chart);
  VectorXd direction = state->chart.packAtoms(mode);
  if (direction.norm() < eonc::safemath::eps) {
    direction.setZero();
    direction[0] = 1.0;
  }
  eonc::safemath::safe_normalize_inplace(direction);
  auto objective =
      std::make_shared<ClimbObjective>(matter, state->chart, params);
  auto optim = eonc::helpers::create::mkOptim(
      objective, params.optimizer_options().method, params);
  const size_t callsBefore = matter.getPotential()->forceCallCounter.load();
  const double tol = params.saddle_search_options().converged_force > 0.0
                         ? params.saddle_search_options().converged_force
                         : params.optimizer_options().converged_force;
  const double maxEnergy = params.saddle_search_options().max_energy;
  state->status = MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS;
  long steps = 0;
  for (long iter = 0; iter < std::max(1L, maxIterations); ++iter) {
    const ModeEstimate estimate =
        lowestMode(matter, state->chart, direction, params);
    direction = estimate.mode;
    state->curvature = estimate.curvature;
    objective->setMode(direction, estimate.curvature);
    const double residual = state->chart.generalizedForce(matter).norm();
    if (residual < tol) {
      state->status = estimate.curvature < 0.0
                          ? MinModeSaddleSearch::STATUS_GOOD
                          : MinModeSaddleSearch::STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE;
      steps = iter;
      break;
    }
    if (iter + 1 >= maxIterations) {
      steps = iter + 1;
      break;
    }
    const VectorXd saved = state->chart.coordinates(matter);
    const int stepStatus = optim->step(params.optimizer_options().max_move);
    ++steps;
    if (stepStatus < 0) {
      state->status = MinModeSaddleSearch::STATUS_OPTIMIZER_ERROR;
      break;
    }
    const double energy = matter.getPotentialEnergy();
    if (!std::isfinite(energy)) {
      state->chart.setCoordinates(matter, saved);
      state->status = MinModeSaddleSearch::STATUS_POTENTIAL_FAILED;
      break;
    }
    if (energy - reactantEnergy > maxEnergy) {
      state->chart.setCoordinates(matter, saved);
      state->status = MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY;
      break;
    }
    if (objective->isConverged() && estimate.curvature < 0.0) {
      state->status = MinModeSaddleSearch::STATUS_GOOD;
      break;
    }
  }
  state->mode = direction;
  state->iteration = static_cast<int>(steps);
  state->forceCalls = static_cast<int>(
      matter.getPotential()->forceCallCounter.load() - callsBefore);
  return state;
}

void placeAlongMode(const SearchState &state, Matter &image, double scale) {
  VectorXd x = state.chart.coordinates(image);
  x += scale * state.mode;
  state.chart.setCoordinates(image, x);
}

void appendElasticDipole(JobResultEnvelope &env, Matter &matter,
                         std::string_view prefix) {
  Matrix3d sigma;
  double volume = 0.0;
  if (!readStress(matter, sigma, volume)) {
    return;
  }
  appendTensor(env, prefix, "stress", sigma);
  env.extras.emplace_back(std::format("{}cell_volume", prefix), volume);
  appendTensor(env, prefix, "elastic_dipole", volume * sigma);
}

void writeElasticDipole(std::ostream &out, Matter &matter,
                        std::string_view prefix) {
  Matrix3d sigma;
  double volume = 0.0;
  if (!readStress(matter, sigma, volume)) {
    return;
  }
  writeTensor(out, prefix, "stress", sigma);
  out << std::format("{:.12e} {}cell_volume\n", volume, prefix);
  writeTensor(out, prefix, "elastic_dipole", volume * sigma);
}

} // namespace eonc::solid

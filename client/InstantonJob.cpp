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
#include "eon/InstantonJob.h"
#include "eon/ConFileIO.h"
#include "eon/EonLogger.h"
#include "eon/Hessian.h"
#include "eon/JobResult.h"
#include "eon/Matter.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"
#include "eon/Tunneling.h"

#include <Eigen/QR>
#include <Eigen/SVD>

#include <array>
#include <cmath>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc {

namespace {

/// One unit of imaginary time, sqrt(amu Angstrom^2 / eV), in fs.
constexpr double kTimeUnitFs = 10.180505717871193;

/// ||H r|| / (||H||_F ||r||) at or below this is a rotational zero mode.
/// A real curvature sits near the scale of ||H||; a finite-difference null
/// vector does not.
constexpr double kRotationZero = 1e-2;

/// Mass-weighted coordinates over the free atoms, measured from a reference
/// structure under its minimum image.
class MassWeighted {
public:
  explicit MassWeighted(const Matter &reference)
      : ref_(reference) {
    for (long i = 0; i < reference.numberOfAtoms(); ++i) {
      if (reference.getFixed(i)) {
        continue;
      }
      const double m = reference.getMass(i);
      if (!(m > 0.0)) {
        throw std::invalid_argument("instanton: a free atom without a mass");
      }
      free_.push_back(i);
      sqrtMass_.push_back(std::sqrt(m));
    }
    if (free_.empty()) {
      throw std::invalid_argument("instanton: every atom is fixed");
    }
  }
  long dimension() const { return 3 * static_cast<long>(free_.size()); }
  VectorXi freeAtoms() const {
    VectorXi out(static_cast<long>(free_.size()));
    for (size_t k = 0; k < free_.size(); ++k) {
      out(static_cast<long>(k)) = static_cast<int>(free_[k]);
    }
    return out;
  }
  VectorXd toQ(const Matter &m) const {
    const AtomMatrix d = ref_.pbc(m.getPositions() - ref_.getPositions());
    VectorXd q(dimension());
    for (size_t k = 0; k < free_.size(); ++k) {
      for (int c = 0; c < 3; ++c) {
        q(static_cast<long>(3 * k) + c) = sqrtMass_[k] * d(free_[k], c);
      }
    }
    return q;
  }
  void place(const VectorXd &q, Matter &m) const {
    AtomMatrix r = ref_.getPositions();
    for (size_t k = 0; k < free_.size(); ++k) {
      for (int c = 0; c < 3; ++c) {
        r(free_[k], c) += q(static_cast<long>(3 * k) + c) / sqrtMass_[k];
      }
    }
    m.setPositions(r);
  }
  /// Three translations, then three rotations about the centre of mass, in
  /// these coordinates. A rotation column is zero for a linear molecule.
  MatrixXd rigidGenerators(const Matter &m) const {
    const long n = dimension();
    MatrixXd b = MatrixXd::Zero(n, 6);
    const AtomMatrix r = m.getPositions();
    Eigen::RowVector3d com = Eigen::RowVector3d::Zero();
    double total = 0.0;
    for (size_t k = 0; k < free_.size(); ++k) {
      const double w = sqrtMass_[k] * sqrtMass_[k];
      com += w * r.row(free_[k]);
      total += w;
    }
    com /= total;
    for (size_t k = 0; k < free_.size(); ++k) {
      const long i = static_cast<long>(3 * k);
      const Eigen::Vector3d x = (r.row(free_[k]) - com).transpose();
      for (int c = 0; c < 3; ++c) {
        b(i + c, c) = sqrtMass_[k];
        Eigen::Vector3d e = Eigen::Vector3d::Zero();
        e(c) = 1.0;
        b.block(i, 3 + c, 3, 1) = sqrtMass_[k] * e.cross(x);
      }
    }
    return b;
  }
  /// Which of the three rotation generators are zero modes of hess.
  /// ||H r|| small against ||H||, not whether the cell is periodic: a
  /// cluster in a box is periodic and still free to rotate.
  void markRotationZeroModes(const MatrixXd &hess, const MatrixXd &generators,
                             std::array<bool, 3> &keep,
                             std::array<double, 3> &residual) const {
    const MatrixXd h = 0.5 * (hess + hess.transpose());
    const double hn = h.norm();
    for (int c = 0; c < 3; ++c) {
      const VectorXd r = generators.col(3 + c);
      const double rn = r.norm();
      if (!(rn > 0.0)) {
        keep[static_cast<size_t>(c)] = false;
        residual[static_cast<size_t>(c)] = 0.0;
        continue;
      }
      const double rel = hn > 0.0 ? (h * r).norm() / (hn * rn) : 0.0;
      residual[static_cast<size_t>(c)] = rel;
      keep[static_cast<size_t>(c)] = rel <= kRotationZero;
    }
  }
  /// With no atom fixed, an orthonormal basis of the rigid motions of m:
  /// three translations, plus each rotation `rotations` marks as a zero
  /// mode. Empty when an atom is fixed. Rank-revealing, so a linear
  /// molecule keeps two rotations.
  MatrixXd rigidBasis(const Matter &m,
                      const std::array<bool, 3> &rotations) const {
    if (static_cast<long>(free_.size()) != m.numberOfAtoms()) {
      return {};
    }
    const MatrixXd g = rigidGenerators(m);
    const long n = g.rows();
    std::vector<int> cols{0, 1, 2};
    for (int c = 0; c < 3; ++c) {
      if (rotations[static_cast<size_t>(c)]) {
        cols.push_back(3 + c);
      }
    }
    MatrixXd b(n, static_cast<long>(cols.size()));
    for (size_t k = 0; k < cols.size(); ++k) {
      b.col(static_cast<long>(k)) = g.col(cols[k]);
    }
    const Eigen::ColPivHouseholderQR<MatrixXd> qr(b);
    const long rank = qr.rank();
    if (rank <= 0) {
      return {};
    }
    return qr.householderQ() * MatrixXd::Identity(n, rank);
  }
  VectorXd gradient(const AtomMatrix &forces) const {
    VectorXd g(dimension());
    for (size_t k = 0; k < free_.size(); ++k) {
      for (int c = 0; c < 3; ++c) {
        g(static_cast<long>(3 * k) + c) = -forces(free_[k], c) / sqrtMass_[k];
      }
    }
    return g;
  }

private:
  const Matter &ref_;
  std::vector<long> free_;
  std::vector<double> sqrtMass_;
};

/// With no atom fixed, removes the rigid part of m's displacement from
/// ref: the mass-weighted mean (translation), and for a cluster the
/// mass-weighted best rotation (Kabsch). A path that carried either would
/// pay kinetic action for motion that costs no energy.
void alignRigid(const Matter &ref, Matter &m) {
  const long n = ref.numberOfAtoms();
  for (long i = 0; i < n; ++i) {
    if (ref.getFixed(i)) {
      return;
    }
  }
  AtomMatrix d = ref.pbc(m.getPositions() - ref.getPositions());
  VectorXd w(n);
  for (long i = 0; i < n; ++i) {
    w(i) = ref.getMass(i);
  }
  const double total = w.sum();
  const Eigen::RowVector3d shift = (w.transpose() * d) / total;
  d.rowwise() -= shift;
  if (!ref.getPeriodic()) {
    const Eigen::RowVector3d com = (w.transpose() * ref.getPositions()) / total;
    AtomMatrix x = ref.getPositions();
    x.rowwise() -= com;
    const AtomMatrix y = x + d;
    const Eigen::Matrix3d h = y.transpose() * w.asDiagonal() * x;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(h, Eigen::ComputeFullU |
                                                 Eigen::ComputeFullV);
    Eigen::Matrix3d fix = Eigen::Matrix3d::Identity();
    fix(2, 2) = (svd.matrixV() * svd.matrixU().transpose()).determinant() < 0.0
                    ? -1.0
                    : 1.0;
    const Eigen::Matrix3d rot = svd.matrixV() * fix * svd.matrixU().transpose();
    d = (y * rot.transpose()) - x;
  }
  m.setPositions(ref.getPositions() + d);
}

} // namespace

namespace {

/// Mode rate: the ring-polymer instanton through the saddle out of the
/// reactant, and its thermal rate.
std::vector<std::string>
runRate(const Parameters &params, const std::shared_ptr<Potential> &pot,
        const Matter &reactant, const MassWeighted &mw,
        const tunneling::BatchPotential &evaluate,
        const std::function<MatrixXd(const VectorXd &)> &hessianAt,
        const std::array<bool, 3> &rotationZero,
        const std::array<double, 3> &rotationResidual) {
  const auto &o = params.instanton_options();
  const std::string resultsFile = "results.dat";
  const std::string pathFile = "instanton.con";
  std::vector<std::string> returnFiles{resultsFile};

  Matter saddle(pot, params);
  if (!io::io_ok(saddle.con2matter(o.saddle_filename))) {
    throw std::runtime_error("instanton: cannot read " + o.saddle_filename);
  }
  if (saddle.numberOfAtoms() != reactant.numberOfAtoms()) {
    throw std::runtime_error(
        "instanton: the saddle and the reactant differ in atom count");
  }
  if (!(o.temperature > 0.0)) {
    throw std::invalid_argument("instanton: mode rate needs [Instanton] "
                                "temperature in K");
  }
  alignRigid(reactant, saddle);
  const long n = mw.dimension();
  const VectorXd qSaddle = mw.toQ(saddle);
  const double vReactant = Matter(reactant).getPotentialEnergy();
  const double vSaddle = saddle.getPotentialEnergy();
  const MatrixXd hReactant = hessianAt(VectorXd::Zero(n));
  const MatrixXd hSaddle = hessianAt(qSaddle);
  const double tc = tunneling::crossoverTemperature(hSaddle);
  const long rigidModes = mw.rigidBasis(reactant, rotationZero).cols();
  const double beta = 1.0 / (tunneling::kBoltzmann * o.temperature);
  EONC_LOG_INFO("[Instanton] rate: {} beads at {:.4g} K, crossover {:.4g} K, "
                "barrier {:.6f} eV, {} degrees of freedom, {} rigid modes",
                o.beads, o.temperature, tc, vSaddle - vReactant, n, rigidModes);
  EONC_LOG_INFO("[Instanton] rotation residuals {:.3g}, {:.3g}, {:.3g}",
                rotationResidual[0], rotationResidual[1], rotationResidual[2]);

  std::vector<std::pair<std::string, double>> extras{
      {"instanton_temperature_K", o.temperature},
      {"instanton_crossover_K", tc},
      {"barrier_classical", vSaddle - vReactant}};
  auto write = [&](RunStatus status) {
    auto env = JobResultEnvelope::fromMinimization(
        status, params.potential_options().potential,
        PotRegistry::get().total_force_calls(), false, 0.0);
    env.job_type = "instanton";
    env.extras.emplace_back("force_calls",
                            static_cast<double>(env.force_calls));
    for (const auto &kv : extras) {
      env.extras.push_back(kv);
    }
    env.writeResultsDat(resultsFile);
  };

  if (!(o.temperature < tc)) {
    // The ring has collapsed onto the saddle. Above the crossover the
    // rate is the parabolic barrier correction, which this job does not
    // evaluate.
    EONC_LOG_ERROR("[Instanton] {:.4g} K is at or above the crossover "
                   "temperature {:.4g} K; the parabolic barrier correction "
                   "applies and this job does not evaluate it",
                   o.temperature, tc);
    write(RunStatus::FAIL_POTENTIAL_FAILED);
    return returnFiles;
  }

  tunneling::RateInstantonOptions ro;
  ro.beads = o.beads;
  ro.maxIterations = o.max_iterations;
  ro.forceTolerance = o.force_tolerance;
  tunneling::RateInstanton inst = tunneling::optimizeRateInstanton(
      qSaddle, hSaddle, beta, {}, evaluate, ro);
  EONC_LOG_INFO("[Instanton] ring U_N {:.6f} eV after {} iterations{}",
                inst.ringPotential, inst.iterations,
                inst.converged ? "" : " (not converged)");

  bool rateOk = false;
  if (inst.converged) {
    // Bead Hessians on every stride-th bead of the ring, linear in
    // between, wrapping from the last anchor back to bead 0.
    const long stride = std::max<long>(1, o.hessian_stride);
    const long N = o.beads;
    std::map<long, MatrixXd> anchors;
    auto anchor = [&](long j) -> const MatrixXd & {
      auto it = anchors.find(j);
      if (it == anchors.end()) {
        it = anchors.emplace(j, hessianAt(inst.beads[static_cast<size_t>(j)]))
                 .first;
      }
      return it->second;
    };
    auto beadHessian = [&](long j, const VectorXd &) -> MatrixXd {
      const long lo = (j / stride) * stride;
      if (j == lo) {
        return anchor(lo);
      }
      const long hi = lo + stride < N ? lo + stride : 0;
      const long span = (hi == 0 ? N : hi) - lo;
      const double t = static_cast<double>(j - lo) / static_cast<double>(span);
      return (1.0 - t) * anchor(lo) + t * anchor(hi);
    };
    try {
      tunneling::instantonRate(inst, beadHessian, hReactant, vReactant, hSaddle,
                               vSaddle, rigidModes);
      rateOk = std::isfinite(inst.logRate) && inst.negativeModes == 1;
      if (inst.negativeModes != 1) {
        EONC_LOG_ERROR("[Instanton] the ring Hessian has {} negative modes, "
                       "not one: the ring is not a first-order saddle of U_N",
                       inst.negativeModes);
      }
    } catch (const std::runtime_error &ex) {
      EONC_LOG_ERROR("[Instanton] {}", ex.what());
    }
  }

  // ln(k s): the rate itself underflows a double for deep tunnelling.
  const double logSecond = std::log(tunneling::kTimeUnitSeconds);
  Matter frame(reactant);
  for (size_t j = 0; j < inst.beads.size(); ++j) {
    mw.place(inst.beads[j], frame);
    io::ConFrameMetadata meta;
    meta.frame_index = static_cast<uint64_t>(j);
    meta.energy = inst.energies[j];
    meta.write_con_forces = false;
    meta.scalars = {{"imaginary_time_fs", static_cast<double>(j) * inst.betaN *
                                              tunneling::kHbar * kTimeUnitFs}};
    if (j == 0) {
      meta.scalars.push_back({"instanton_temperature_K", o.temperature});
      meta.scalars.push_back({"instanton_crossover_K", tc});
      meta.scalars.push_back(
          {"instanton_converged", inst.converged ? 1.0 : 0.0});
      if (rateOk) {
        meta.scalars.push_back(
            {"rate_instanton_log", inst.logRate - logSecond});
        meta.scalars.push_back(
            {"barrier_effective_instanton", inst.effectiveBarrier});
      }
    }
    if (!io::io_ok(frame.matter2con(pathFile, j > 0, &meta))) {
      throw std::runtime_error("instanton: cannot write " + pathFile);
    }
  }
  returnFiles.push_back(pathFile);

  extras.emplace_back("instanton_iterations",
                      static_cast<double>(inst.iterations));
  extras.emplace_back("instanton_ring_potential", inst.ringPotential);
  extras.emplace_back("instanton_bN", inst.bN);
  if (rateOk) {
    extras.emplace_back("rate_instanton", inst.rate);
    extras.emplace_back("rate_instanton_log", inst.logRate - logSecond);
    extras.emplace_back("rate_htst", inst.classicalRate);
    extras.emplace_back("rate_htst_log", inst.classicalLogRate - logSecond);
    extras.emplace_back("barrier_effective_instanton", inst.effectiveBarrier);
    extras.emplace_back("instanton_negative_modes",
                        static_cast<double>(inst.negativeModes));
    extras.emplace_back("instanton_zero_mode", inst.zeroEigenvalue);
  }
  write(rateOk ? RunStatus::GOOD
               : (inst.converged ? RunStatus::FAIL_POTENTIAL_FAILED
                                 : RunStatus::FAIL_MAX_ITERATIONS));
  return returnFiles;
}

} // namespace

std::vector<std::string> InstantonJob::run(void) {
  const auto &o = params.instanton_options();
  std::vector<std::string> returnFiles;
  const std::string resultsFile = "results.dat";
  const std::string pathFile = "instanton.con";
  returnFiles.push_back(resultsFile);

  auto reactant = std::make_unique<Matter>(pot, params);
  if (!io::io_ok(reactant->con2matter(o.reactant_filename))) {
    throw std::runtime_error("instanton: cannot read " + o.reactant_filename);
  }
  const MassWeighted mw(*reactant);
  const long n = mw.dimension();
  const VectorXd qStart = VectorXd::Zero(n);

  // Bead evaluations: one batch per call, through forceBatch when the
  // potential spreads a batch over calculators.
  std::vector<std::unique_ptr<Matter>> pool;
  auto evaluate = [&](const std::vector<VectorXd> &q, std::vector<double> &v,
                      std::vector<VectorXd> &grad) {
    while (pool.size() < q.size()) {
      pool.push_back(std::make_unique<Matter>(*reactant));
    }
    for (size_t j = 0; j < q.size(); ++j) {
      mw.place(q[j], *pool[j]);
    }
    v.resize(q.size());
    grad.resize(q.size());
    if (pot->supportsBatchEvaluation() && q.size() > 1) {
      const long atoms = reactant->numberOfAtoms();
      std::vector<VectorXi> nrs(q.size());
      std::vector<Matrix3d> boxes(q.size());
      std::vector<const double *> posPtr, boxPtr;
      std::vector<const int *> nrsPtr;
      std::vector<double *> frcPtr;
      for (size_t j = 0; j < q.size(); ++j) {
        nrs[j] = pool[j]->getAtomicNrs();
        boxes[j] = pool[j]->getPeriodic() ? pool[j]->getCell()
                                          : Matrix3d::Zero().eval();
      }
      for (size_t j = 0; j < q.size(); ++j) {
        posPtr.push_back(pool[j]->getPositions().data());
        nrsPtr.push_back(nrs[j].data());
        frcPtr.push_back(pool[j]->forcesData());
        boxPtr.push_back(boxes[j].data());
      }
      std::vector<double> energies(q.size()), variances(q.size());
      pot->forceBatch(static_cast<long>(q.size()), atoms, posPtr.data(),
                      nrsPtr.data(), frcPtr.data(), energies.data(),
                      variances.data(), boxPtr.data());
      for (size_t j = 0; j < q.size(); ++j) {
        pool[j]->setComputedPotential(energies[j], variances[j]);
      }
    }
    for (size_t j = 0; j < q.size(); ++j) {
      v[j] = pool[j]->getPotentialEnergy();
      grad[j] = mw.gradient(pool[j]->getForces());
    }
  };

  // The first call is the reactant. Its Hessian decides which rotations
  // are zero modes; later beads reuse that decision.
  std::array<bool, 3> rotationZero{{false, false, false}};
  std::array<double, 3> rotationResidual{{0.0, 0.0, 0.0}};
  bool rotationsKnown = false;
  auto hessianAt = [&](const VectorXd &q) {
    Matter m(*reactant);
    mw.place(q, m);
    Hessian h(params, &m);
    h.writeHessianFile(false);
    MatrixXd out = h.getHessian(&m, mw.freeAtoms());
    if (out.rows() != n) {
      throw std::runtime_error("instanton: a bead Hessian failed");
    }
    if (!rotationsKnown) {
      mw.markRotationZeroModes(out, mw.rigidGenerators(m), rotationZero,
                               rotationResidual);
      rotationsKnown = true;
    }
    // A finite-difference Hessian of a free structure has small nonzero
    // rigid eigenvalues of either sign; project them to zero.
    const MatrixXd rigid = mw.rigidBasis(m, rotationZero);
    if (rigid.cols() > 0) {
      const MatrixXd p = MatrixXd::Identity(n, n) - rigid * rigid.transpose();
      out = p * out * p;
    }
    return out;
  };

  if (o.mode == "rate") {
    return runRate(params, pot, *reactant, mw, evaluate, hessianAt,
                   rotationZero, rotationResidual);
  }

  auto product = std::make_unique<Matter>(pot, params);
  if (!io::io_ok(product->con2matter(o.product_filename))) {
    throw std::runtime_error("instanton: cannot read " + o.product_filename);
  }
  if (reactant->numberOfAtoms() != product->numberOfAtoms()) {
    throw std::runtime_error("instanton: the minima differ in atom count");
  }
  alignRigid(*reactant, *product);
  const VectorXd qEnd = mw.toQ(*product);

  // Starting path: a band from file, else the straight line.
  std::vector<VectorXd> guess;
  if (!o.initial_path.empty()) {
    const auto frames = readcon::read_all_frames(o.initial_path);
    for (const auto &frame : frames) {
      Matter m(*reactant);
      if (!io::io_ok(io::con2matter(m, frame))) {
        throw std::runtime_error("instanton: cannot read " + o.initial_path);
      }
      alignRigid(*reactant, m);
      guess.push_back(mw.toQ(m));
    }
    if (guess.size() >= 2) {
      guess.front() = qStart;
      guess.back() = qEnd;
    }
  }

  const MatrixXd hStart = hessianAt(qStart);
  const MatrixXd hEnd = hessianAt(qEnd);
  const double omega = tunneling::pathOmega(hStart, hEnd, qStart, qEnd);
  const double betaHbar = o.beta_hbar_omega / omega;
  tunneling::InstantonOptions opt;
  opt.beads = o.beads;
  opt.betaHbarOmega = o.beta_hbar_omega;
  opt.maxIterations = o.max_iterations;
  opt.forceTolerance = o.force_tolerance;
  EONC_LOG_INFO("[Instanton] {} beads over beta hbar = {:.4f} fs, {} degrees "
                "of freedom",
                o.beads, betaHbar * kTimeUnitFs, n);

  tunneling::Instanton inst = tunneling::optimizeInstanton(
      qStart, qEnd, betaHbar, guess, evaluate, opt);
  EONC_LOG_INFO("[Instanton] action {:.6f} after {} iterations{}", inst.action,
                inst.iterations, inst.converged ? "" : " (not converged)");

  bool splitOk = false;
  std::string failure;
  // beta |delta|: the propagator ratio reads delta0 only when the wells
  // lie within a small fraction of kB T of each other.
  const double betaAsymmetry =
      std::abs(inst.asymmetry) * betaHbar / tunneling::kHbar;
  if (inst.converged && !inst.symmetricEnough) {
    EONC_LOG_WARNING("[Instanton] beta |delta| = {:.3g}: the wells differ by "
                     "{:.4g} eV, too far for the splitting; the path and "
                     "action are written, the splitting is not",
                     betaAsymmetry, inst.asymmetry);
  }
  if (inst.converged && inst.symmetricEnough) {
    const long stride = std::max<long>(1, o.hessian_stride);
    const long P = o.beads;
    std::map<long, MatrixXd> anchors;
    auto anchor = [&](long j) -> const MatrixXd & {
      auto it = anchors.find(j);
      if (it == anchors.end()) {
        it = anchors.emplace(j, hessianAt(inst.path[static_cast<size_t>(j)]))
                 .first;
      }
      return it->second;
    };
    auto beadHessian = [&](long j, const VectorXd &) -> MatrixXd {
      if (stride == 1) {
        return anchor(j);
      }
      const long lo = 1 + ((j - 1) / stride) * stride;
      const long hi = std::min(lo + stride, P - 1);
      if (j == lo || hi == lo) {
        return anchor(lo);
      }
      const double t =
          static_cast<double>(j - lo) / static_cast<double>(hi - lo);
      return (1.0 - t) * anchor(lo) + t * anchor(hi);
    };
    try {
      tunneling::instantonSplitting(inst, beadHessian, hStart, hEnd);
      splitOk = std::isfinite(inst.delta0);
    } catch (const std::runtime_error &ex) {
      failure = ex.what();
      EONC_LOG_ERROR("[Instanton] {}", failure);
    }
  }

  // The path, one frame per bead; the splitting on the first frame.
  const double kelvin = tunneling::kHbar / (tunneling::kBoltzmann * betaHbar);
  Matter frame(*reactant);
  for (size_t j = 0; j < inst.path.size(); ++j) {
    mw.place(inst.path[j], frame);
    io::ConFrameMetadata meta;
    meta.frame_index = static_cast<uint64_t>(j);
    meta.energy = inst.energies[j];
    meta.write_con_forces = false;
    meta.scalars = {{"imaginary_time_fs",
                     static_cast<double>(j) * inst.dtau * kTimeUnitFs}};
    if (j == 0) {
      meta.scalars.push_back({"instanton_action", inst.action});
      meta.scalars.push_back(
          {"instanton_beta_hbar_fs", betaHbar * kTimeUnitFs});
      meta.scalars.push_back({"instanton_temperature_K", kelvin});
      meta.scalars.push_back(
          {"instanton_converged", inst.converged ? 1.0 : 0.0});
      meta.scalars.push_back({"tunnel_asymmetry", inst.asymmetry});
      if (splitOk) {
        meta.scalars.push_back({"tunnel_splitting_instanton", inst.delta0});
        meta.scalars.push_back(
            {"tls_energy_instanton", std::hypot(inst.asymmetry, inst.delta0)});
        meta.scalars.push_back({"instanton_zero_mode", inst.zeroMode});
        meta.scalars.push_back(
            {"instanton_mode_separation", inst.modeSeparation});
        meta.scalars.push_back(
            {"instanton_symmetric", inst.symmetricEnough ? 1.0 : 0.0});
      }
    }
    if (!io::io_ok(frame.matter2con(pathFile, j > 0, &meta))) {
      throw std::runtime_error("instanton: cannot write " + pathFile);
    }
  }
  returnFiles.push_back(pathFile);

  // A converged path between wells too far apart is a result, not a
  // failure: the flags say why no splitting was written.
  const bool good = splitOk || (inst.converged && !inst.symmetricEnough);
  const auto status = good ? RunStatus::GOOD
                           : (inst.converged ? RunStatus::FAIL_POTENTIAL_FAILED
                                             : RunStatus::FAIL_MAX_ITERATIONS);
  auto env = JobResultEnvelope::fromMinimization(
      status, params.potential_options().potential,
      PotRegistry::get().total_force_calls(), false, 0.0);
  env.job_type = "instanton";
  env.extras.emplace_back("force_calls", static_cast<double>(env.force_calls));
  env.extras.emplace_back("instanton_iterations",
                          static_cast<double>(inst.iterations));
  env.extras.emplace_back("instanton_action", inst.action);
  env.extras.emplace_back("instanton_temperature_K", kelvin);
  env.extras.emplace_back("tunnel_asymmetry", inst.asymmetry);
  env.extras.emplace_back("instanton_beta_asymmetry", betaAsymmetry);
  env.extras.emplace_back("instanton_symmetric",
                          inst.symmetricEnough ? 1.0 : 0.0);
  if (splitOk) {
    env.extras.emplace_back("tunnel_splitting_instanton", inst.delta0);
    env.extras.emplace_back("instanton_mode_separation", inst.modeSeparation);
  }
  env.writeResultsDat(resultsFile);
  return returnFiles;
}

} // namespace eonc

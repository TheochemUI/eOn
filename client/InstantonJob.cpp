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
#include "eon/PathIntegral.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"
#include "eon/Tunneling.h"

#include <Eigen/QR>
#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
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

/// One bead between every neighbour. A ring of N becomes a ring of 2N.
std::vector<VectorXd> doubledRing(const std::vector<VectorXd> &ring) {
  const size_t n = ring.size();
  std::vector<VectorXd> fine(2 * n);
  for (size_t j = 0; j < n; ++j) {
    fine[2 * j] = ring[j];
    fine[2 * j + 1] = 0.5 * (ring[j] + ring[(j + 1) % n]);
  }
  return fine;
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
  if (o.hessian_final != "recomputed") {
    throw std::invalid_argument("instanton: hessian_final must be recomputed");
  }
  std::vector<double> temperatures = o.temperatures;
  if (temperatures.empty()) {
    temperatures.push_back(o.temperature);
  }
  for (const double t : temperatures) {
    if (!(t > 0.0)) {
      throw std::invalid_argument(
          "instanton: mode rate needs [Instanton] temperature or "
          "temperatures in K");
    }
  }
  // Highest first, so each ring can start the colder one.
  std::sort(temperatures.begin(), temperatures.end(), std::greater<>());
  alignRigid(reactant, saddle);
  const long n = mw.dimension();
  const VectorXd qSaddle = mw.toQ(saddle);
  const double vReactant = Matter(reactant).getPotentialEnergy();
  const double vSaddle = saddle.getPotentialEnergy();
  const MatrixXd hReactant = hessianAt(VectorXd::Zero(n));
  const MatrixXd hSaddle = hessianAt(qSaddle);
  const double tc = tunneling::crossoverTemperature(hSaddle);
  const long rigidModes = mw.rigidBasis(reactant, rotationZero).cols();
  if (temperatures.size() == 1) {
    EONC_LOG_INFO("[Instanton] rate: {} beads at {:.4g} K, crossover {:.4g} K, "
                  "barrier {:.6f} eV, {} degrees of freedom, {} rigid modes",
                  o.beads, temperatures.front(), tc, vSaddle - vReactant, n,
                  rigidModes);
  } else {
    EONC_LOG_INFO("[Instanton] rate: {} beads at {} temperatures from {:.4g} K "
                  "down to {:.4g} K, crossover {:.4g} K, barrier {:.6f} eV, "
                  "{} degrees of freedom, {} rigid modes",
                  o.beads, temperatures.size(), temperatures.front(),
                  temperatures.back(), tc, vSaddle - vReactant, n, rigidModes);
  }
  EONC_LOG_INFO("[Instanton] rotation residuals {:.3g}, {:.3g}, {:.3g}",
                rotationResidual[0], rotationResidual[1], rotationResidual[2]);

  std::vector<std::pair<std::string, double>> extras{
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

  // A band over the barrier seeds the ring and carries the
  // one-dimensional WKB rate. An empty path leaves the saddle-mode seed.
  std::vector<VectorXd> pathQ;
  std::vector<double> pathV;
  std::unique_ptr<tunneling::Profile> profile;
  double hwPath = 0.0;
  if (!o.initial_path.empty()) {
    const auto frames = readcon::read_all_frames(o.initial_path);
    bool haveEnergies = true;
    for (const auto &frame : frames) {
      Matter image(reactant);
      if (!io::io_ok(io::con2matter(image, frame))) {
        throw std::runtime_error("instanton: cannot read " + o.initial_path);
      }
      if (image.numberOfAtoms() != reactant.numberOfAtoms()) {
        throw std::runtime_error("instanton: " + o.initial_path +
                                 " differs in atom count");
      }
      alignRigid(reactant, image);
      pathQ.push_back(mw.toQ(image));
      const auto energy = frame.energy_opt();
      haveEnergies = haveEnergies && energy.has_value();
      pathV.push_back(energy.value_or(0.0));
    }
    if (pathQ.size() < 3) {
      throw std::runtime_error("instanton: " + o.initial_path +
                               " holds fewer than three frames");
    }
    if (!haveEnergies) {
      std::vector<VectorXd> grads;
      evaluate(pathQ, pathV, grads);
    }
    std::vector<double> arc(pathQ.size(), 0.0);
    for (size_t k = 1; k < pathQ.size(); ++k) {
      arc[k] = arc[k - 1] + (pathQ[k] - pathQ[k - 1]).norm();
    }
    try {
      profile = std::make_unique<tunneling::Profile>(std::move(arc), pathV);
    } catch (const std::invalid_argument &ex) {
      EONC_LOG_WARNING("[Instanton] {}; seeding from the saddle mode",
                       ex.what());
    }
    if (profile) {
      try {
        hwPath = tunneling::hbarOmega(tunneling::wellCurvature(*profile, true));
      } catch (const std::exception &ex) {
        EONC_LOG_WARNING("[Instanton] no reactant frequency along the path: {}",
                         ex.what());
      }
    }
  }

  const std::string tableFile = "rate_instanton.dat";
  std::ofstream table(tableFile);
  if (!table) {
    throw std::runtime_error("instanton: cannot write " + tableFile);
  }
  table << "# T_K T_c_K beads converged iterations U_N_eV negative_modes "
           "ln_k_per_s k_per_s ln_k_htst_per_s barrier_effective_eV "
           "ln_k_wkb_path_per_s ln_k_parabolic_per_s parabolic_factor\n";
  table << std::setprecision(10);
  returnFiles.push_back(tableFile);
  const double logSecond = std::log(tunneling::kTimeUnitSeconds);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  auto wkbAt = [&](double beta) {
    if (!(profile && hwPath > 0.0)) {
      return nan;
    }
    try {
      return tunneling::wkbLogRateAlongPath(*profile, beta, hwPath) - logSecond;
    } catch (const std::exception &ex) {
      EONC_LOG_WARNING("[Instanton] no WKB rate along the path: {}", ex.what());
      return nan;
    }
  };

  std::vector<VectorXd> ring;
  RunStatus status = RunStatus::GOOD;
  bool rateFailed = false;
  for (size_t ti = 0; ti < temperatures.size(); ++ti) {
    const double temperature = temperatures[ti];
    const double beta = 1.0 / (tunneling::kBoltzmann * temperature);
    const bool last = ti + 1 == temperatures.size();
    const double wkbLog = wkbAt(beta);
    if (!(temperature < tc)) {
      // The ring collapses onto the saddle. Above T_c the rate is the
      // parabolic factor times harmonic TST. At T_c the factor diverges.
      bool wrote = false;
      if (temperature > tc) {
        try {
          const double factor = tunneling::parabolicFactor(temperature, tc);
          const double logHtst = tunneling::harmonicTstLogRate(
              hReactant, hSaddle, beta, vSaddle - vReactant, rigidModes);
          const double logPar = logHtst + std::log(factor);
          const double kPar = std::exp(logPar) / tunneling::kTimeUnitSeconds;
          const double kHtst = std::exp(logHtst) / tunneling::kTimeUnitSeconds;
          EONC_LOG_INFO("[Instanton] {:.4g} K is above the crossover {:.4g} "
                        "K; parabolic factor {:.6g}, ln(k s) = {:.4f}",
                        temperature, tc, factor, logPar - logSecond);
          if (factor > 10.0) {
            EONC_LOG_WARNING(
                "[Instanton] parabolic factor {:.6g} is large; this close "
                "to the crossover a uniform theory is the finite rate",
                factor);
          }
          table << temperature << ' ' << tc << ' ' << o.beads
                << " 0 0 nan 0 nan nan " << (logHtst - logSecond) << " nan "
                << wkbLog << ' ' << (logPar - logSecond) << ' ' << factor
                << '\n';
          if (last) {
            extras.emplace_back("instanton_temperature_K", temperature);
            extras.emplace_back("parabolic_factor", factor);
            extras.emplace_back("rate_parabolic", kPar);
            extras.emplace_back("rate_parabolic_log", logPar - logSecond);
            extras.emplace_back("rate_htst", kHtst);
            extras.emplace_back("rate_htst_log", logHtst - logSecond);
            if (std::isfinite(wkbLog)) {
              extras.emplace_back("rate_wkb_path_log", wkbLog);
            }
          }
          wrote = true;
        } catch (const std::exception &ex) {
          EONC_LOG_ERROR("[Instanton] {}", ex.what());
        }
      }
      if (!wrote) {
        if (!(temperature > tc)) {
          EONC_LOG_ERROR("[Instanton] {:.4g} K is at the crossover temperature "
                         "{:.4g} K; the parabolic factor diverges there",
                         temperature, tc);
        }
        table << temperature << ' ' << tc << ' ' << o.beads
              << " 0 0 nan 0 nan nan nan nan " << wkbLog << " nan nan\n";
        rateFailed = true;
        status = RunStatus::FAIL_POTENTIAL_FAILED;
        if (last) {
          extras.emplace_back("instanton_temperature_K", temperature);
          if (std::isfinite(wkbLog)) {
            extras.emplace_back("rate_wkb_path_log", wkbLog);
          }
        }
      }
      continue;
    }

    tunneling::RateInstantonOptions ro;
    ro.beads = o.beads;
    ro.maxIterations = o.max_iterations;
    ro.forceTolerance = o.force_tolerance;
    ro.halfRing = o.half_ring;
    ro.energyShift = o.energy_shift;
    std::vector<VectorXd> guess = ring;
    if (guess.empty() && profile) {
      try {
        guess = tunneling::ringFromPath(pathQ, pathV, beta * tunneling::kHbar,
                                        o.beads);
        EONC_LOG_INFO("[Instanton] ring seeded from {} by the period condition",
                      o.initial_path);
      } catch (const std::invalid_argument &ex) {
        EONC_LOG_WARNING("[Instanton] {}; seeding from the saddle mode instead",
                         ex.what());
        guess.clear();
      }
    }
    long ladderIterations = 0;
    if (guess.empty() && o.bead_ladder && o.beads >= 16) {
      long coarse = o.beads / 4;
      if (coarse % 2 != 0) {
        ++coarse;
      }
      if (coarse < 4) {
        coarse = 4;
      }
      std::vector<VectorXd> rung;
      for (long nb = coarse; nb < o.beads; nb *= 2) {
        tunneling::RateInstantonOptions step = ro;
        step.beads = nb;
        const tunneling::RateInstanton rungInst =
            tunneling::optimizeRateInstanton(qSaddle, hSaddle, beta, rung,
                                             evaluate, step);
        ladderIterations += rungInst.iterations;
        EONC_LOG_INFO("[Instanton] ladder rung {} beads: U_N {:.6f} eV after "
                      "{} iterations{}",
                      nb, rungInst.ringPotential, rungInst.iterations,
                      rungInst.converged ? "" : " (not converged)");
        rung = rungInst.beads;
        if (static_cast<long>(rung.size()) != nb) {
          rung.clear();
          break;
        }
        rung = doubledRing(rung);
        if (static_cast<long>(rung.size()) > o.beads) {
          rung.clear();
          break;
        }
      }
      if (static_cast<long>(rung.size()) == o.beads) {
        guess = std::move(rung);
      }
    }
    tunneling::RateInstanton inst = tunneling::optimizeRateInstanton(
        qSaddle, hSaddle, beta, guess, evaluate, ro);
    inst.iterations += ladderIterations;
    EONC_LOG_INFO("[Instanton] {:.4g} K: ring U_N {:.6f} eV after {} "
                  "iterations{}",
                  temperature, inst.ringPotential, inst.iterations,
                  inst.converged ? "" : " (not converged)");

    bool rateOk = false;
    if (inst.converged) {
      // Bead Hessians on every stride-th bead of the ring, linear in
      // between, wrapping from the last anchor back to bead 0.
      const long stride = std::max<long>(1, o.hessian_stride);
      const long nBeads = o.beads;
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
        const long hi = lo + stride < nBeads ? lo + stride : 0;
        const long span = (hi == 0 ? nBeads : hi) - lo;
        const double t =
            static_cast<double>(j - lo) / static_cast<double>(span);
        return (1.0 - t) * anchor(lo) + t * anchor(hi);
      };
      try {
        tunneling::instantonRate(inst, beadHessian, hReactant,
                                 vReactant - o.energy_shift, hSaddle,
                                 vSaddle - o.energy_shift, rigidModes);
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
    if (rateOk) {
      EONC_LOG_INFO("[Instanton] {:.4g} K: ln(k s) = {:.4f}, harmonic TST "
                    "{:.4f}, effective barrier {:.4f} eV",
                    temperature, inst.logRate - logSecond,
                    inst.classicalLogRate - logSecond, inst.effectiveBarrier);
    }
    table << temperature << ' ' << tc << ' ' << o.beads << ' '
          << (inst.converged ? 1 : 0) << ' ' << inst.iterations << ' '
          << inst.ringPotential << ' ' << inst.negativeModes << ' ';
    if (rateOk) {
      table << inst.logRate - logSecond << ' ' << inst.rate << ' '
            << inst.classicalLogRate - logSecond << ' '
            << inst.effectiveBarrier;
    } else {
      table << "nan nan nan nan";
    }
    table << ' ' << wkbLog << " nan nan\n";

    std::vector<std::string> files;
    if (last) {
      files.push_back(pathFile);
    }
    if (temperatures.size() > 1) {
      std::ostringstream name;
      name << std::defaultfloat << std::setprecision(6) << "instanton_"
           << temperature << "K.con";
      files.push_back(name.str());
    }
    if (!inst.beads.empty()) {
      Matter frame(reactant);
      for (const auto &file : files) {
        for (size_t j = 0; j < inst.beads.size(); ++j) {
          mw.place(inst.beads[j], frame);
          io::ConFrameMetadata meta;
          meta.frame_index = static_cast<uint64_t>(j);
          meta.energy = inst.energies[j];
          meta.write_con_forces = false;
          meta.scalars = {
              {"imaginary_time_fs", static_cast<double>(j) * inst.betaN *
                                        tunneling::kHbar * kTimeUnitFs}};
          if (j == 0) {
            meta.scalars.push_back({"instanton_temperature_K", temperature});
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
          if (!io::io_ok(frame.matter2con(file, j > 0, &meta))) {
            throw std::runtime_error("instanton: cannot write " + file);
          }
        }
        returnFiles.push_back(file);
      }
    }
    ring = inst.beads;
    if (!last) {
      continue;
    }
    extras.emplace_back("instanton_temperature_K", temperature);
    extras.emplace_back("instanton_iterations",
                        static_cast<double>(inst.iterations));
    extras.emplace_back("instanton_ring_potential", inst.ringPotential);
    extras.emplace_back("instanton_bN", inst.bN);
    if (std::isfinite(wkbLog)) {
      extras.emplace_back("rate_wkb_path_log", wkbLog);
    }
    if (rateOk) {
      extras.emplace_back("rate_instanton", inst.rate);
      extras.emplace_back("rate_instanton_log", inst.logRate - logSecond);
      extras.emplace_back("rate_htst", inst.classicalRate);
      extras.emplace_back("rate_htst_log", inst.classicalLogRate - logSecond);
      extras.emplace_back("barrier_effective_instanton", inst.effectiveBarrier);
      extras.emplace_back("instanton_negative_modes",
                          static_cast<double>(inst.negativeModes));
      extras.emplace_back("instanton_zero_mode", inst.zeroEigenvalue);
    } else {
      rateFailed = true;
      status = inst.converged ? RunStatus::FAIL_POTENTIAL_FAILED
                              : RunStatus::FAIL_MAX_ITERATIONS;
    }
  }
  if (!rateFailed) {
    status = RunStatus::GOOD;
  }
  write(status);
  return returnFiles;
}

} // namespace

std::vector<std::string> InstantonJob::run(void) {
  const auto &o = params.instanton_options();
  pathintegral::requireTrotterSprings(o.springs, "instanton");
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

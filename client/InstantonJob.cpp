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
#include "eon/PIQTST.h"
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
  /// The mean of a set of beads, and the root-mean-square spread of each
  /// free atom's position about it along x, y and z in Angstrom (row-major
  /// over every atom of m, zero for a fixed atom).
  std::vector<double> spreadAbout(const std::vector<VectorXd> &beads,
                                  VectorXd &centroid, long atoms) const {
    centroid = VectorXd::Zero(dimension());
    for (const auto &q : beads) {
      centroid += q;
    }
    centroid /= static_cast<double>(std::max<size_t>(1, beads.size()));
    std::vector<double> out(static_cast<size_t>(3 * atoms), 0.0);
    for (size_t k = 0; k < free_.size(); ++k) {
      for (int c = 0; c < 3; ++c) {
        const long i = static_cast<long>(3 * k) + c;
        double ss = 0.0;
        for (const auto &q : beads) {
          const double d = q(i) - centroid(i);
          ss += d * d;
        }
        out[static_cast<size_t>(3 * free_[k] + c)] =
            std::sqrt(ss /
                      static_cast<double>(std::max<size_t>(1, beads.size()))) /
            sqrtMass_[k];
      }
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
  /// Which of the three rotation generators are zero modes of hess, by
  /// tunneling::rotationZeroModes: the rotation's curvature against the
  /// softest vibration, not whether the cell is periodic (a cluster in a box
  /// is periodic and still free to rotate).
  void markRotationZeroModes(const MatrixXd &hess, const MatrixXd &generators,
                             std::array<bool, 3> &keep,
                             std::array<double, 3> &residual) const {
    const tunneling::RotationZeroModes z =
        tunneling::rotationZeroModes(hess, generators);
    keep = z.zero;
    residual = z.residual;
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
  const std::vector<double> &sqrtMasses() const { return sqrtMass_; }
  /// Cartesian positions of the free atoms of the reference, 3 per atom.
  VectorXd referenceFree() const {
    const AtomMatrix r = ref_.getPositions();
    VectorXd out(dimension());
    for (size_t k = 0; k < free_.size(); ++k) {
      for (int c = 0; c < 3; ++c) {
        out(static_cast<long>(3 * k) + c) = r(free_[k], c);
      }
    }
    return out;
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

/// One frame at the beads' centroid with the readcon spreads section: the
/// delocalised configuration as centroid plus per-atom root-mean-square
/// spread, which a path-integral trajectory writes the same way.
void writeCentroid(const std::string &file, const std::vector<VectorXd> &beads,
                   const MassWeighted &mw, const Matter &reactant,
                   std::vector<io::ConMetadataValue> scalars) {
  VectorXd centroid;
  Matter frame(reactant);
  io::ConFrameMetadata meta;
  meta.spreads = mw.spreadAbout(beads, centroid, reactant.numberOfAtoms());
  mw.place(centroid, frame);
  meta.frame_index = 0;
  meta.write_con_forces = false;
  double largest = 0.0;
  for (const double s : meta.spreads) {
    largest = std::max(largest, s);
  }
  scalars.push_back({"beads", static_cast<double>(beads.size())});
  scalars.push_back({"spread_max", largest});
  meta.scalars = std::move(scalars);
  if (!io::io_ok(frame.matter2con(file, false, &meta))) {
    throw std::runtime_error("instanton: cannot write " + file);
  }
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

/// Steepest-descent path in mass-weighted coordinates out of the saddle
/// along both signs of its unstable mode, ordered from one end through the
/// saddle to the other, reactant end first. Each step is q - alpha g with alpha
/// from backtracking, starting at the inverse of the stiffest saddle curvature,
/// and no atom moves more than `cartStep` Angstrom. A side ends where the
/// gradient falls below `gradTol` or no step lowers the energy. It stands in
/// for a band when none is given.
void steepestDescentPath(const VectorXd &qSaddle, double vSaddle,
                         const MatrixXd &hSaddle,
                         const std::vector<double> &sqrtMass,
                         const tunneling::BatchPotential &evaluate,
                         std::vector<VectorXd> &pathQ,
                         std::vector<double> &pathV) {
  constexpr double cartStep = 0.01;
  constexpr double gradTol = 1e-3;
  constexpr long maxSteps = 4000;
  const long n = qSaddle.size();
  const Eigen::SelfAdjointEigenSolver<MatrixXd> es(
      0.5 * (hSaddle + hSaddle.transpose()));
  const VectorXd mode = es.eigenvectors().col(0);
  const double stiff = es.eigenvalues().cwiseAbs().maxCoeff();
  // The largest Cartesian move of a mass-weighted displacement.
  auto cartesian = [&](const VectorXd &dq) {
    double big = 0.0;
    for (long i = 0; i < n; ++i) {
      const size_t a = static_cast<size_t>(i / 3);
      const double m = a < sqrtMass.size() ? sqrtMass[a] : 1.0;
      big = std::max(big, std::abs(dq(i)) / m);
    }
    return big;
  };
  auto capped = [&](VectorXd dq) {
    const double big = cartesian(dq);
    if (big > cartStep) {
      dq *= cartStep / big;
    }
    return dq;
  };
  auto at = [&](const VectorXd &q, double &v, VectorXd &g) {
    std::vector<VectorXd> one{q};
    std::vector<double> vs;
    std::vector<VectorXd> gs;
    evaluate(one, vs, gs);
    if (vs.empty() || gs.empty() || !std::isfinite(vs[0]) ||
        !gs[0].array().isFinite().all()) {
      return false;
    }
    v = vs[0];
    g = gs[0];
    return true;
  };
  std::vector<std::vector<VectorXd>> sideQ(2);
  std::vector<std::vector<double>> sideV(2);
  for (int side = 0; side < 2; ++side) {
    const double big0 = cartesian(mode);
    VectorXd q = qSaddle + (side == 0 ? -1.0 : 1.0) * (cartStep / big0) * mode;
    double v = 0.0;
    VectorXd g;
    if (!at(q, v, g) || !(v < vSaddle)) {
      continue;
    }
    sideQ[side].push_back(q);
    sideV[side].push_back(v);
    double alpha = stiff > 0.0 ? 1.0 / stiff : 1.0;
    for (long k = 0; k < maxSteps && g.norm() > gradTol; ++k) {
      bool lowered = false;
      for (int halving = 0; halving < 30 && !lowered; ++halving) {
        const VectorXd trial = q + capped(-alpha * g);
        double vt = 0.0;
        VectorXd gt;
        if (at(trial, vt, gt) && vt < v) {
          q = trial;
          v = vt;
          g = gt;
          lowered = true;
          alpha *= 1.2;
        } else {
          alpha *= 0.5;
        }
      }
      if (!lowered) {
        break;
      }
      sideQ[side].push_back(q);
      sideV[side].push_back(v);
    }
  }
  pathQ.assign(sideQ[0].rbegin(), sideQ[0].rend());
  pathV.assign(sideV[0].rbegin(), sideV[0].rend());
  pathQ.push_back(qSaddle);
  pathV.push_back(vSaddle);
  pathQ.insert(pathQ.end(), sideQ[1].begin(), sideQ[1].end());
  pathV.insert(pathV.end(), sideV[1].begin(), sideV[1].end());
  // The reactant sits at q = 0, and a path starts at the reactant end.
  if (pathQ.back().norm() < pathQ.front().norm()) {
    std::reverse(pathQ.begin(), pathQ.end());
    std::reverse(pathV.begin(), pathV.end());
  }
}

/// Mode rate: the ring-polymer instanton through the saddle out of the
/// reactant, and its thermal rate.
std::vector<std::string>
runRate(const Parameters &params, const std::shared_ptr<Potential> &pot,
        const Matter &reactant, const MassWeighted &mw,
        const tunneling::BatchPotential &evaluate,
        const std::function<MatrixXd(const VectorXd &)> &hessianAt,
        const std::function<MatrixXd(const VectorXd &)> &beadHessianAt,
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
  const VectorXd unstableMode = Eigen::SelfAdjointEigenSolver<MatrixXd>(
                                    0.5 * (hSaddle + hSaddle.transpose()))
                                    .eigenvectors()
                                    .col(0);
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
  EONC_LOG_INFO("[Instanton] rotation curvatures over the softest vibration "
                "{:.3g}, {:.3g}, {:.3g} (zero modes at or below {:.3g})",
                rotationResidual[0], rotationResidual[1], rotationResidual[2],
                tunneling::kRotationZeroFraction);

  std::vector<std::pair<std::string, double>> extras{
      {"instanton_crossover_K", tc},
      {"barrier_classical", vSaddle - vReactant}};
  auto write = [&](RunStatus status) {
    auto env = JobResultEnvelope::fromMinimization(
        status, params.potential_options().potential,
        PotRegistry::get().total_force_calls(), false, 0.0);
    env.provenance = provenanceForJob(params);
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
  long sdPathForceCalls = 0;
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
      const long before = PotRegistry::get().total_force_calls();
      std::vector<VectorXd> grads;
      evaluate(pathQ, pathV, grads);
      sdPathForceCalls = PotRegistry::get().total_force_calls() - before;
    } else {
      sdPathForceCalls = 0;
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
  // No band: the steepest-descent path out of the saddle carries the ring
  // seed by the period condition. Cooling a cosine from the crossover finds
  // the ring only where it grows continuously out of the saddle; where it
  // does not, the search walks to a neighbouring saddle.
  if (!profile) {
    pathQ.clear();
    pathV.clear();
    const long before = PotRegistry::get().total_force_calls();
    steepestDescentPath(qSaddle, vSaddle, hSaddle, mw.sqrtMasses(), evaluate,
                        pathQ, pathV);
    sdPathForceCalls = PotRegistry::get().total_force_calls() - before;
    std::vector<double> arc(pathQ.size(), 0.0);
    for (size_t k = 1; k < pathQ.size(); ++k) {
      arc[k] = arc[k - 1] + (pathQ[k] - pathQ[k - 1]).norm();
    }
    if (pathQ.size() >= 3 && pathQ.size() == pathV.size()) {
      std::vector<std::shared_ptr<Matter>> images;
      std::vector<io::ConFrameMetadata> metas;
      images.reserve(pathQ.size());
      for (size_t k = 0; k < pathQ.size(); ++k) {
        auto image = std::make_shared<Matter>(reactant);
        mw.place(pathQ[k], *image);
        io::ConFrameMetadata meta;
        meta.frame_index = static_cast<uint64_t>(k);
        meta.energy = pathV[k];
        meta.scalars.push_back({"arc_length", arc[k]});
        metas.push_back(std::move(meta));
        images.push_back(std::move(image));
      }
      const auto frames = io::buildNebPathFrames(images, metas);
      if (frames.empty() ||
          !io::io_ok(io::writeConFrames("instanton_sd_path.con", frames))) {
        throw std::runtime_error(
            "instanton: cannot write instanton_sd_path.con");
      }
      returnFiles.push_back("instanton_sd_path.con");
    }
    try {
      profile = std::make_unique<tunneling::Profile>(std::move(arc), pathV);
      hwPath = tunneling::hbarOmega(tunneling::wellCurvature(*profile, true));
      EONC_LOG_INFO("[Instanton] steepest-descent path of {} points, {} "
                    "force calls",
                    pathQ.size(),
                    PotRegistry::get().total_force_calls() - before);
      for (size_t k = 0; k < pathQ.size(); ++k) {
        EONC_LOG_DEBUG("[Instanton] path {} s {:.6f} V - V_reactant {:.8f}", k,
                       profile->s()[k], pathV[k] - vReactant);
      }
    } catch (const std::exception &ex) {
      EONC_LOG_WARNING("[Instanton] steepest-descent path unusable: {}; "
                       "seeding from the saddle mode",
                       ex.what());
      if (!profile) {
        pathQ.clear();
        pathV.clear();
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
      // The ring is not searched. Above T_c the rate is the parabolic
      // barrier factor times quantum harmonic TST. At T_c the factor
      // diverges and that temperature records no rate.
      bool wrote = false;
      if (temperature > tc) {
        try {
          const double factor = tunneling::parabolicFactor(temperature, tc);
          const double logQhtst = tunneling::quantumHarmonicTstLogRate(
              hReactant, hSaddle, beta, vSaddle - vReactant, rigidModes);
          const double logPar = logQhtst + std::log(factor);
          const double kPar = std::exp(logPar) / tunneling::kTimeUnitSeconds;
          EONC_LOG_INFO("[Instanton] {:.4g} K is above the crossover {:.4g} "
                        "K; parabolic factor {:.6g}, ln(k s) = {:.4f}",
                        temperature, tc, factor, logPar - logSecond);
          table << temperature << ' ' << tc << ' ' << o.beads
                << " 0 0 nan 0 nan nan nan nan " << wkbLog << ' '
                << (logPar - logSecond) << ' ' << factor << '\n';
          if (last) {
            extras.emplace_back("instanton_temperature_K", temperature);
            extras.emplace_back("parabolic_factor", factor);
            extras.emplace_back("rate_parabolic", kPar);
            extras.emplace_back("rate_parabolic_log", logPar - logSecond);
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
        EONC_LOG_ERROR("[Instanton] {:.4g} K is at the crossover temperature "
                       "{:.4g} K; the parabolic factor diverges there",
                       temperature, tc);
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
    ro.initialHessians = o.initial_hessians;
    ro.energyShift = o.energy_shift;
    ro.discretization = o.discretization;
    ro.friction = o.friction != "none";
    ro.frictionExplicit = o.friction == "explicit";
    ro.frictionEta = o.friction_eta;
    ro.frictionEtaBeads = o.friction_eta_beads;
    if (static_cast<long>(mw.sqrtMasses().size()) == reactant.numberOfAtoms()) {
      ro.rigidSqrtMasses = mw.sqrtMasses();
      ro.rigidReference = mw.referenceFree();
      ro.rigidRotations = rotationZero;
    }
    std::vector<VectorXd> guess = ring;
    if (guess.empty() && profile) {
      try {
        tunneling::RingSeed seed;
        guess = tunneling::ringFromPath(pathQ, pathV, beta * tunneling::kHbar,
                                        o.beads, &seed);
        EONC_LOG_INFO("[Instanton] seed orbit {:.6f} eV above the reactant, "
                      "period {:.4g} against beta hbar {:.4g}; the path's "
                      "higher end sits {:.6f} eV above the reactant",
                      seed.energy - vReactant, seed.period,
                      beta * tunneling::kHbar, seed.pathLow - vReactant);
        EONC_LOG_INFO("[Instanton] ring seeded from {} by the period condition",
                      o.initial_path.empty() ? "the steepest-descent path"
                                             : o.initial_path);
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

    if (inst.collapsed) {
      EONC_LOG_ERROR("[Instanton] {:.4g} K: the ring collapsed after {} "
                     "iterations (B_N {:.3e}, every bead at one point near "
                     "s = {:.4f}); the search left the bounce for a "
                     "stationary point of V, no rate",
                     temperature, inst.iterations, inst.bN,
                     inst.beads.empty()
                         ? 0.0
                         : (inst.beads.front() - qSaddle).dot(unstableMode));
    }
    const tunneling::RingChannel channel =
        tunneling::ringChannel(inst.beads, qSaddle, unstableMode);
    EONC_LOG_INFO("[Instanton] {:.4g} K: ring spans s = {:.4f} to {:.4f} "
                  "amu^0.5 A along the unstable mode, chord overlap {:.3f}, "
                  "dividing-plane crossing {:.4f} amu^0.5 A off the saddle",
                  temperature, channel.sMin, channel.sMax, channel.chordOverlap,
                  channel.crossingOffset);
    if (inst.converged && !channel.belongs) {
      EONC_LOG_ERROR("[Instanton] {:.4g} K: the ring does not pass through "
                     "the seeded saddle's channel (it must straddle the "
                     "dividing plane, cross it within its own span of the "
                     "saddle, and run within 60 degrees of the unstable "
                     "mode); it belongs to another saddle, no rate",
                     temperature);
      inst.converged = false;
    }
    bool rateOk = false;
    if (inst.converged) {
      // Bead Hessians on every stride-th bead of the ring, linear in
      // between, wrapping from the last anchor back to bead 0.
      const long stride = std::max<long>(1, o.hessian_stride);
      const long nBeads = o.beads;
      std::map<long, MatrixXd> anchors;
      // Bead N - j mirrors bead j on an out-and-back ring and shares its
      // Hessian.
      auto anchor = [&](long j) -> const MatrixXd & {
        const long n = static_cast<long>(inst.beads.size());
        if (j > n / 2 && (inst.beads[static_cast<size_t>(j)] -
                          inst.beads[static_cast<size_t>(n - j)])
                                 .norm() <= 1e-10) {
          j = n - j;
        }
        auto it = anchors.find(j);
        if (it == anchors.end()) {
          it =
              anchors
                  .emplace(j, beadHessianAt(inst.beads[static_cast<size_t>(j)]))
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
        tunneling::RingRigidBodies bodies;
        bodies.sqrtMasses = ro.rigidSqrtMasses;
        bodies.reference = ro.rigidReference;
        bodies.rotations = ro.rigidRotations;
        tunneling::instantonRate(
            inst, beadHessian, hReactant, vReactant - o.energy_shift, hSaddle,
            vSaddle - o.energy_shift, rigidModes, 4096, bodies);
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
    {
      const std::string centroidFile =
          last ? "instanton_centroid.con"
               : "instanton_centroid_" +
                     files.back().substr(std::string("instanton_").size());
      writeCentroid(centroidFile, inst.beads, mw, reactant,
                    {{"instanton_temperature_K", temperature},
                     {"instanton_crossover_K", tc},
                     {"instanton_converged", inst.converged ? 1.0 : 0.0}});
      returnFiles.push_back(centroidFile);
    }

    ring = inst.beads;
    if (!last) {
      continue;
    }
    extras.emplace_back("instanton_temperature_K", temperature);
    extras.emplace_back("sd_path_force_calls",
                        static_cast<double>(sdPathForceCalls));
    extras.emplace_back("instanton_iterations",
                        static_cast<double>(inst.iterations));
    extras.emplace_back("instanton_ring_potential", inst.ringPotential);
    extras.emplace_back("instanton_bN", inst.bN);
    extras.emplace_back("instanton_collapsed", inst.collapsed ? 1.0 : 0.0);
    extras.emplace_back("instanton_s_min", channel.sMin);
    extras.emplace_back("instanton_s_max", channel.sMax);
    extras.emplace_back("instanton_chord_overlap", channel.chordOverlap);
    extras.emplace_back("instanton_crossing_offset", channel.crossingOffset);
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
  if (o.pi_planes > 0) {
    const auto planeFiles = piqtst::runAfterInstanton(
        params, *pot, reactant, saddle, hSaddle, pathQ, temperatures, extras);
    returnFiles.insert(returnFiles.end(), planeFiles.begin(), planeFiles.end());
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
  auto hessianWith = [&](const VectorXd &q, bool projectRotations) {
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
    // rigid eigenvalues of either sign; project them to zero. At a point
    // that is not stationary, H r for a rotation r is the rotated gradient,
    // which on a ring bead balances the springs; ring beads keep it and
    // lose only their translations.
    const std::array<bool, 3> none{{false, false, false}};
    const MatrixXd rigid =
        mw.rigidBasis(m, projectRotations ? rotationZero : none);
    if (rigid.cols() > 0) {
      const MatrixXd p = MatrixXd::Identity(n, n) - rigid * rigid.transpose();
      out = p * out * p;
    }
    return out;
  };
  auto hessianAt = [&](const VectorXd &q) { return hessianWith(q, true); };
  auto beadHessianAt = [&](const VectorXd &q) {
    if (!rotationsKnown) {
      hessianWith(VectorXd::Zero(n), true);
    }
    return hessianWith(q, false);
  };

  if (o.mode == "rate") {
    return runRate(params, pot, *reactant, mw, evaluate, hessianAt,
                   beadHessianAt, rotationZero, rotationResidual);
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
  writeCentroid("instanton_centroid.con", inst.path, mw, *reactant,
                {{"instanton_temperature_K", kelvin},
                 {"instanton_converged", inst.converged ? 1.0 : 0.0}});
  returnFiles.push_back("instanton_centroid.con");

  // A converged path between wells too far apart is a result, not a
  // failure: the flags say why no splitting was written.
  const bool good = splitOk || (inst.converged && !inst.symmetricEnough);
  const auto status = good ? RunStatus::GOOD
                           : (inst.converged ? RunStatus::FAIL_POTENTIAL_FAILED
                                             : RunStatus::FAIL_MAX_ITERATIONS);
  auto env = JobResultEnvelope::fromMinimization(
      status, params.potential_options().potential,
      PotRegistry::get().total_force_calls(), false, 0.0);
  env.provenance = provenanceForJob(params);
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

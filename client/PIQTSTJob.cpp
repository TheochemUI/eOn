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
#include "eon/ConFileIO.h"
#include "eon/EonLogger.h"
#include "eon/Matter.h"
#include "eon/PIQTST.h"
#include "eon/Parameters.h"
#include "eon/Tunneling.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace eonc::piqtst {

void validateOptions(const instanton_options_t &o) {
  if (o.pi_planes < 0) {
    throw std::invalid_argument("[Instanton] pi_planes must not be negative");
  }
  if (o.pi_planes == 0) {
    return;
  }
  if (o.mode != "rate") {
    throw std::invalid_argument("[Instanton] pi_planes needs mode = rate");
  }
  if (o.pi_planes < 2) {
    throw std::invalid_argument(
        "[Instanton] pi_planes must be 0 or at least 2");
  }
  if (o.pi_beads < 1) {
    throw std::invalid_argument("[Instanton] pi_beads must be positive");
  }
  if (o.pi_equilibration_steps < 0 || o.pi_sampling_steps < 20) {
    throw std::invalid_argument(
        "[Instanton] pi_equilibration_steps must not be negative and "
        "pi_sampling_steps must be at least 20 (ten blocks of two)");
  }
  if (!(o.pi_time_step > 0.0) || !(o.pi_pile_tau > 0.0) ||
      !(o.pi_pile_scale > 0.0)) {
    throw std::invalid_argument("[Instanton] pi_time_step, pi_pile_tau and "
                                "pi_pile_scale must be positive");
  }
  if (o.pi_seed < 0) {
    throw std::invalid_argument("[Instanton] pi_seed must not be negative");
  }
  if (o.pi_thermostat != "pile" && o.pi_thermostat != "piglet") {
    throw std::invalid_argument(
        "[Instanton] pi_thermostat must be pile or piglet, not " +
        o.pi_thermostat);
  }
  if (o.pi_thermostat == "piglet" && o.pi_beads > 1 && o.pi_gle_file.empty()) {
    throw std::invalid_argument(
        "[Instanton] pi_thermostat = piglet needs pi_gle_file");
  }
  if (o.pi_direction != "mode" && o.pi_direction != "line") {
    throw std::invalid_argument(
        "[Instanton] pi_direction must be mode or line, not " + o.pi_direction);
  }
  if (!(o.pi_reactant_extent >= 0.0)) {
    throw std::invalid_argument(
        "[Instanton] pi_reactant_extent must not be negative");
  }
  if (o.pi_recrossing_parents < 0 || o.pi_recrossing_parents == 1) {
    throw std::invalid_argument(
        "[Instanton] pi_recrossing_parents must be 0 or at least 2");
  }
  if (o.pi_recrossing_parents == 0) {
    return;
  }
  if (o.pi_recrossing_children < 1 || o.pi_recrossing_spacing < 1) {
    throw std::invalid_argument("[Instanton] pi_recrossing_children and "
                                "pi_recrossing_spacing must be positive");
  }
  if (!(o.pi_recrossing_time >= 4.0 * o.pi_time_step)) {
    throw std::invalid_argument(
        "[Instanton] pi_recrossing_time must be at least four pi_time_step");
  }
}

namespace {

/// The instanton's mass-weighted coordinates: free atoms in index order,
/// three components each, measured from the reactant.
struct Embedding {
  std::vector<long> freeAtoms;
  std::vector<double> sqrtMass;
  long atoms{0};

  explicit Embedding(const Matter &reactant)
      : atoms(reactant.numberOfAtoms()) {
    for (long i = 0; i < atoms; ++i) {
      if (!reactant.getFixed(i)) {
        freeAtoms.push_back(i);
        sqrtMass.push_back(std::sqrt(reactant.getMass(i)));
      }
    }
  }
  long dimension() const { return 3 * static_cast<long>(freeAtoms.size()); }
  /// A mass-weighted vector over the free atoms as one over every atom.
  VectorXd full(const VectorXd &q) const {
    VectorXd out = VectorXd::Zero(3 * atoms);
    for (size_t k = 0; k < freeAtoms.size(); ++k) {
      for (int c = 0; c < 3; ++c) {
        out(3 * freeAtoms[k] + c) = q(static_cast<long>(3 * k) + c);
      }
    }
    return out;
  }
  /// Cartesian positions, 3 * atoms, at mass-weighted displacement q.
  VectorXd cartesian(const VectorXd &reference, const VectorXd &q) const {
    VectorXd x = reference;
    for (size_t k = 0; k < freeAtoms.size(); ++k) {
      for (int c = 0; c < 3; ++c) {
        x(3 * freeAtoms[k] + c) +=
            q(static_cast<long>(3 * k) + c) / sqrtMass[k];
      }
    }
    return x;
  }
  VectorXd toQ(const Matter &reactant, const Matter &m) const {
    const AtomMatrix d =
        reactant.pbc(m.getPositions() - reactant.getPositions());
    VectorXd q(dimension());
    for (size_t k = 0; k < freeAtoms.size(); ++k) {
      for (int c = 0; c < 3; ++c) {
        q(static_cast<long>(3 * k) + c) = sqrtMass[k] * d(freeAtoms[k], c);
      }
    }
    return q;
  }
};

VectorXd flatten(const AtomMatrix &m) {
  VectorXd x(3 * m.rows());
  for (long i = 0; i < m.rows(); ++i) {
    for (int c = 0; c < 3; ++c) {
      x(3 * i + c) = m(i, c);
    }
  }
  return x;
}

std::string kelvinTag(double t) {
  std::ostringstream name;
  name << std::defaultfloat << std::setprecision(6) << t << "K";
  return name.str();
}

} // namespace

std::vector<std::string>
runAfterInstanton(const Parameters &params, Potential &pot,
                  const Matter &reactant, const Matter &saddle,
                  const MatrixXd &hSaddle, const std::vector<VectorXd> &pathQ,
                  const std::vector<double> &temperatures,
                  std::vector<std::pair<std::string, double>> &extras) {
  const auto &o = params.instanton_options();
  validateOptions(o);
  const double timeUnit = params.constants().timeUnit;
  const Embedding emb(reactant);
  const VectorXd qSaddle = emb.toQ(reactant, saddle);
  const double qNorm = qSaddle.norm();
  if (!(qNorm > 0.0)) {
    throw std::runtime_error("piqtst: the saddle sits on the reactant");
  }
  const VectorXd lineDir = qSaddle / qNorm;

  // The plane normal. The saddle's unstable mode is the dividing surface
  // the rate is evaluated on; one normal for every plane keeps s a linear
  // coordinate, so the mean force integrates to its free energy.
  VectorXd dir = lineDir;
  std::string chosen = "line";
  if (o.pi_direction == "mode") {
    const Eigen::SelfAdjointEigenSolver<MatrixXd> es(
        0.5 * (hSaddle + hSaddle.transpose()));
    if (es.info() == Eigen::Success && es.eigenvalues()(0) < 0.0) {
      VectorXd v = es.eigenvectors().col(0);
      if (v.dot(lineDir) < 0.0) {
        v = -v;
      }
      if (v.dot(lineDir) >= 0.5) {
        dir = v;
        chosen = "mode";
      } else {
        EONC_LOG_WARNING("[Instanton] piqtst: the unstable mode makes {:.3f} "
                         "rad with the reactant-saddle line; planes follow "
                         "the line",
                         std::acos(v.dot(lineDir)));
      }
    } else {
      EONC_LOG_WARNING("[Instanton] piqtst: the saddle Hessian has no "
                       "negative eigenvalue; planes follow the line");
    }
  }
  if (emb.freeAtoms.size() == static_cast<size_t>(reactant.numberOfAtoms()) &&
      !reactant.getPeriodic()) {
    EONC_LOG_WARNING("[Instanton] piqtst: no atom is fixed; the planes are "
                     "fixed in the reactant's frame and a rotation of the "
                     "cluster moves the centroid between them");
  }
  const double sStar = dir.dot(qSaddle);
  const long nPlanes = o.pi_planes;
  const double s0 = -o.pi_reactant_extent * sStar;
  std::vector<double> sPlanes;
  for (long j = 0; j < nPlanes; ++j) {
    sPlanes.push_back(s0 + (sStar - s0) * static_cast<double>(j) /
                               static_cast<double>(nPlanes - 1));
  }

  Coordinate c;
  c.atoms = reactant.numberOfAtoms();
  c.masses.resize(static_cast<size_t>(c.atoms));
  c.numbers.resize(static_cast<size_t>(c.atoms));
  c.free.assign(static_cast<size_t>(3 * c.atoms), 0);
  const VectorXi z = reactant.getAtomicNrs();
  for (long i = 0; i < c.atoms; ++i) {
    c.masses[static_cast<size_t>(i)] = reactant.getMass(i);
    c.numbers[static_cast<size_t>(i)] = z(i);
  }
  for (const long i : emb.freeAtoms) {
    for (int a = 0; a < 3; ++a) {
      c.free[static_cast<size_t>(3 * i + a)] = 1;
    }
  }
  c.reference = flatten(reactant.getPositions());
  c.direction = emb.full(dir);
  const Matrix3d cell =
      reactant.getPeriodic() ? reactant.getCell() : Matrix3d::Zero().eval();
  c.box = cell.data();

  // Seeds: the band where it crosses the plane, else the straight line
  // through the reactant and the saddle.
  std::vector<double> bandS;
  for (const auto &q : pathQ) {
    bandS.push_back(dir.dot(q));
  }
  auto seed = [&](double s) -> VectorXd {
    for (size_t k = 0; k + 1 < pathQ.size(); ++k) {
      const double lo = bandS[k];
      const double hi = bandS[k + 1];
      if ((lo <= s && s <= hi) || (hi <= s && s <= lo)) {
        const double t = hi != lo ? (s - lo) / (hi - lo) : 0.0;
        return emb.cartesian(c.reference,
                             (1.0 - t) * pathQ[k] + t * pathQ[k + 1]);
      }
    }
    return emb.cartesian(c.reference, (s / sStar) * qSaddle);
  };

  const double vReactant = Matter(reactant).getPotentialEnergy();
  const double vSaddle = Matter(saddle).getPotentialEnergy();
  double effective = std::numeric_limits<double>::quiet_NaN();
  for (const auto &kv : extras) {
    if (kv.first == "barrier_effective_instanton") {
      effective = kv.second;
    }
  }
  EONC_LOG_INFO("[Instanton] piqtst: {} planes from s = {:.4f} to s* = {:.4f} "
                "amu^0.5 A along the {}, {} beads, {} + {} steps of {:.4g} fs",
                nPlanes, s0, sStar, chosen == "mode" ? "unstable mode" : "line",
                o.pi_beads, o.pi_equilibration_steps, o.pi_sampling_steps,
                o.pi_time_step);

  const std::string tableFile = "rate_piqtst.dat";
  std::ofstream table(tableFile);
  if (!table) {
    throw std::runtime_error("piqtst: cannot write " + tableFile);
  }
  const bool kappaOn = o.pi_recrossing_parents > 0;
  if (kappaOn && o.pi_thermostat == "piglet") {
    EONC_LOG_WARNING("[Instanton] piqtst: the transmission factor's parents "
                     "take pile; RPMD needs the ring polymer's own "
                     "distribution, which piglet does not sample");
  }
  table << "# T_K s_amu05A dF_ds_eV_per_amu05A dF_ds_error F_eV F_error_eV "
           "spread_max_A";
  if (kappaOn) {
    table << " kappa kappa_error ln_k_rpmd_s ln_k_rpmd_s_error";
  }
  table << '\n';
  table << std::setprecision(10);
  std::vector<std::string> files{tableFile};

  std::vector<double> sorted = temperatures;
  std::sort(sorted.begin(), sorted.end(), std::greater<>());
  const double logSecond = std::log(tunneling::kTimeUnitSeconds);
  for (size_t ti = 0; ti < sorted.size(); ++ti) {
    const double t = sorted[ti];
    const bool last = ti + 1 == sorted.size();
    const double beta = 1.0 / (tunneling::kBoltzmann * t);
    ScanOptions so;
    so.planes = sPlanes;
    so.equilibration = o.pi_equilibration_steps;
    so.production = o.pi_sampling_steps;
    so.blocks = 10;
    so.ring.beads = o.pi_beads;
    so.ring.temperature = t;
    so.ring.kB = tunneling::kBoltzmann;
    so.ring.hbar = tunneling::kHbar;
    so.ring.dt = o.pi_time_step / timeUnit;
    so.ring.thermostat = o.pi_thermostat == "piglet"
                             ? pathintegral::Thermostat::Piglet
                             : pathintegral::Thermostat::Pile;
    so.ring.pileTau = o.pi_pile_tau / timeUnit;
    so.ring.pileScale = o.pi_pile_scale;
    so.ring.gleFile = o.pi_gle_file;
    so.ring.seed = static_cast<std::uint64_t>(o.pi_seed);
    so.seed = seed;
    const std::vector<Plane> planes = scan(pot, c, so);
    const Rate r = rate(planes, beta);
    const double lnPerSecond = r.logRate - logSecond;

    Recrossing kappa;
    double lnRpmd = std::numeric_limits<double>::quiet_NaN();
    double lnRpmdError = std::numeric_limits<double>::quiet_NaN();
    if (kappaOn) {
      RecrossingOptions ro;
      ro.s = planes.back().s;
      ro.equilibration = o.pi_equilibration_steps;
      ro.parents = o.pi_recrossing_parents;
      ro.spacing = o.pi_recrossing_spacing;
      ro.children = o.pi_recrossing_children;
      ro.steps = std::lround(o.pi_recrossing_time / o.pi_time_step);
      ro.ring = so.ring;
      ro.seed = seed;
      kappa = recrossing(pot, c, ro);
      if (kappa.plateau > 0.0) {
        lnRpmd = lnPerSecond + std::log(kappa.plateau);
        lnRpmdError =
            std::hypot(r.logRateError, kappa.plateauError / kappa.plateau);
      } else {
        EONC_LOG_WARNING("[Instanton] piqtst {:.4g} K: the transmission "
                         "factor {:.4g} +- {:.2g} is not positive; ln k_RPMD "
                         "is undefined (raise pi_recrossing_parents)",
                         t, kappa.plateau, kappa.plateauError);
      }
      EONC_LOG_INFO("[Instanton] piqtst {:.4g} K: transmission factor "
                    "kappa = {:.4f} +- {:.4f} from {} trajectories, "
                    "ln(k_RPMD s) = {:.4f} +- {:.4f}",
                    t, kappa.plateau, kappa.plateauError, kappa.trajectories,
                    lnRpmd, lnRpmdError);
      std::vector<std::string> curveFiles;
      if (last) {
        curveFiles.push_back("kappa_piqtst.dat");
      }
      if (sorted.size() > 1) {
        curveFiles.push_back("kappa_piqtst_" + kelvinTag(t) + ".dat");
      }
      for (const auto &file : curveFiles) {
        std::ofstream curve(file);
        if (!curve) {
          throw std::runtime_error("piqtst: cannot write " + file);
        }
        curve << "# t_fs kappa\n" << std::setprecision(10);
        for (size_t i = 0; i < kappa.time.size(); ++i) {
          curve << kappa.time[i] * timeUnit << ' ' << kappa.kappa[i] << '\n';
        }
        files.push_back(file);
      }
    }

    std::vector<std::string> conFiles;
    if (last) {
      conFiles.push_back("piqtst_planes.con");
    }
    if (sorted.size() > 1) {
      conFiles.push_back("piqtst_planes_" + kelvinTag(t) + ".con");
    }
    for (size_t j = 0; j < planes.size(); ++j) {
      const Plane &p = planes[j];
      double largest = 0.0;
      for (const double sp : p.spread) {
        largest = std::max(largest, sp);
      }
      table << t << ' ' << p.s << ' ' << p.meanForce << ' ' << p.meanForceError
            << ' ' << p.freeEnergy << ' ' << p.freeEnergyError << ' '
            << largest;
      if (kappaOn) {
        table << ' ' << kappa.plateau << ' ' << kappa.plateauError << ' '
              << lnRpmd << ' ' << lnRpmdError;
      }
      table << '\n';
      Matter frame(reactant);
      AtomMatrix pos = frame.getPositions();
      for (long i = 0; i < c.atoms; ++i) {
        for (int a = 0; a < 3; ++a) {
          pos(i, a) = p.centroid(3 * i + a);
        }
      }
      frame.setPositions(pos);
      io::ConFrameMetadata meta;
      meta.frame_index = static_cast<uint64_t>(j);
      meta.write_con_forces = false;
      meta.spreads = p.spread;
      meta.scalars = {{"piqtst_temperature_K", t},
                      {"piqtst_s", p.s},
                      {"piqtst_dF_ds", p.meanForce},
                      {"piqtst_dF_ds_error", p.meanForceError},
                      {"piqtst_F", p.freeEnergy},
                      {"piqtst_F_error", p.freeEnergyError},
                      {"beads", static_cast<double>(o.pi_beads)},
                      {"spread_max", largest}};
      for (const auto &file : conFiles) {
        if (!io::io_ok(frame.matter2con(file, j > 0, &meta))) {
          throw std::runtime_error("piqtst: cannot write " + file);
        }
      }
    }
    for (const auto &file : conFiles) {
      files.push_back(file);
    }

    const Plane &top = planes.back();
    EONC_LOG_INFO("[Instanton] piqtst {:.4g} K: free-energy barrier {:.5f} "
                  "+- {:.5f} eV (classical {:.5f} eV, instanton effective "
                  "{:.5f} eV), ln(k s) = {:.4f} +- {:.4f}",
                  t, r.barrier, r.barrierError, vSaddle - vReactant, effective,
                  lnPerSecond, r.logRateError);
    if (!(r.barrier > 0.0)) {
      EONC_LOG_WARNING("[Instanton] piqtst {:.4g} K: F(s*) is {:.4g} eV "
                       "against the lowest plane before it; s* is no barrier "
                       "on the centroid free energy",
                       t, r.barrier);
    }
    if (r.firstPlaneHeight < 5.0) {
      EONC_LOG_WARNING("[Instanton] piqtst {:.4g} K: the first plane is "
                       "{:.3g} kT above the reactant minimum of F; the "
                       "reactant integral is cut there (raise "
                       "pi_reactant_extent)",
                       t, r.firstPlaneHeight);
    }
    if (std::abs(top.meanForce) > 3.0 * top.meanForceError) {
      EONC_LOG_WARNING("[Instanton] piqtst {:.4g} K: dF/ds at the saddle "
                       "plane is {:.4g} +- {:.2g} eV / (amu^0.5 A); the "
                       "maximum of F is not at s*",
                       t, top.meanForce, top.meanForceError);
    }
    if (!last) {
      continue;
    }
    extras.emplace_back("piqtst_temperature_K", t);
    extras.emplace_back("piqtst_planes", static_cast<double>(nPlanes));
    extras.emplace_back("piqtst_beads", static_cast<double>(o.pi_beads));
    extras.emplace_back("piqtst_s_star", sStar);
    extras.emplace_back("piqtst_dF_ds_star", top.meanForce);
    extras.emplace_back("piqtst_dF_ds_star_error", top.meanForceError);
    extras.emplace_back("piqtst_first_plane_kT", r.firstPlaneHeight);
    extras.emplace_back("barrier_piqtst", r.barrier);
    extras.emplace_back("barrier_piqtst_error", r.barrierError);
    extras.emplace_back("rate_piqtst", std::exp(lnPerSecond));
    extras.emplace_back("rate_piqtst_log", lnPerSecond);
    extras.emplace_back("rate_piqtst_log_error", r.logRateError);
    if (kappaOn) {
      extras.emplace_back("piqtst_kappa", kappa.plateau);
      extras.emplace_back("piqtst_kappa_error", kappa.plateauError);
      extras.emplace_back("rate_rpmd", std::exp(lnRpmd));
      extras.emplace_back("rate_rpmd_log", lnRpmd);
      extras.emplace_back("rate_rpmd_log_error", lnRpmdError);
    }
  }
  return files;
}

} // namespace eonc::piqtst

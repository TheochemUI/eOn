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

#include "eon/EonLogger.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include "eon/BaseStructures.h"
#include "eon/BasinHoppingJob.h"
#include "eon/Dynamics.h"
#include "eon/HelperFunctions.h"
#include "eon/JobResult.h"
#include "eon/LBFGS.h"
#include "eon/ObjectiveFunction.h"
#include "eon/Optimizer.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"

namespace eonc {

std::vector<std::string> BasinHoppingJob::run() {
  bool swapMove;
  double swap_accept = 0.0;
  jump_count = 0; // count of jump movies
  swap_count = 0; // count of swap moves
  disp_count = 0; // count of displacement moves
  // Quench tail may not run when stop_energy breaks first.
  int quench_displacements = 0;
  int consecutive_rejected_trials = 0;
  double totalAccept = 0.0;
  std::unique_ptr<Matter> minTrial = std::make_unique<Matter>(pot, params);
  std::unique_ptr<Matter> swapTrial = std::make_unique<Matter>(pot, params);

  std::string conFilename =
      eonc::helpers::getRelevantFile(params.main_options().conFilename);
  if (!eonc::io::io_ok(current->con2matter(conFilename))) {
    QUILL_LOG_CRITICAL(log, "Failed to load {}", conFilename);
    throw std::runtime_error("failed to load " + conFilename);
  }

  // Sanity Check
  std::vector<long> Elements;
  Elements = getElements(current.get());
  if (params.basin_hopping_options().swap_probability > 0 &&
      Elements.size() == 1) {
    log = eonc::log::traceback();
    QUILL_LOG_CRITICAL(log,
                       "error: [Basin Hopping] swap move probability must be "
                       "zero if there is only one element type\n");
    throw std::invalid_argument(
        "[Basin Hopping] swap_probability must be zero with one element type");
  }

  double randomProb =
      params.basin_hopping_options().initial_random_structure_probability;
  if (randomProb > 0.0) {
    QUILL_LOG_DEBUG(log, "generating random structure with probability {:.4f}",
                    randomProb);
  }
  double u = eonc::rng::random();
  if (u < params.basin_hopping_options().initial_random_structure_probability) {
    AtomMatrix randomPositions = current->getPositionsFree();
    for (int i = 0; i < current->numberOfFreeAtoms(); i++) {
      for (int j = 0; j < 3; j++) {
        randomPositions(i, j) = eonc::rng::random();
      }
    }
    randomPositions *= current->getCell();
    current->setPositionsFree(randomPositions);

    eonc::geometry::pushApart(
        current, params.basin_hopping_options().push_apart_distance);
  }

  *trial = *current;
  *minTrial = *current;

  stock_ = configureStockQuench();
  temperature_ = params.main_options().temperature;
  initialDisplacement_ = params.basin_hopping_options().displacement;
  angularDisplacement_ = initialAngular_;
  if (stock_) {
    prepareRigidWater();
    if (!rigidWater_ && !current->getPeriodic()) {
      bool fixed = false;
      for (long i = 0; i < current->numberOfAtoms(); ++i) {
        if (current->getFixed(i)) {
          fixed = true;
          break;
        }
      }
      if (!fixed) {
        ParametersLoadAccess::optimizer_options(params).lbfgs.project_rigid =
            true;
      }
    }
    quench(*current, true, std::numeric_limits<double>::infinity());
  } else {
    current->relax(true);
  }

  double currentEnergy = current->getPotentialEnergy();
  double minimumEnergy = currentEnergy;

  auto minimumEnergyStructure = std::make_shared<Matter>(pot, params);
  *minimumEnergyStructure = *current;
  int nsteps = params.basin_hopping_options().steps +
               params.basin_hopping_options().quenching_steps;

  QUILL_LOG_DEBUG(
      log, "[Basin Hopping] {:4s} {:12s} {:12s} {:12s} {:4s} {:5s} {:5s}",
      "step", "current", "trial", "global min", "fc", "ar", "md");
  QUILL_LOG_DEBUG(
      log, "[Basin Hopping] {:4s} {:12s} {:12s} {:12s} {:4s} {:5s} {:5s}",
      "----", "-------", "-----", "----------", "--", "--", "--");

  int recentAccept = 0;
  double curDisplacement = params.basin_hopping_options().displacement;

  for (int step = 0; step < nsteps; step++) {

    // Swap or displace
    const double swapP = swapProbability(Elements);
    if (eonc::rng::randomDouble(1.0) < swapP &&
        step < params.basin_hopping_options().steps) {
      *swapTrial = *current;
      randomSwap(swapTrial.get());
      swapMove = true;
      *minTrial = *swapTrial;
    } else if (rigidWater_) {
      if (step >= params.basin_hopping_options().steps) {
        quench_displacements++;
      }
      *trial = *current;
      displaceRigid(trial.get(), curDisplacement, angularDisplacement_);
      swapMove = false;
      *minTrial = *trial;
    } else {
      AtomMatrix displacement;
      displacement = displaceRandom(curDisplacement);
      if (step >= params.basin_hopping_options().steps) {
        quench_displacements++;
      }

      trial->setPositions(current->getPositions() + displacement);
      swapMove = false;
      eonc::geometry::pushApart(
          trial, params.basin_hopping_options().push_apart_distance);

      *minTrial = *trial;
    }

    if (params.debug_options().write_movies) {
      if (!eonc::io::io_ok(trial->matter2con("trials", true))) {
        QUILL_LOG_WARNING(log, "Failed to append trials movie frame");
      }
    }

    const bool converged =
        quench(*minTrial, false, minimumEnergy);
    const double trialEnergy = minTrial->getPotentialEnergy();
    const bool finiteTrial = std::isfinite(trialEnergy);
    double deltaE = finiteTrial ? trialEnergy - currentEnergy : 0.0;
    double p = 0.0;
    if (converged && finiteTrial) {
      if (step >= params.basin_hopping_options().steps) {
        if (deltaE <= 0.0) {
          p = 1.0;
        }
      } else {
        // Divide only for an uphill hop at positive temperature.
        p = metropolisProbability(deltaE, params.constants().kB,
                                  hopTemperature());
      }
    }

    bool accepted = false;
    if (converged && finiteTrial && eonc::rng::randomDouble(1.0) < p) {
      accepted = true;
      if (params.basin_hopping_options().significant_structure) {
        *current = *minTrial;
      } else if (swapMove) {
        *current = *swapTrial;
      } else {
        *current = *trial;
      }
      if (swapMove) {
        swap_accept += 1;
      }
      if (step < params.basin_hopping_options().steps) {
        totalAccept += 1;
        recentAccept += 1;
      }

      currentEnergy = minTrial->getPotentialEnergy();

      if (currentEnergy < minimumEnergy) {
        minimumEnergy = currentEnergy;
        *minimumEnergyStructure = *minTrial;
        if (!eonc::io::io_ok(minimumEnergyStructure->matter2con("min.con"))) {
          QUILL_LOG_WARNING(log, "Failed to write min.con");
        }
      }

      if (params.basin_hopping_options().write_unique) {
        bool newStructure = true;
        for (unsigned int i = 0; i < uniqueEnergies.size(); i++) {
          // if minTrial has a different energy or a different structure
          // it is new, otherwise it is old
          if (std::fabs(currentEnergy - uniqueEnergies[i]) <
              params.structure_comparison_options().energy_difference) {
            Matter probe = *current;
            if (probe.compare(*uniqueStructures[i],
                              params.structure_comparison_options()
                                  .indistinguishable_atoms)) {
              newStructure = false;
            }
          }
        }

        if (newStructure) {
          uniqueEnergies.push_back(currentEnergy);
          auto currentCopy = std::make_shared<Matter>(pot, params);
          *currentCopy = *current;
          uniqueStructures.push_back(currentCopy);

          char fname[128];
          snprintf(fname, 128, "min_%.5i.con", step + 1);
          if (!eonc::io::io_ok(current->matter2con(fname))) {
            QUILL_LOG_WARNING(log, "Failed to write {}", fname);
          }
          returnFiles.push_back(fname);

          snprintf(fname, 128, "energy_%.5i.dat", step + 1);
          returnFiles.push_back(fname);
          {
            std::ofstream fh(fname);
            if (fh)
              fh << std::format("{:.10e}\n", currentEnergy);
          }
        }
      }

      consecutive_rejected_trials = 0; // STC: I think this should go here.
    } else {
      if (stock_ && converged && finiteTrial && deltaE > 0.0 &&
          step < params.basin_hopping_options().steps) {
        liftTemperature(deltaE);
      }
      consecutive_rejected_trials++;
    }

    if (params.debug_options().write_movies) {
      if (!eonc::io::io_ok(minTrial->matter2con("movie", true))) {
        QUILL_LOG_WARNING(log, "Failed to append basin-hopping movie frame");
      }
    }

    if (minimumEnergy < params.basin_hopping_options().stop_energy) {
      break;
    }

    if (consecutive_rejected_trials ==
            params.basin_hopping_options().jump_max &&
        step < params.basin_hopping_options().steps) {
      consecutive_rejected_trials = 0;
      AtomMatrix jump;
      const long nJump = jumpStepCount();
      for (int j = 0; j < nJump; j++) {
        jump_count++;
        if (rigidWater_) {
          displaceRigid(current.get(), curDisplacement, angularDisplacement_);
        } else {
          jump = displaceRandom(curDisplacement);
          current->setPositions(current->getPositions() + jump);
        }
        // Only a minimized jump is a basin. The raw geometry must not become
        // the Metropolis reference or the stored global minimum.
        if (params.basin_hopping_options().significant_structure) {
          if (!rigidWater_) {
            eonc::geometry::pushApart(
                current, params.basin_hopping_options().push_apart_distance);
          }
          if (stock_) {
            quench(*current, false, minimumEnergy);
          } else {
            current->relax(true);
          }
          currentEnergy = current->getPotentialEnergy();
          if (std::isfinite(currentEnergy) && currentEnergy < minimumEnergy) {
            minimumEnergy = currentEnergy;
            *minimumEnergyStructure = *current;
          }
        }
      }
    }

    int nadjust = params.basin_hopping_options().adjust_period;
    double adjustFraction = params.basin_hopping_options().adjust_fraction;
    // A zero period is not an interval; the modulo would divide by zero.
    if (nadjust != 0 && (step + 1) % nadjust == 0 &&
        params.basin_hopping_options().adjust_displacement) {
      if (stock_) {
        adaptStockStep(curDisplacement, recentAccept, nadjust);
      } else {
        double recentRatio =
            static_cast<double>(recentAccept) / static_cast<double>(nadjust);
        if (recentRatio > params.basin_hopping_options().target_ratio) {
          curDisplacement *= 1.0 + adjustFraction;
        } else {
          curDisplacement *= 1.0 - adjustFraction;
        }
      }

      recentAccept = 0;
    }
  }
  /* Save Results */

  std::string resultsFilename("results.dat");

  if (params.debug_options().write_movies) {
    std::string movieFilename("movie.con");
    returnFiles.push_back(movieFilename);
  }

  {
    auto env = JobResultEnvelope::fromMinimization(
        RunStatus::GOOD, params.potential_options().potential,
        pots().total_force_calls(), true, minimumEnergy);
    env.provenance = provenanceForJob(params);
    env.job_type = "basin_hopping";
    env.random_seed = params.main_options().randomSeed;
    env.extras.emplace_back("minimum_energy", minimumEnergy);
    const double nsteps_ratio = params.basin_hopping_options().steps;
    env.extras.emplace_back("acceptance_ratio",
                            nsteps_ratio ? totalAccept / nsteps_ratio : 0.0);
    if (params.basin_hopping_options().swap_probability > 0 || swap_count > 0) {
      env.extras.emplace_back(
          "swap_acceptance_ratio",
          swap_count ? swap_accept / static_cast<double>(swap_count) : 0.0);
    }
    env.extras.emplace_back(
        "total_normal_displacement_steps",
        static_cast<double>(disp_count - jump_count - quench_displacements));
    env.extras.emplace_back("total_jump_steps",
                            static_cast<double>(jump_count));
    env.extras.emplace_back("total_swap_steps",
                            static_cast<double>(swap_count));
    if (stock_) {
      env.extras.emplace_back("stock_temperature", temperature_);
      env.extras.emplace_back("stock_displacement", curDisplacement);
    }
    env.writeResultsDat(resultsFilename);
    returnFiles.push_back(resultsFilename);
  }

  std::string productFilename("min.con");
  if (eonc::io::io_ok(minimumEnergyStructure->matter2con(productFilename))) {
    returnFiles.push_back(productFilename);
  } else {
    QUILL_LOG_ERROR(log, "Failed to write {}", productFilename);
  }

  // minTrial and swapTrial automatically cleaned up by unique_ptr
  return returnFiles;
}

namespace {

class RigidBodyObjective : public ObjectiveFunction {
public:
  RigidBodyObjective(Matter &matter, const Parameters &parameters,
                     const std::vector<BasinHoppingJob::RigidBody> &bodies)
      : ObjectiveFunction(parameters), matter_(matter), bodies_(bodies) {
    x_ = VectorXd::Zero(6 * static_cast<Eigen::Index>(bodies_.size()));
    reference_.resize(bodies_.size());
    const AtomMatrix pos = matter_.getPositions();
    for (size_t b = 0; b < bodies_.size(); ++b) {
      Vector3d com = Vector3d::Zero();
      for (int site = 0; site < 3; ++site) {
        const long atom = bodies_[b].atom[site];
        com[0] += pos(atom, 0);
        com[1] += pos(atom, 1);
        com[2] += pos(atom, 2);
      }
      com /= 3.0;
      x_.segment<3>(static_cast<Eigen::Index>(6 * b)) = com;
      for (int site = 0; site < 3; ++site) {
        const long atom = bodies_[b].atom[site];
        reference_[b][site] =
            Vector3d(pos(atom, 0), pos(atom, 1), pos(atom, 2)) - com;
      }
    }
    apply(x_);
  }

  double getEnergy() override { return matter_.getPotentialEnergy(); }

  VectorXd getGradient(bool /*fdstep*/ = false) override {
    const AtomMatrix forces = matter_.getForces();
    VectorXd gradient =
        VectorXd::Zero(6 * static_cast<Eigen::Index>(bodies_.size()));
    for (size_t b = 0; b < bodies_.size(); ++b) {
      const Vector3d omega =
          x_.segment<3>(static_cast<Eigen::Index>(6 * b + 3));
      Vector3d dCom = Vector3d::Zero();
      Vector3d dOmega = Vector3d::Zero();
      for (int site = 0; site < 3; ++site) {
        const long atom = bodies_[b].atom[site];
        const Vector3d dEdr(-forces(atom, 0), -forces(atom, 1),
                            -forces(atom, 2));
        dCom += dEdr;
        const Matrix3d jacobian =
            BasinHoppingJob::rotateAngleAxisJacobian(reference_[b][site], omega);
        dOmega += jacobian.transpose() * dEdr;
      }
      gradient.segment<3>(static_cast<Eigen::Index>(6 * b)) = dCom;
      gradient.segment<3>(static_cast<Eigen::Index>(6 * b + 3)) = dOmega;
    }
    return gradient;
  }

  void setPositions(const VectorXd &x) override {
    x_ = x;
    apply(x_);
  }
  VectorXd getPositions() override { return x_; }
  int degreesOfFreedom() override { return static_cast<int>(x_.size()); }
  bool isConverged() override {
    return getConvergence() < params.optimizer_options().converged_force;
  }
  double getConvergence() override {
    const VectorXd gradient = getGradient();
    const double n = static_cast<double>(gradient.size());
    return n > 0.0 ? gradient.norm() / std::sqrt(n) : 0.0;
  }
  VectorXd difference(const VectorXd &a, const VectorXd &b) override {
    return a - b;
  }
  bool getPeriodic() const override { return false; }

private:
  void apply(const VectorXd &x) {
    AtomMatrix pos = matter_.getPositions();
    for (size_t b = 0; b < bodies_.size(); ++b) {
      const Vector3d com = x.segment<3>(static_cast<Eigen::Index>(6 * b));
      const Vector3d omega =
          x.segment<3>(static_cast<Eigen::Index>(6 * b + 3));
      for (int site = 0; site < 3; ++site) {
        const Vector3d placed =
            com + BasinHoppingJob::rotateAngleAxis(reference_[b][site], omega);
        const long atom = bodies_[b].atom[site];
        pos(atom, 0) = placed[0];
        pos(atom, 1) = placed[1];
        pos(atom, 2) = placed[2];
      }
    }
    matter_.setPositions(pos);
  }

  Matter &matter_;
  std::vector<BasinHoppingJob::RigidBody> bodies_;
  std::vector<std::array<Vector3d, 3>> reference_;
  VectorXd x_;
};

} // namespace

bool BasinHoppingJob::configureStockQuench() {
  auto &opt = ParametersLoadAccess::optimizer_options(params);
  if (opt.method != OptType::CG || opt.max_iterations != 1000 ||
      opt.converged_force != 0.01 || opt.max_move != 0.2) {
    return false;
  }
  opt.method = OptType::LBFGS;
  opt.lbfgs.accept = "energy";
  if (opt.lbfgs.memory < 20) {
    opt.lbfgs.memory = 20;
  }
  opt.lbfgs.auto_scale = true;
  opt.lbfgs.curvature = "damped";
  // RMS so the same tolerance means the same force per coordinate at
  // every cluster size. The norm of the force vector tightens as N grows.
  opt.convergence_metric = "rms";
  // A larger trial step reaches the basin in fewer iterations. Energy
  // acceptance still rejects a step that raises the energy.
  opt.max_move = 0.4;
  return true;
}

bool BasinHoppingJob::prepareRigidWater() {
  rigidBodies_.clear();
  rigidWater_ = false;
  const long nAtoms = current->numberOfAtoms();
  if (nAtoms < 3 || nAtoms % 3 != 0) {
    return false;
  }
  std::vector<long> oxygens;
  std::vector<long> hydrogens;
  oxygens.reserve(static_cast<size_t>(nAtoms / 3));
  hydrogens.reserve(static_cast<size_t>(2 * nAtoms / 3));
  for (long atom = 0; atom < nAtoms; ++atom) {
    if (current->getFixed(atom)) {
      return false;
    }
    const long z = current->getAtomicNr(atom);
    if (z == 8) {
      oxygens.push_back(atom);
    } else if (z == 1) {
      hydrogens.push_back(atom);
    } else {
      return false;
    }
  }
  if (oxygens.empty() || hydrogens.size() != oxygens.size() * 2) {
    return false;
  }
  const AtomMatrix pos = current->getPositions();
  // A contact can sit closer than a real O-H bond. Pair by water geometry:
  // both distances near 0.9572 Å and the angle near 104.52 degrees.
  constexpr double kRoh = 0.9572;
  constexpr double kTheta = 104.52 * M_PI / 180.0;
  constexpr double kCut = 1.6;
  constexpr double kAngleWeight = 0.05;
  struct Candidate {
    double cost;
    long oxygen;
    long first;
    long second;
  };
  std::vector<Candidate> candidates;
  candidates.reserve(oxygens.size() * hydrogens.size() * hydrogens.size() / 2);
  for (long oxygen : oxygens) {
    const Vector3d origin(pos(oxygen, 0), pos(oxygen, 1), pos(oxygen, 2));
    for (size_t left = 0; left < hydrogens.size(); ++left) {
      const long h1 = hydrogens[left];
      const Vector3d bond1(pos(h1, 0) - origin[0], pos(h1, 1) - origin[1],
                           pos(h1, 2) - origin[2]);
      const double d1 = bond1.norm();
      if (!(d1 > 1.0e-8) || d1 >= kCut) {
        continue;
      }
      for (size_t right = left + 1; right < hydrogens.size(); ++right) {
        const long h2 = hydrogens[right];
        const Vector3d bond2(pos(h2, 0) - origin[0], pos(h2, 1) - origin[1],
                             pos(h2, 2) - origin[2]);
        const double d2 = bond2.norm();
        if (!(d2 > 1.0e-8) || d2 >= kCut) {
          continue;
        }
        const double cosine = std::clamp(bond1.dot(bond2) / (d1 * d2), -1.0, 1.0);
        const double angle = std::acos(cosine);
        const double cost = (d1 - kRoh) * (d1 - kRoh) + (d2 - kRoh) * (d2 - kRoh) +
                            kAngleWeight * (angle - kTheta) * (angle - kTheta);
        const long first = std::min(h1, h2);
        const long second = std::max(h1, h2);
        candidates.push_back(Candidate{cost, oxygen, first, second});
      }
    }
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate &left, const Candidate &right) {
              if (left.cost != right.cost) {
                return left.cost < right.cost;
              }
              if (left.oxygen != right.oxygen) {
                return left.oxygen < right.oxygen;
              }
              if (left.first != right.first) {
                return left.first < right.first;
              }
              return left.second < right.second;
            });
  std::vector<char> oxygenUsed(static_cast<size_t>(nAtoms), 0);
  std::vector<char> hydrogenUsed(static_cast<size_t>(nAtoms), 0);
  std::vector<std::array<long, 2>> chosen(static_cast<size_t>(nAtoms));
  size_t paired = 0;
  for (const Candidate &candidate : candidates) {
    if (oxygenUsed[static_cast<size_t>(candidate.oxygen)] ||
        hydrogenUsed[static_cast<size_t>(candidate.first)] ||
        hydrogenUsed[static_cast<size_t>(candidate.second)]) {
      continue;
    }
    oxygenUsed[static_cast<size_t>(candidate.oxygen)] = 1;
    hydrogenUsed[static_cast<size_t>(candidate.first)] = 1;
    hydrogenUsed[static_cast<size_t>(candidate.second)] = 1;
    chosen[static_cast<size_t>(candidate.oxygen)] = {candidate.first, candidate.second};
    paired++;
    if (paired == oxygens.size()) {
      break;
    }
  }
  if (paired != oxygens.size()) {
    rigidBodies_.clear();
    return false;
  }
  for (long oxygen : oxygens) {
    RigidBody body;
    body.atom[0] = oxygen;
    body.atom[1] = chosen[static_cast<size_t>(oxygen)][0];
    body.atom[2] = chosen[static_cast<size_t>(oxygen)][1];
    rigidBodies_.push_back(body);
  }
  rigidWater_ = rigidBodies_.size() == oxygens.size();
  return rigidWater_;
}

void BasinHoppingJob::separateWaters(AtomMatrix &positions) const {
  // Keep atoms of different molecules from sitting on each other. A
  // hydrogen bond is longer than this, so the push only breaks clashes.
  constexpr double kMinPair = 1.6;
  constexpr int kPasses = 12;
  for (int pass = 0; pass < kPasses; ++pass) {
    bool moved = false;
    for (size_t left = 0; left < rigidBodies_.size(); ++left) {
      for (size_t right = left + 1; right < rigidBodies_.size(); ++right) {
        double best = kMinPair;
        Vector3d delta = Vector3d::Zero();
        bool hit = false;
        for (int siteA = 0; siteA < 3; ++siteA) {
          const long atomA = rigidBodies_[left].atom[siteA];
          const Vector3d ra(positions(atomA, 0), positions(atomA, 1),
                            positions(atomA, 2));
          for (int siteB = 0; siteB < 3; ++siteB) {
            const long atomB = rigidBodies_[right].atom[siteB];
            const Vector3d rb(positions(atomB, 0), positions(atomB, 1),
                              positions(atomB, 2));
            const Vector3d step = ra - rb;
            const double dist = step.norm();
            if (dist < best) {
              best = dist;
              delta = step;
              hit = true;
            }
          }
        }
        if (!hit) {
          continue;
        }
        const double dist = delta.norm();
        const Vector3d dir =
            dist > 1.0e-8 ? (delta / dist) : Vector3d(1.0, 0.0, 0.0);
        const double push = 0.5 * (kMinPair - dist) + 1.0e-6;
        for (int site = 0; site < 3; ++site) {
          const long a = rigidBodies_[left].atom[site];
          const long b = rigidBodies_[right].atom[site];
          for (int axis = 0; axis < 3; ++axis) {
            positions(a, axis) += push * dir[axis];
            positions(b, axis) -= push * dir[axis];
          }
        }
        moved = true;
      }
    }
    if (!moved) {
      break;
    }
  }
}

void BasinHoppingJob::displaceRigid(Matter *matter, double translation,
                                    double angle) {
  disp_count++;
  AtomMatrix positions = matter->getPositions();
  int only = -1;
  if (params.basin_hopping_options().single_atom_displace &&
      !rigidBodies_.empty()) {
    only = static_cast<int>(eonc::rng::randomInt(
        0, static_cast<int>(rigidBodies_.size()) - 1));
  }
  for (size_t body = 0; body < rigidBodies_.size(); ++body) {
    if (only >= 0 && static_cast<int>(body) != only) {
      continue;
    }
    Vector3d com = Vector3d::Zero();
    std::array<Vector3d, 3> local;
    for (int site = 0; site < 3; ++site) {
      const long atom = rigidBodies_[body].atom[site];
      local[static_cast<size_t>(site)] = Vector3d(
          positions(atom, 0), positions(atom, 1), positions(atom, 2));
      com += local[static_cast<size_t>(site)];
    }
    com /= 3.0;
    Vector3d shift;
    Vector3d omega;
    for (int axis = 0; axis < 3; ++axis) {
      shift[axis] = eonc::rng::randomDouble(2.0 * translation) - translation;
      omega[axis] = eonc::rng::randomDouble(2.0 * angle) - angle;
    }
    for (int site = 0; site < 3; ++site) {
      local[static_cast<size_t>(site)] -= com;
    }
    com += shift;
    for (int site = 0; site < 3; ++site) {
      const Vector3d placed =
          com + rotateAngleAxis(local[static_cast<size_t>(site)], omega);
      const long atom = rigidBodies_[body].atom[site];
      positions(atom, 0) = placed[0];
      positions(atom, 1) = placed[1];
      positions(atom, 2) = placed[2];
    }
  }
  separateWaters(positions);
  matter->setPositions(positions);
}

bool BasinHoppingJob::quench(Matter &matter, bool wantTight,
                             double bestEnergy) {
  if (!stock_) {
    matter.relax(true);
    return std::isfinite(matter.getPotentialEnergy());
  }
  if (rigidWater_) {
    return quenchRigid(matter, wantTight, bestEnergy);
  }
  return quenchCartesian(matter, wantTight, bestEnergy);
}

namespace {

bool quenchIsUsable(const Matter &matter, bool converged) {
  if (!std::isfinite(matter.getPotentialEnergy())) {
    return false;
  }
  if (converged) {
    return true;
  }
  const double n = std::max(
      1.0, 3.0 * static_cast<double>(matter.numberOfAtoms()));
  const double rms = matter.getForces().norm() / std::sqrt(n);
  return std::isfinite(rms) && rms < 0.05;
}

bool energyWorthTightening(double energy, double bestEnergy) {
  if (!std::isfinite(energy) || !std::isfinite(bestEnergy)) {
    return false;
  }
  // A loose quench that has reached the basin floor is already within
  // a few hundredths of its minimum. Only those candidates pay for a
  // second, tighter minimization.
  const double window = std::max(0.02, 1.0e-4 * std::abs(bestEnergy));
  return energy < bestEnergy + window;
}

} // namespace

bool BasinHoppingJob::quenchCartesian(Matter &matter, bool wantTight,
                                      double bestEnergy) {
  auto &opt = ParametersLoadAccess::optimizer_options(params);
  const double savedTol = opt.converged_force;
  const size_t savedIter = opt.max_iterations;
  auto stage = [&](double tolerance, size_t iterations) {
    opt.converged_force = tolerance;
    opt.max_iterations = iterations;
    return matter.relax(true);
  };
  // Loose RMS first, enough iterations to finish a cluster quench.
  // A tight RMS follows only for a hop that might improve the best
  // energy, so the other hops stay cheap.
  bool converged = stage(0.01, 400);
  double energy = matter.getPotentialEnergy();
  if (wantTight || energyWorthTightening(energy, bestEnergy)) {
    converged = stage(1.0e-3, 500);
    energy = matter.getPotentialEnergy();
    if (!converged) {
      converged = stage(1.0e-2, 200);
      energy = matter.getPotentialEnergy();
    }
  }
  opt.converged_force = savedTol;
  opt.max_iterations = savedIter;
  return quenchIsUsable(matter, converged);
}

bool BasinHoppingJob::quenchRigid(Matter &matter, bool wantTight,
                                  double bestEnergy) {
  auto &opt = ParametersLoadAccess::optimizer_options(params);
  const double savedTol = opt.converged_force;
  const size_t savedIter = opt.max_iterations;
  const bool savedProject = opt.lbfgs.project_rigid;
  // The coordinate chart is already translation and rotation. Projecting
  // Cartesian rigid modes out of those six numbers is the wrong null space.
  opt.lbfgs.project_rigid = false;
  auto stage = [&](double tolerance, size_t iterations) {
    opt.converged_force = tolerance;
    opt.max_iterations = iterations;
    auto objective =
        std::make_shared<RigidBodyObjective>(matter, params, rigidBodies_);
    LBFGS solver(objective, params);
    return solver.run(iterations, opt.max_move) > 0;
  };
  bool converged = stage(0.01, 400);
  double energy = matter.getPotentialEnergy();
  if (wantTight || energyWorthTightening(energy, bestEnergy)) {
    converged = stage(1.0e-3, 500);
    energy = matter.getPotentialEnergy();
    if (!converged) {
      converged = stage(1.0e-2, 200);
      energy = matter.getPotentialEnergy();
    }
  }
  opt.converged_force = savedTol;
  opt.max_iterations = savedIter;
  opt.lbfgs.project_rigid = savedProject;
  return quenchIsUsable(matter, converged);
}

void BasinHoppingJob::adaptStockStep(double &step, int recentAccept,
                                     int period) {
  if (period <= 0) {
    return;
  }
  const double ratio =
      static_cast<double>(recentAccept) / static_cast<double>(period);
  const double fraction = params.basin_hopping_options().adjust_fraction;
  const double target = params.basin_hopping_options().target_ratio;
  constexpr double kDispMax = 1.5;
  constexpr double kAngMax = 1.5;
  // Hold the step near its initial length. Growing it to the hard cap
  // while the temperature is high randomises the cluster every hop.
  const double stepAim = std::min(kDispMax, std::max(initialDisplacement_, 1.0e-6));
  const double angAim = std::min(kAngMax, std::max(initialAngular_, 1.0e-6));
  const double stepFloor = std::max(0.15, 0.5 * stepAim);
  const double angFloor = std::max(0.15, 0.5 * angAim);
  temperatureHeld_ = ratio + 1.0e-12 >= target && target > 0.0;
  if (ratio == 0.0) {
    temperature_ = std::min(temperature_ * 2.0, 1.0e8);
    step = stepAim;
    angularDisplacement_ = angAim;
  } else if (ratio + 1.0e-12 >= target) {
    if (step < stepAim * 0.98) {
      step = std::min(stepAim, step * (1.0 + fraction));
      angularDisplacement_ = std::min(angAim, angularDisplacement_ * (1.0 + fraction));
    } else {
      temperature_ = std::max(temperature_ / 1.25, 1.0);
      if (step > stepAim) {
        step = std::max(stepAim, step * (1.0 - fraction));
      }
      if (angularDisplacement_ > angAim) {
        angularDisplacement_ =
            std::max(angAim, angularDisplacement_ * (1.0 - fraction));
      }
    }
  } else if (step > stepFloor * 1.001) {
    step = std::max(stepFloor, step * (1.0 - fraction));
    angularDisplacement_ =
        std::max(angFloor, angularDisplacement_ * (1.0 - fraction));
  } else {
    temperature_ = std::min(temperature_ * 1.5, 1.0e8);
    step = stepAim;
    angularDisplacement_ = angAim;
  }
}

void BasinHoppingJob::liftTemperature(double deltaE) {
  const double kB = params.constants().kB;
  if (temperatureHeld_) {
    return;
  }
  if (!(deltaE > 0.0) || !(kB > 0.0) || !(temperature_ > 0.0)) {
    return;
  }
  const double kT = kB * temperature_;
  // kB is in eV/K. A surface in another unit rejects ordinary uphill
  // hops at 300 K. Double kT only while the gap is still enormous
  // compared with kT. A clash at the working temperature must not
  // ratchet kT upward; the acceptance period owns that adjustment.
  if (deltaE > 20.0 * kT) {
    temperature_ = std::min((kT * 2.0) / kB, 1.0e8);
  }
}

double BasinHoppingJob::hopTemperature() const {
  return stock_ ? temperature_ : params.main_options().temperature;
}

double BasinHoppingJob::swapProbability(const std::vector<long> &elements) const {
  const double configured = params.basin_hopping_options().swap_probability;
  if (stock_ && configured == 0.0 && elements.size() > 1) {
    return 0.5;
  }
  return configured;
}

long BasinHoppingJob::jumpStepCount() const {
  const long configured = params.basin_hopping_options().jump_steps;
  if (stock_ && configured == 0) {
    return 5;
  }
  return configured;
}

Vector3d BasinHoppingJob::rotateAngleAxis(const Vector3d &point,
                                          const Vector3d &omega) {
  const double alpha2 = omega.squaredNorm();
  const double alpha = std::sqrt(alpha2);
  double cosine = 1.0;
  double c3 = -0.5;
  double sineOver = 1.0;
  if (alpha >= 1.0e-8) {
    cosine = std::cos(alpha);
    c3 = (cosine - 1.0) / alpha2;
    sineOver = std::sin(alpha) / alpha;
  }
  const double projection = omega.dot(point);
  const Vector3d crossed = omega.cross(point);
  return cosine * point - c3 * projection * omega - sineOver * crossed;
}

Matrix3d BasinHoppingJob::rotateAngleAxisJacobian(const Vector3d &point,
                                                 const Vector3d &omega) {
  const double alpha2 = omega.squaredNorm();
  const double alpha = std::sqrt(alpha2);
  const bool small = alpha < 1.0e-8;
  double c3 = -0.5;
  double sineOver = 1.0;
  double dCosine = 0.0;
  double dC3 = 0.0;
  double dSine = 0.0;
  Vector3d dAlpha = Vector3d::Zero();
  if (!small) {
    const double cosine = std::cos(alpha);
    const double sine = std::sin(alpha);
    c3 = (cosine - 1.0) / alpha2;
    sineOver = sine / alpha;
    dCosine = -sine;
    dC3 = (-sine * alpha2 - (cosine - 1.0) * 2.0 * alpha) / (alpha2 * alpha2);
    dSine = (cosine * alpha - sine) / alpha2;
    dAlpha = omega / alpha;
  }
  const double projection = omega.dot(point);
  const Vector3d crossed = omega.cross(point);
  Matrix3d jacobian;
  for (int axis = 0; axis < 3; ++axis) {
    Vector3d unit = Vector3d::Zero();
    unit[axis] = 1.0;
    const double scaleCosine = small ? 0.0 : dCosine * dAlpha[axis];
    const double scaleC3 = small ? 0.0 : dC3 * dAlpha[axis];
    const double scaleSine = small ? 0.0 : dSine * dAlpha[axis];
    Vector3d term = scaleCosine * point;
    term -= scaleC3 * projection * omega;
    term -= c3 * point[axis] * omega;
    term -= c3 * projection * unit;
    term -= scaleSine * crossed;
    term -= sineOver * unit.cross(point);
    jacobian.col(axis) = term;
  }
  return jacobian;
}

double BasinHoppingJob::metropolisProbability(double de, double kB,
                                              double temperature) {
  if (!(de > 0.0)) {
    return 1.0;
  }
  if (!(temperature > 0.0) || !(kB > 0.0)) {
    return 0.0;
  }
  return std::exp(-de / (kB * temperature));
}

AtomMatrix BasinHoppingJob::displaceRandom(double curDisplacement) {
  disp_count++;
  // Create a random displacement
  AtomMatrix displacement;
  displacement.resize(trial->numberOfAtoms(), 3);
  displacement.setZero();
  VectorXd distvec = calculateDistanceFromCenter(current.get());
  // Coincident atoms share a zero radius. Dividing by it writes NaN and
  // the hop cannot leave the center. Those atoms are the outer shell, so
  // the step stays the unscaled displacement.
  const double radius = distvec.size() > 0 ? distvec.maxCoeff() : 0.0;
  int num = trial->numberOfAtoms();
  int m = 0;
  if (params.basin_hopping_options().single_atom_displace) {
    m = eonc::rng::randomInt(0, trial->numberOfAtoms() - 1);
    num = m + 1;
  }

  for (int i = m; i < num; i++) {
    double dist = distvec(i);
    double disp = 0.0; // displacement size, possibly scaled

    if (!trial->getFixed(i)) {
      const std::string &algorithm =
          params.basin_hopping_options().displacement_algorithm;
      if (algorithm == "standard") {
        disp = curDisplacement;
      }
      // scale displacement linearly with the particle radius
      else if (algorithm == "linear") {
        disp =
            radius == 0.0 ? curDisplacement : curDisplacement * (dist / radius);
      }
      // scale displacement quadratically with the particle radius
      else if (algorithm == "quadratic") {
        if (radius == 0.0) {
          disp = curDisplacement;
        } else {
          const double scale = dist / radius;
          disp = curDisplacement * scale * scale;
        }
      } else {
        if (quill::Logger *trace = eonc::log::traceback()) {
          log = trace;
          QUILL_LOG_CRITICAL(log, "Unknown displacement_algorithm\n");
        }
        throw std::invalid_argument(
            std::format("[Basin Hopping] unknown displacement_algorithm: {}",
                        params.basin_hopping_options().displacement_algorithm));
      }
      for (int j = 0; j < 3; j++) {
        if (params.basin_hopping_options().displacement_distribution ==
            "uniform") {
          displacement(i, j) = eonc::rng::randomDouble(2 * disp) - disp;
        } else if (params.basin_hopping_options().displacement_distribution ==
                   "gaussian") {
          displacement(i, j) = eonc::rng::gaussRandom(0.0, disp);
        } else {
          if (quill::Logger *trace = eonc::log::traceback()) {
            log = trace;
            QUILL_LOG_CRITICAL(log, "Unknown displacement_distribution\n");
          }
          throw std::invalid_argument(std::format(
              "[Basin Hopping] unknown displacement_distribution: {}",
              params.basin_hopping_options().displacement_distribution));
        }
      }
    }
  }
  return displacement;
}

void BasinHoppingJob::randomSwap(Matter *matter) {
  swap_count++;
  std::vector<long> Elements;
  Elements = getElements(matter);

  long ela;
  long elb;
  long ia = eonc::rng::randomInt(0, Elements.size() - 1);
  ela = Elements.at(ia);
  Elements.erase(Elements.begin() + ia);

  long ib = eonc::rng::randomInt(0, Elements.size() - 1);
  elb = Elements.at(ib);

  int changera = 0;
  int changerb = 0;

  changera = eonc::rng::randomInt(0, matter->numberOfAtoms() - 1);
  int guard = 0;
  while ((matter->getAtomicNr(changera) != ela || matter->getFixed(changera)) &&
         guard < 10000) {
    changera = eonc::rng::randomInt(0, matter->numberOfAtoms() - 1);
    guard++;
  }
  changerb = eonc::rng::randomInt(0, matter->numberOfAtoms() - 1);
  guard = 0;
  while ((matter->getAtomicNr(changerb) != elb || matter->getFixed(changerb) ||
          changerb == changera) &&
         guard < 10000) {
    changerb = eonc::rng::randomInt(0, matter->numberOfAtoms() - 1);
    guard++;
  }
  if (matter->getAtomicNr(changera) != ela || matter->getFixed(changera) ||
      matter->getAtomicNr(changerb) != elb || matter->getFixed(changerb)) {
    return;
  }

  double posax = matter->getPosition(changera, 0);
  double posay = matter->getPosition(changera, 1);
  double posaz = matter->getPosition(changera, 2);

  matter->setPosition(changera, 0, matter->getPosition(changerb, 0));
  matter->setPosition(changera, 1, matter->getPosition(changerb, 1));
  matter->setPosition(changera, 2, matter->getPosition(changerb, 2));

  matter->setPosition(changerb, 0, posax);
  matter->setPosition(changerb, 1, posay);
  matter->setPosition(changerb, 2, posaz);
}

std::vector<long> BasinHoppingJob::getElements(Matter *matter) {
  // Z is 0..118. A 118-slot table stores 118 one past the end.
  constexpr int kElementSlots = 119;
  std::array<int, kElementSlots> allElements{};
  std::vector<long> elements;

  for (long y = 0; y < matter->numberOfAtoms(); ++y) {
    if (!matter->getFixed(y)) {
      const long z = matter->getAtomicNr(y);
      if (z >= 0 && z < kElementSlots) {
        allElements[static_cast<size_t>(z)] = 1;
      }
    }
  }

  for (int i = 0; i < kElementSlots; ++i) {
    if (allElements[static_cast<size_t>(i)] != 0) {
      elements.push_back(i);
    }
  }

  return elements;
}

VectorXd BasinHoppingJob::calculateDistanceFromCenter(Matter *matter) {
  AtomMatrix pos = matter->getPositions();
  Vector3d cen(0, 0, 0);
  int num = matter->numberOfAtoms();

  cen = pos.colwise().sum() / static_cast<double>(num);

  VectorXd dist(num);

  for (int n = 0; n < num; n++) {
    pos.row(n) -= cen;
    dist(n) = pos.row(n).norm();
  }

  return dist;
}

} // namespace eonc

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
#include "eon/MinModeSaddleSearch.h"
#include "eon/ConjugateGradients.h"
#include "eon/EigenmodeStrategy.h"
#include "eon/EonLogger.h"
#include "eon/EpiCenters.h"
#include "eon/GeometryAnalysis.h"
#include "eon/HelperFunctions.h"
#include "eon/ObjectiveFunction.h"
#include "eon/PrfoRefine.h"
#include "eon/SaddleSearchJob.h"
#include "eon/SafeMath.h"
#include "eon/eonExceptions.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

namespace eonc {

namespace {

// Twenty halvings leave 2^-20 of the attempted step, about 10^-6.
// For the default max_move of 0.2 Å that remainder is 2e-7 Å, which
// is still above the motion floor, so this ceiling is what stops a
// long step. The floor stops a shorter one first.
constexpr int kMaxEnergyHalvings = 20;

// Each stall reseeds and halves the cap, and the seed sign alternates.
// Six stalls is three tries of each sign, and the cap has dropped by
// 64. Both seeds have been repeated at three step sizes.
constexpr int kMaxEnergyStalls = 6;

// A kept fraction below 1/64 has been halved more than six times.
// The mode was estimated for the full step. A piece that small is a
// different displacement, so the rotational CG or L-BFGS state is
// dropped and the next estimate starts from the original tangent.
constexpr double kMaxEnergyReseedAlpha = 1.0 / 64.0;

// Smallest translation cap, in Å. The default finite-difference
// distance is 1e-2 Å and the saddle ini files use 1e-3 Å. A cap
// below 1e-4 Å is smaller than those probes. Holding it here keeps
// the next attempt a real move, and it stays above the motion floor.
constexpr double kMaxEnergyStepFloor = 1.0e-4;

// Max-atom motion, in Å, below which a probe repeats a geometry
// already evaluated. L-BFGS uses the same 1e-8 scale to drop a
// secant pair whose overlap is numerical noise.
constexpr double kMaxEnergyMotionFloor = 1.0e-8;

// One doubling per accepted under-limit step, the inverse of the
// probe halving. Each legal step gives back one factor the wall
// removed. The default max_move of 0.2 Å returns from the step
// floor in ceil(log2(0.2/1e-4)) = 11 legal steps.
constexpr double kMaxEnergyCapRegrowth = 2.0;

// Crossings of max_energy in one search, whether a probe was found
// or the centre was restored. Six consecutive stalls already stop a
// search that never finds a legal point. A found probe, or one
// under-limit step, resets that counter and the cap doubles back,
// so a search sitting on the limit crosses once a regrowth cycle
// until max_iterations. Sixteen is more than one full regrowth of
// a 0.2 Å cap from the step floor (11) and more than the stall
// budget. The LJ38 dimer start this recovery finishes crosses
// once, so the cap does not cut that search off. Past sixteen the
// search is pinned and ends with the high-energy status. The
// default max_iterations is 1000, and each crossing can spend 20
// probes.
constexpr int kMaxEnergyWallEvents = 16;

} // namespace

class MinModeObjectiveFunction : public ObjectiveFunction {
private:
  std::shared_ptr<Matter> matter;
  AtomMatrix eigenvector;
  std::shared_ptr<EigenmodeStrategy> minModeMethod;
  int iteration{0};
  double reactantEnergy{0.0};

public:
  MinModeObjectiveFunction(
      std::shared_ptr<Matter> matterPassed,
      std::shared_ptr<EigenmodeStrategy> minModeMethodPassed,
      AtomMatrix modePassed, double reactantEnergyPassed,
      const Parameters &paramsPassed)
      : ObjectiveFunction(paramsPassed),
        matter{std::move(matterPassed)},
        eigenvector{std::move(modePassed)},
        minModeMethod{minModeMethodPassed},
        reactantEnergy{reactantEnergyPassed} {}

  ~MinModeObjectiveFunction() override = default;

  // |C_tau| is the curvature of the potential along the mode. The
  // translated force is F - 2 (F·tau) tau, and that force has a
  // different curvature. Improved dimer and Lanczos do not report the
  // mode scale, so their first L-BFGS step measures F_eff. That probe
  // is what the LJ13 and Pt heptamer force-call counts follow.
  //
  // The classic dimer (improved = false) is the adatom AKMC catalogue
  // path. It reports the negative mode curvature. L-BFGS keeps the
  // capped scale for a mild mode whose uncapped step exceeds max_move,
  // and probes F_eff for a steep overrun.
  //
  // With rgsaddle the rotation of every dimer and Lanczos search is
  // one solver type, so the C++ class does not say which search this
  // is. The scale follows the dimer setting.
  std::optional<double> knownCurvature() const override {
    if (iteration == 0 || !minModeMethod) {
      return std::nullopt;
    }
    if (params.saddle_search_options().minmode_method !=
            LowestEigenmode::MINMODE_DIMER ||
        params.dimer_options().improved) {
      return std::nullopt;
    }
    const double c = eonc::eigenmodeGetEigenvalue(*minModeMethod);
    if (!std::isfinite(c) || c >= 0.0) {
      return std::nullopt;
    }
    return -c;
  }

  // The translation force has its component along the mode reversed, so
  // the search climbs the potential energy that getEnergy() returns.
  bool gradientMatchesEnergy() const override { return false; }

  // A new centre and the mode solver's first displaced image along the
  // current mode are independent evaluations; a batching potential takes
  // them together, and compute() then finds the image evaluated.
  // Finite-difference gradient probes never call compute(), so they
  // evaluate the centre alone.
  void evaluateCentre() {
    if (matter->needsForceUpdate()) {
      minModeMethod->evaluateWithProbe(matter, eigenvector);
    }
  }

  // The centre (and a batching solver's first image) is evaluated
  // here. FIRE then skips getGradient() when this is true, so the
  // mode is not rotated at a point the search will pull back.
  bool pastClimbCeiling() override {
    evaluateCentre();
    return getEnergy() - reactantEnergy >
           params.saddle_search_options().max_energy;
  }

  VectorXd getGradient(bool fdstep = false) {
    if (!fdstep || iteration == 0) {
      evaluateCentre();
    }
    AtomMatrix force = matter->getForces();

    if (!fdstep || iteration == 0) {
      eonc::eigenmodeCompute(*minModeMethod, matter, eigenvector);
      // Check if ImprovedDimer lost the mode
      auto *dimer = eonc::asImprovedDimer(*minModeMethod);
      if (dimer && !dimer->rotationDidConverge) {
        if (dimer->getEigenvalue() < 0.0) {
          eigenvector = eonc::eigenmodeGetEigenvector(*minModeMethod);
          EONC_LOG_DEBUG(
              "[MinMode] Dimer restored to best state with C_tau={:.4f}",
              dimer->getEigenvalue());
          throw eonc::DimerModeRestoredException();
        } else {
          throw eonc::DimerModeLostException();
        }
      }
      iteration++;
    }

    eigenvector = eonc::eigenmodeGetEigenvector(*minModeMethod);
    double eigenvalue = eonc::eigenmodeGetEigenvalue(*minModeMethod);
    // OCI-NEB returns a positive-curvature image to the band before
    // translating it away from the current path.
    if (eigenvalue > 0.0 &&
        params.neb_options().climbing_image.ocineb.use_mmf) {
      throw eonc::DimerModeLostException();
    }

    AtomMatrix proj = matDot(force, eigenvector) *
                      eonc::safemath::safe_normalized(eigenvector);

    if (eigenvalue > 0.0) {
      if (params.saddle_search_options().perp_force_ratio > 0.0) {
        double d = params.saddle_search_options().perp_force_ratio;
        force = d * force - (1.0 + d) * proj;
      } else if (params.saddle_search_options().confine_positive.enabled) {
        if (params.saddle_search_options().confine_positive.bowl_breakout) {
          AtomMatrix forceTemp = matter->getForces();
          const long nAtoms = matter->numberOfAtoms();
          int nBowlActive = static_cast<int>(std::min<long>(
              params.saddle_search_options().confine_positive.bowl_active,
              nAtoms));
          if (nBowlActive <= 0) {
            force.setZero();
          } else {
            std::vector<int> indices_max(nBowlActive);

            // Find the nBowlActive atoms with largest forces
            for (int j = 0; j < nBowlActive; j++) {
              double f_max = forceTemp.row(0).norm();
              int i_max = 0;
              for (long i = 0; i < matter->numberOfAtoms(); i++) {
                if (f_max < forceTemp.row(i).norm()) {
                  f_max = forceTemp.row(i).norm();
                  i_max = static_cast<int>(i);
                }
              }
              forceTemp.row(i_max).setZero();
              indices_max[j] = i_max;
            }
            forceTemp.setZero();
            for (int j = 0; j < nBowlActive; j++) {
              forceTemp.row(indices_max[j]) = -proj.row(indices_max[j]);
            }
            force = forceTemp;
          }
        } else {
          int sufficientForce = 0;
          double minForce =
              params.saddle_search_options().confine_positive.min_force;
          const long maxBoostTries = std::max(
              3 * matter->numberOfAtoms(),
              params.saddle_search_options().confine_positive.min_active);
          long boostTries = 0;
          while (
              sufficientForce <
                  params.saddle_search_options().confine_positive.min_active &&
              boostTries < maxBoostTries) {
            sufficientForce = 0;
            force = matter->getForces();
            for (long i = 0; i < matter->numberOfAtoms(); i++) {
              for (int k = 0; k < 3; k++) {
                if (std::abs(force(i, k)) < minForce) {
                  force(i, k) = 0;
                } else {
                  sufficientForce++;
                  force(i, k) =
                      -params.saddle_search_options().confine_positive.boost *
                      proj(i, k);
                }
              }
            }
            minForce *=
                params.saddle_search_options().confine_positive.scale_ratio;
            boostTries++;
          }
        }
      } else {
        force = -proj;
      }
    } else {
      force += -2.0 * proj;
    }

    VectorXd forceV = VectorXd::Map(force.data(), 3 * matter->numberOfAtoms());
    return -forceV;
  }

  // No prefetch here: an energy-accepting optimizer reads trial points
  // that may be rejected.
  double getEnergy() { return matter->getPotentialEnergy(); }
  void setPositions(const VectorXd &x) { matter->setPositionsV(x); }
  VectorXd getPositions() { return matter->getPositionsV(); }
  int degreesOfFreedom() { return 3 * matter->numberOfAtoms(); }
  bool isConverged() {
    return getConvergence() < params.saddle_search_options().converged_force;
  }

  double getConvergence() {
    evaluateCentre();
    if (params.optimizer_options().convergence_metric == "norm") {
      return matter->getForcesFreeV().norm();
    } else if (params.optimizer_options().convergence_metric == "max_atom") {
      return matter->maxForce();
    } else if (params.optimizer_options().convergence_metric ==
               "max_component") {
      return matter->getForces().cwiseAbs().maxCoeff();
    } else {
      EONC_LOG_CRITICAL("[MinModeSaddleSearch] unknown convergence metric: {}",
                        params.optimizer_options().convergence_metric);
      throw std::invalid_argument(
          std::format("[MinModeSaddleSearch] unknown convergence_metric: {}",
                      params.optimizer_options().convergence_metric));
    }
  }

  VectorXd difference(const VectorXd &a, const VectorXd &b) {
    return matter->pbcV(a - b);
  }

  // Drop the rotation memory and start the next estimate from seed.
  // A translation that has walked into a wall keeps the old rotational
  // CG or L-BFGS state, and that state repeats the same climb.
  void reseedMode(const AtomMatrix &seed) {
    eigenvector = seed;
    if (auto *dimer = eonc::asImprovedDimer(*minModeMethod)) {
      dimer->init_cg = true;
      dimer->init_lbfgs = true;
      dimer->s.clear();
      dimer->y.clear();
      dimer->rho.clear();
      dimer->gradients.clear();
      dimer->positions.clear();
    }
  }
};

MinModeSaddleSearch::MinModeSaddleSearch(std::shared_ptr<Matter> matterPassed,
                                         AtomMatrix modePassed,
                                         double reactantEnergyPassed,
                                         const Parameters &parametersPassed,
                                         std::shared_ptr<Potential> potPassed)
    : SaddleSearchMethod(potPassed, parametersPassed),
      matter{matterPassed} {
  eonc::helpers::requireKnownConvergenceMetric(
      params.optimizer_options().convergence_metric, "[MinModeSaddleSearch]");
  reactantEnergy = reactantEnergyPassed;
  mode = modePassed;
  initialTangent_ = modePassed;
  status = STATUS_GOOD;
  iteration = 0;

  minModeMethod = eonc::buildEigenmodeStrategy(matter, params, pot);

  // Set reference mode for ImprovedDimer (prevents mode switching)
  if (auto *dimer = eonc::asImprovedDimer(*minModeMethod)) {
    VectorXd refVec = VectorXd::Map(mode.data(), 3 * matter->numberOfAtoms());
    refVec = refVec.array() * matter->getFreeV().array();
    dimer->setReferenceMode(refVec);
  }
}

int MinModeSaddleSearch::run() {
  return run(params.saddle_search_options().max_iterations);
}

int MinModeSaddleSearch::runRetainFrames(long max_iterations_override) {
  retain_climb_frames_ = true;
  climb_frames_.clear();
  const long maxIter = max_iterations_override < 0
                           ? params.saddle_search_options().max_iterations
                           : max_iterations_override;
  const int st = run(maxIter);
  retain_climb_frames_ = false;
  return st;
}

int MinModeSaddleSearch::run(long max_iterations_override) {
  long effectiveMaxIter = max_iterations_override;
  QUILL_LOG_DEBUG(
      log, "Saddle point search started from reactant with energy {} eV.",
      reactantEnergy);

  int optStatus;
  bool firstIteration = true;
  const char *forceLabel =
      params.optimizer_options().convergence_metric_label.c_str();

  if (params.saddle_search_options().minmode_method ==
      LowestEigenmode::MINMODE_GPRDIMER) {
    QUILL_LOG_DEBUG(
        log, "================= Using the GP Dimer Library =================");
    eonc::eigenmodeCompute(*minModeMethod, matter, mode);
    if (eonc::eigenmodeGetEigenvalue(*minModeMethod) > 0) {
      QUILL_LOG_DEBUG(log, "GPR eigenvalue: {}",
                      eonc::eigenmodeGetEigenvalue(*minModeMethod));
      status = STATUS_NONNEGATIVE_ABORT;
      return status;
    }
    if (getEigenvalue() > 0.0 && status == STATUS_GOOD) {
      QUILL_LOG_DEBUG(log, "[MinModeSaddleSearch] eigenvalue not negative");
      status = STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE;
    }
    if (std::abs(eonc::eigenmodeGetEigenvalue(*minModeMethod)) <
        params.saddle_search_options().zero_mode_abort_curvature) {
      QUILL_LOG_DEBUG(log, "Zero mode eigenvalue: {}",
                      eonc::eigenmodeGetEigenvalue(*minModeMethod));
      status = STATUS_ZEROMODE_ABORT;
    }
    iteration = eonc::eigenmodeTotalIterations(*minModeMethod);
    forcecalls = eonc::eigenmodeTotalForceCalls(*minModeMethod);
  } else {

    if (params.saddle_search_options().minmode_method ==
        LowestEigenmode::MINMODE_DIMER) {
      QUILL_LOG_INFO(log,
                     "[Dimer]  {:9s}   {:9s}   {:10s}   {:18s}   {:9s}   "
                     "{:7s}   {:6s}   {:4s}   {:5s}\n",
                     "Step", "Step Size", "Delta E", forceLabel, "Curvature",
                     "Torque", "Angle", "Rots", "Align");
    } else if (params.saddle_search_options().minmode_method ==
               LowestEigenmode::MINMODE_LANCZOS) {
      QUILL_LOG_INFO(
          log,
          "[Lanczos]  {:9s} {:9s} {:10s} {:18s} {:9s} {:10s} {:7s} {:5s}\n",
          "Step", "Step Size", "Delta E", forceLabel, "Curvature", "Rel Change",
          "Angle", "Iters");
    } else if (params.saddle_search_options().minmode_method ==
               LowestEigenmode::MINMODE_GPRDIMER) {
      QUILL_LOG_INFO(log,
                     "[GPRDimer]  {:9s}   {:9s}   {:10s}   {:18s}   {:9s}  "
                     " {:7s}   {:6s}   {:4s}\n",
                     "Step", "Step Size", "Delta E", forceLabel, "Curvature",
                     "Torque", "Angle", "Rots");
    }

    std::string climbLabel = "climb";
    std::string climbDatFilename = "climb.dat";

    AtomMatrix initialPosition = matter->getPositions();

    auto objf = std::make_shared<MinModeObjectiveFunction>(
        matter, minModeMethod, mode, reactantEnergy, params);
    auto write_climb_frame = [&](uint64_t frameIndex, bool append,
                                 double stepSize, double de, double conv,
                                 double eigenval, double torque, double angle,
                                 long rotations) {
      eonc::io::ConFrameMetadata metadata;
      metadata.frame_index = frameIndex;
      metadata.energy = matter->getPotentialEnergy();
      metadata.scalars.push_back({"step_size", stepSize});
      metadata.scalars.push_back({"delta_e", de});
      metadata.scalars.push_back({"convergence", conv});
      metadata.scalars.push_back({"eigenvalue", eigenval});
      metadata.scalars.push_back({"torque", torque});
      metadata.scalars.push_back({"angle", angle});
      metadata.scalars.push_back({"rotations", static_cast<double>(rotations)});
      if (retain_climb_frames_) {
        climb_frames_.push_back(eonc::io::matterToConFrame(*matter, &metadata));
      }
      if (params.debug_options().write_movies) {
        if (!eonc::io::io_ok(
                matter->matter2con(climbLabel, append, &metadata))) {
          QUILL_LOG_WARNING(log, "Failed to write climb movie frame {}",
                            climbLabel);
        }
        eonc::helpers::saveMode(std::format("mode_{:03}.dat", frameIndex),
                                matter,
                                eonc::eigenmodeGetEigenvector(*minModeMethod));
      }

      if (params.debug_options().write_deprecated_outs) {
        std::ofstream climbDat(climbDatFilename,
                               append ? (std::ios::binary | std::ios::app)
                                      : std::ios::binary);
        if (climbDat) {
          if (!append) {
            climbDat << "iteration\tstep_size\tdelta_e\tconvergence"
                        "\teigenvalue\ttorque\tangle\trotations\n";
          }
          climbDat << std::format("{}\t{:.7e}\t{:.6f}\t{:.5e}\t{:.6f}"
                                  "\t{:.6f}\t{:.4f}\t{}\n",
                                  frameIndex, stepSize, de, conv, eigenval,
                                  torque, angle, rotations);
        }
      }
    };
    if (params.debug_options().write_movies || retain_climb_frames_) {
      write_climb_frame(0, false, 0.0, 0.0, objf->getConvergence(),
                        eonc::eigenmodeGetEigenvalue(*minModeMethod), 0.0, 0.0,
                        0);
    }
    if (params.saddle_search_options().nonnegative_displacement_abort) {
      objf->getGradient();
      if (eonc::eigenmodeGetEigenvalue(*minModeMethod) > 0) {
        QUILL_LOG_DEBUG(log, "Nonnegative eigenvalue: {}",
                        eonc::eigenmodeGetEigenvalue(*minModeMethod));
        status = STATUS_NONNEGATIVE_ABORT;
        return status;
      }
    }

    auto optim = eonc::helpers::create::mkOptim(
        objf, params.optimizer_options().method, params);

    // Stock max_move stays the cap until a translation crosses max_energy.
    // A search that never crosses it takes the same steps as before.
    const double stockMove = params.optimizer_options().max_move;
    double translationCap = stockMove;
    int wallStalls = 0;
    int wallEvents = 0;
    int reseeds = 0;
    int probeForceCalls = 0;
    auto probedEnergy = [&]() {
      const long beforeCalls = matter->getForceCalls();
      const double energy = objf->getEnergy();
      probeForceCalls +=
          static_cast<int>(matter->getForceCalls() - beforeCalls);
      return energy;
    };
    auto pullUnderEnergy = [&](const VectorXd &before) -> bool {
      const double limit = params.saddle_search_options().max_energy;
      const VectorXd after = objf->getPositions();
      const VectorXd delta = objf->difference(after, before);
      double alpha = 0.5;
      bool found = false;
      for (int probe = 0; probe < kMaxEnergyHalvings; ++probe) {
        if (eonc::geometry::maxAtomMotionV(alpha * delta) <
            kMaxEnergyMotionFloor) {
          break;
        }
        objf->setPositions(before + alpha * delta);
        if (probedEnergy() - reactantEnergy <= limit) {
          found = true;
          break;
        }
        alpha *= 0.5;
      }
      auto freshMode = [&]() {
        AtomMatrix seed = initialTangent_;
        if (reseeds % 2 == 1) {
          seed = -seed;
        }
        objf->reseedMode(seed);
        ++reseeds;
      };
      if (found) {
        wallStalls = 0;
        ++wallEvents;
        const double accepted = eonc::geometry::maxAtomMotionV(alpha * delta);
        translationCap =
            std::min(stockMove, std::max(accepted, kMaxEnergyStepFloor));
        if (alpha < kMaxEnergyReseedAlpha) {
          freshMode();
          QUILL_LOG_DEBUG(
              log,
              "[MinMode] energy above max_energy {:.4f}, shortened "
              "the translation to {:.3e} and re-initialised the mode",
              limit, accepted);
        } else {
          QUILL_LOG_DEBUG(log,
                          "[MinMode] energy above max_energy {:.4f}, shortened "
                          "the translation to {:.3e}",
                          limit, accepted);
        }
        if (wallEvents >= kMaxEnergyWallEvents) {
          QUILL_LOG_DEBUG(log,
                          "[MinMode] {} crossings of max_energy {:.4f}, "
                          "ending the search",
                          wallEvents, limit);
          return false;
        }
        return true;
      }
      objf->setPositions(before);
      if (probedEnergy() - reactantEnergy > limit) {
        return false;
      }
      ++wallEvents;
      ++wallStalls;
      translationCap = std::min(
          stockMove, std::max(translationCap * 0.5, kMaxEnergyStepFloor));
      freshMode();
      QUILL_LOG_DEBUG(log,
                      "[MinMode] translation stays above max_energy {:.4f}, "
                      "restored the centre and re-initialised the mode",
                      limit);
      if (wallEvents >= kMaxEnergyWallEvents) {
        QUILL_LOG_DEBUG(log,
                        "[MinMode] {} crossings of max_energy {:.4f}, "
                        "ending the search",
                        wallEvents, limit);
        return false;
      }
      return wallStalls < kMaxEnergyStalls;
    };

    while (!objf->isConverged() || iteration == 0) {
      matter->pollCancel("saddle");

      if (!firstIteration) {

        if (params.saddle_search_options().nonlocal_count_abort != 0) {
          long nm = eonc::geometry::numAtomsMoved(
              initialPosition - matter->getPositions(),
              params.saddle_search_options().nonlocal_distance_abort);
          if (nm >= params.saddle_search_options().nonlocal_count_abort) {
            status = STATUS_NONLOCAL_ABORT;
            break;
          }
        }

        if (std::abs(eonc::eigenmodeGetEigenvalue(*minModeMethod)) <
            params.saddle_search_options().zero_mode_abort_curvature) {
          QUILL_LOG_DEBUG(log, "Zero mode eigenvalue: {}",
                          eonc::eigenmodeGetEigenvalue(*minModeMethod));
          status = STATUS_ZEROMODE_ABORT;
          break;
        }
      }
      firstIteration = false;

      if (iteration >= effectiveMaxIter) {
        status = STATUS_BAD_MAX_ITERATIONS;
        break;
      }

      const VectorXd posBefore = objf->getPositions();
      AtomMatrix pos = matter->getPositions();

      try {
        if (params.saddle_search_options().confine_positive.bowl_breakout &&
            eonc::eigenmodeGetEigenvalue(*minModeMethod) > 0 &&
            params.optimizer_options().method == OptType::CG) {
          optStatus = optim->step(-translationCap);
        } else {
          optStatus = optim->step(translationCap);
        }
      } catch (const eonc::DimerModeRestoredException &) {
        QUILL_LOG_DEBUG(
            log, "Dimer restored to best state. Checking convergence...");
        status = objf->isConverged() ? STATUS_GOOD : STATUS_DIMER_RESTORED_BEST;
        break;
      } catch (const eonc::DimerModeLostException &) {
        QUILL_LOG_WARNING(log, "Dimer lost mode completely. Aborting.");
        status = STATUS_DIMER_LOST_MODE;
        break;
      }

      if (optStatus < 0) {
        status = STATUS_OPTIMIZER_ERROR;
        break;
      }

      double de = objf->getEnergy() - reactantEnergy;
      const double energyLimit = params.saddle_search_options().max_energy;
      if (de > energyLimit) {
        const bool kept = pullUnderEnergy(posBefore);
        if (kept) {
          optim->onTranslationPulledBack();
        }
        if (!kept || objf->getEnergy() - reactantEnergy > energyLimit) {
          status = STATUS_BAD_HIGH_ENERGY;
          // Count the attempt so a first-step abort does not fall through
          // to the iteration == 0 mode estimate below.
          ++iteration;
          break;
        }
        de = objf->getEnergy() - reactantEnergy;
      } else if (translationCap < stockMove) {
        translationCap =
            std::min(stockMove, translationCap * kMaxEnergyCapRegrowth);
        wallStalls = 0;
      }

      // Melander, Laasonen, Jonsson, JCTC 11(3), 1055-1062, 2015
      if (params.saddle_search_options().remove_rotation) {
        eonc::geometry::rotationRemove(pos, matter);
      }
      double stepSize = (matter->pbc(matter->getPositions() - pos)).norm();

      iteration++;

      // Logging
      double eigenval = eonc::eigenmodeGetEigenvalue(*minModeMethod);
      double torque = eonc::eigenmodeStatsTorque(*minModeMethod);
      double angle = eonc::eigenmodeStatsAngle(*minModeMethod);
      long rotations = eonc::eigenmodeStatsRotations(*minModeMethod);
      double conv = objf->getConvergence();

      if (params.saddle_search_options().minmode_method ==
          LowestEigenmode::MINMODE_LANCZOS) {
        QUILL_LOG_DEBUG(
            log,
            "[Lanczos]  {:9} {:9.6f} {:10.4f} {:18.5e} {:9.4f} {:10.6f} "
            "{:7.3f} {:5}\n",
            iteration, stepSize, de, conv, eigenval, torque, angle, rotations);
      } else {
        QUILL_LOG_DEBUG(
            log,
            "[Dimer]  {:9}   {:9.7f}   {:10.4f}   {:18.5e}   {:9.4f}   {:7.3f} "
            "  {:6.3f}   {:4}\n",
            iteration, stepSize, de, conv, eigenval, torque, angle, rotations);
      }

      if (params.debug_options().write_movies || retain_climb_frames_) {
        write_climb_frame(static_cast<uint64_t>(iteration), true, stepSize, de,
                          conv, eigenval, torque, angle, rotations);
      }

      if (params.main_options().checkpoint) {
        if (!eonc::io::io_ok(
                matter->matter2con("displacement_cp.con", false))) {
          QUILL_LOG_WARNING(log, "Failed to write displacement_cp.con");
        }
        eonc::helpers::saveMode("mode_cp.dat", matter,
                                eonc::eigenmodeGetEigenvector(*minModeMethod));
      }

      if (atKnownSaddle()) {
        QUILL_LOG_DEBUG(log, "[MinModeSaddleSearch] reached a known saddle");
        status = STATUS_KNOWN_SADDLE;
        break;
      }

      // Check ImprovedDimer mode convergence
      if (auto *dimer = eonc::asImprovedDimer(*minModeMethod)) {
        if (!dimer->rotationDidConverge) {
          status = (dimer->getEigenvalue() < 0.0) ? STATUS_DIMER_RESTORED_BEST
                                                  : STATUS_DIMER_LOST_MODE;
          if (status == STATUS_DIMER_RESTORED_BEST) {
            QUILL_LOG_DEBUG(log, "Dimer restored to valid state. C_tau={:.4f}",
                            dimer->getEigenvalue());
          }
          break;
        }
      }
    }

    // A rejected or restored mode already carries the curvature needed by
    // the caller; solving it again would discard that handoff state.
    if (iteration == 0 && status != STATUS_DIMER_LOST_MODE &&
        status != STATUS_DIMER_RESTORED_BEST) {
      eonc::eigenmodeCompute(*minModeMethod, matter, mode);
    }

    // Never report STATUS_GOOD when the climb objective is unconverged
    // (issue #20: unfeasible systems must not look like success).
    const bool climbConverged = objf->isConverged();
    const int statusBeforeGuard = status;
    status = finalizeClimbStatus(status, climbConverged);
    if (statusBeforeGuard == STATUS_GOOD && status != STATUS_GOOD) {
      QUILL_LOG_WARNING(log, "[MinModeSaddleSearch] objective not converged; "
                             "refusing STATUS_GOOD");
    }

    long refineForceCalls = 0;
    if (params.saddle_search_options().prfo_refine.enabled &&
        status == STATUS_GOOD && getEigenvalue() < 0.0) {
      const auto refined = eonc::prfoRefineSaddle(*matter, getEigenvector(),
                                                  getEigenvalue(), params);
      refineForceCalls = refined.forceCalls;
      if (!refined.restored) {
        // The mode at the refined point, for mode.dat and the sign test.
        eonc::eigenmodeCompute(*minModeMethod, matter, getEigenvector());
      }
    }

    if (getEigenvalue() > 0.0 && status == STATUS_GOOD) {
      QUILL_LOG_DEBUG(log, "[MinModeSaddleSearch] eigenvalue not negative");
      status = STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE;
    }
    forcecalls =
        static_cast<int>(eonc::eigenmodeTotalForceCalls(*minModeMethod) +
                         probeForceCalls + refineForceCalls);
  }

  return status;
}

void MinModeSaddleSearch::setKnownSaddles(std::vector<AtomMatrix> positions,
                                          std::vector<double> energies) {
  if (positions.size() != energies.size()) {
    throw std::invalid_argument(
        "MinModeSaddleSearch::setKnownSaddles: size mismatch");
  }
  knownSaddlePositions_ = std::move(positions);
  knownSaddleEnergies_ = std::move(energies);
}

bool MinModeSaddleSearch::atKnownSaddle() const {
  if (knownSaddleEnergies_.empty() ||
      !(eonc::eigenmodeGetEigenvalue(*minModeMethod) < 0.0)) {
    return false;
  }
  const double energy = matter->getPotentialEnergy();
  const double eps = params.structure_comparison_options().energy_difference;
  for (size_t k = 0; k < knownSaddleEnergies_.size(); ++k) {
    if (std::abs(energy - knownSaddleEnergies_[k]) > eps) {
      continue;
    }
    Matter known(*matter);
    known.setPositions(knownSaddlePositions_[k]);
    if (known.compare(*matter)) {
      return true;
    }
  }
  return false;
}

double MinModeSaddleSearch::getEigenvalue() {
  return eonc::eigenmodeGetEigenvalue(*minModeMethod);
}

AtomMatrix MinModeSaddleSearch::getEigenvector() {
  return eonc::eigenmodeGetEigenvector(*minModeMethod);
}

} // namespace eonc

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
#include "MinModeSaddleSearch.h"
#include "ConjugateGradients.h"
#include "EigenmodeStrategy.h"
#include "EonLogger.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "ObjectiveFunction.h"
#include "SaddleSearchJob.h"
#include "eonExceptions.hpp"

#include <cmath>
#include <format>
#include <fstream>
#include <memory>
#include <string>

using namespace eonc::helpers;

class MinModeObjectiveFunction : public ObjectiveFunction {
private:
  std::shared_ptr<Matter> matter;
  AtomMatrix eigenvector;
  std::shared_ptr<EigenmodeStrategy> minModeMethod;
  int iteration{0};

public:
  MinModeObjectiveFunction(
      std::shared_ptr<Matter> matterPassed,
      std::shared_ptr<EigenmodeStrategy> minModeMethodPassed,
      AtomMatrix modePassed, const Parameters &paramsPassed)
      : ObjectiveFunction(paramsPassed),
        matter{std::move(matterPassed)},
        minModeMethod{minModeMethodPassed},
        eigenvector{std::move(modePassed)} {}

  ~MinModeObjectiveFunction() override = default;

  VectorXd getGradient(bool fdstep = false) {
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

    AtomMatrix proj = matDot(force, eigenvector) * eigenvector.normalized();

    if (eigenvalue > 0.0) {
      if (params.saddle_search_options.perp_force_ratio > 0.0) {
        double d = params.saddle_search_options.perp_force_ratio;
        force = d * force - (1.0 + d) * proj;
      } else if (params.saddle_search_options.confine_positive.enabled) {
        if (params.saddle_search_options.confine_positive.bowl_breakout) {
          AtomMatrix forceTemp = matter->getForces();
          int nBowlActive =
              params.saddle_search_options.confine_positive.bowl_active;
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
        } else {
          int sufficientForce = 0;
          double minForce =
              params.saddle_search_options.confine_positive.min_force;
          while (sufficientForce <
                 params.saddle_search_options.confine_positive.min_active) {
            sufficientForce = 0;
            force = matter->getForces();
            for (long i = 0; i < matter->numberOfAtoms(); i++) {
              for (int k = 0; k < 3; k++) {
                if (std::abs(force(i, k)) < minForce) {
                  force(i, k) = 0;
                } else {
                  sufficientForce++;
                  force(i, k) =
                      -params.saddle_search_options.confine_positive.boost *
                      proj(i, k);
                }
              }
            }
            minForce *=
                params.saddle_search_options.confine_positive.scale_ratio;
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

  double getEnergy() { return matter->getPotentialEnergy(); }
  void setPositions(const VectorXd &x) { matter->setPositionsV(x); }
  VectorXd getPositions() { return matter->getPositionsV(); }
  int degreesOfFreedom() { return 3 * matter->numberOfAtoms(); }
  bool isConverged() {
    return getConvergence() < params.saddle_search_options.converged_force;
  }

  double getConvergence() {
    if (params.optimizer_options.convergence_metric == "norm") {
      return matter->getForcesFreeV().norm();
    } else if (params.optimizer_options.convergence_metric == "max_atom") {
      return matter->maxForce();
    } else if (params.optimizer_options.convergence_metric == "max_component") {
      return matter->getForces().maxCoeff();
    } else {
      EONC_LOG_CRITICAL("[MinModeSaddleSearch] unknown convergence metric: {}",
                        params.optimizer_options.convergence_metric);
      std::exit(1);
    }
  }

  VectorXd difference(const VectorXd &a, const VectorXd &b) {
    return matter->pbcV(a - b);
  }
};

MinModeSaddleSearch::MinModeSaddleSearch(std::shared_ptr<Matter> matterPassed,
                                         AtomMatrix modePassed,
                                         double reactantEnergyPassed,
                                         const Parameters &parametersPassed,
                                         std::shared_ptr<Potential> potPassed)
    : SaddleSearchMethod(potPassed, parametersPassed),
      matter{matterPassed} {
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
  return run(params.saddle_search_options.max_iterations);
}

int MinModeSaddleSearch::run(long max_iterations_override) {
  long effectiveMaxIter = max_iterations_override;
  QUILL_LOG_DEBUG(
      log, "Saddle point search started from reactant with energy {} eV.",
      reactantEnergy);

  int optStatus;
  bool firstIteration = true;
  const char *forceLabel =
      params.optimizer_options.convergence_metric_label.c_str();

  if (params.saddle_search_options.minmode_method ==
      LowestEigenmode::MINMODE_GPRDIMER) {
    QUILL_LOG_DEBUG(
        log, "================= Using the GP Dimer Library =================");
    eonc::eigenmodeCompute(*minModeMethod, matter, mode);
    if (eonc::eigenmodeGetEigenvalue(*minModeMethod) > 0) {
      QUILL_LOG_DEBUG(log, "GPR eigenvalue: {}",
                      eonc::eigenmodeGetEigenvalue(*minModeMethod));
      return STATUS_NONNEGATIVE_ABORT;
    }
    if (getEigenvalue() > 0.0 && status == STATUS_GOOD) {
      QUILL_LOG_DEBUG(log, "[MinModeSaddleSearch] eigenvalue not negative");
      status = STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE;
    }
    if (std::abs(eonc::eigenmodeGetEigenvalue(*minModeMethod)) <
        params.saddle_search_options.zero_mode_abort_curvature) {
      QUILL_LOG_DEBUG(log, "Zero mode eigenvalue: {}",
                      eonc::eigenmodeGetEigenvalue(*minModeMethod));
      status = STATUS_ZEROMODE_ABORT;
    }
    iteration = eonc::eigenmodeTotalIterations(*minModeMethod);
    forcecalls = eonc::eigenmodeTotalForceCalls(*minModeMethod);
  } else {

    if (params.saddle_search_options.minmode_method ==
        LowestEigenmode::MINMODE_DIMER) {
      QUILL_LOG_INFO(log,
                     "[Dimer]  {:9s}   {:9s}   {:10s}   {:18s}   {:9s}   "
                     "{:7s}   {:6s}   {:4s}   {:5s}\n",
                     "Step", "Step Size", "Delta E", forceLabel, "Curvature",
                     "Torque", "Angle", "Rots", "Align");
    } else if (params.saddle_search_options.minmode_method ==
               LowestEigenmode::MINMODE_LANCZOS) {
      QUILL_LOG_INFO(
          log,
          "[Lanczos]  {:9s} {:9s} {:10s} {:18s} {:9s} {:10s} {:7s} {:5s}\n",
          "Step", "Step Size", "Delta E", forceLabel, "Curvature", "Rel Change",
          "Angle", "Iters");
    } else if (params.saddle_search_options.minmode_method ==
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
        matter, minModeMethod, mode, params);
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
      metadata.scalars.push_back({"rotations",
                                  static_cast<double>(rotations)});
      matter->matter2con(climbLabel, append, &metadata);

      if (params.debug_options.write_deprecated_outs) {
        std::ofstream climbDat(
            climbDatFilename,
            append ? (std::ios::binary | std::ios::app) : std::ios::binary);
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
    if (params.debug_options.write_movies) {
      write_climb_frame(0, false, 0.0, 0.0, objf->getConvergence(),
                        eonc::eigenmodeGetEigenvalue(*minModeMethod), 0.0, 0.0,
                        0);
    }
    if (params.saddle_search_options.nonnegative_displacement_abort) {
      objf->getGradient();
      if (eonc::eigenmodeGetEigenvalue(*minModeMethod) > 0) {
        QUILL_LOG_DEBUG(log, "Nonnegative eigenvalue: {}",
                        eonc::eigenmodeGetEigenvalue(*minModeMethod));
        return STATUS_NONNEGATIVE_ABORT;
      }
    }

    auto optim = eonc::helpers::create::mkOptim(
        objf, params.optimizer_options.method, params);

    while (!objf->isConverged() || iteration == 0) {

      if (!firstIteration) {

        if (params.saddle_search_options.nonlocal_count_abort != 0) {
          long nm = numAtomsMoved(
              initialPosition - matter->getPositions(),
              params.saddle_search_options.nonlocal_distance_abort);
          if (nm >= params.saddle_search_options.nonlocal_count_abort) {
            status = STATUS_NONLOCAL_ABORT;
            break;
          }
        }

        if (std::abs(eonc::eigenmodeGetEigenvalue(*minModeMethod)) <
            params.saddle_search_options.zero_mode_abort_curvature) {
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

      AtomMatrix pos = matter->getPositions();

      try {
        if (params.saddle_search_options.confine_positive.bowl_breakout &&
            eonc::eigenmodeGetEigenvalue(*minModeMethod) > 0 &&
            params.optimizer_options.method == OptType::CG) {
          optStatus = optim->step(-params.optimizer_options.max_move);
        } else {
          optStatus = optim->step(params.optimizer_options.max_move);
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

      // Melander, Laasonen, Jonsson, JCTC 11(3), 1055-1062, 2015
      if (params.saddle_search_options.remove_rotation) {
        rotationRemove(pos, matter);
      }
      double stepSize = (matter->pbc(matter->getPositions() - pos)).norm();

      iteration++;

      // Logging
      double eigenval = eonc::eigenmodeGetEigenvalue(*minModeMethod);
      double torque = eonc::eigenmodeStatsTorque(*minModeMethod);
      double angle = eonc::eigenmodeStatsAngle(*minModeMethod);
      long rotations = eonc::eigenmodeStatsRotations(*minModeMethod);
      double conv = objf->getConvergence();

      if (params.saddle_search_options.minmode_method ==
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

      if (params.debug_options.write_movies) {
        write_climb_frame(static_cast<uint64_t>(iteration), true, stepSize, de,
                          conv, eigenval, torque, angle, rotations);
      }

      if (params.main_options.checkpoint) {
        matter->matter2con("displacement_cp.con", false);
        eonc::helpers::saveMode("mode_cp.dat", matter,
                                eonc::eigenmodeGetEigenvector(*minModeMethod));
      }

      if (de > params.saddle_search_options.max_energy) {
        status = STATUS_BAD_HIGH_ENERGY;
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

    if (iteration == 0) {
      eonc::eigenmodeCompute(*minModeMethod, matter, mode);
    }

    if (getEigenvalue() > 0.0 && status == STATUS_GOOD) {
      QUILL_LOG_DEBUG(log, "[MinModeSaddleSearch] eigenvalue not negative");
      status = STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE;
    }
  }

  return status;
}

double MinModeSaddleSearch::getEigenvalue() {
  return eonc::eigenmodeGetEigenvalue(*minModeMethod);
}

AtomMatrix MinModeSaddleSearch::getEigenvector() {
  return eonc::eigenmodeGetEigenvector(*minModeMethod);
}

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
#include "Dimer.h"
#include "HelperFunctions.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "LowestEigenmode.h"
#include "SaddleSearchJob.h"
#include "eonExceptions.hpp"
#include <memory>
#ifdef WITH_GPRD
#include "AtomicGPDimer.h"
#endif
#include "EpiCenters.h"
#include "ObjectiveFunction.h"

using namespace helper_functions;

class MinModeObjectiveFunction : public ObjectiveFunction {
private:
  AtomMatrix eigenvector;
  std::shared_ptr<LowestEigenmode> minModeMethod;
  int iteration;

public:
  MinModeObjectiveFunction(std::shared_ptr<Matter> matterPassed,
                           std::shared_ptr<LowestEigenmode> minModeMethodPassed,
                           AtomMatrix modePassed,
                           const Parameters &paramsPassed)
      : ObjectiveFunction(matterPassed, paramsPassed),
        minModeMethod{minModeMethodPassed} {
    eigenvector = modePassed;
  }

  ~MinModeObjectiveFunction(void) {}

  VectorXd getGradient(bool fdstep = false) {
    AtomMatrix proj;
    AtomMatrix force = matter->getForces();

    if (!fdstep || iteration == 0) {
      minModeMethod->compute(matter, eigenvector);
      // Check if ImprovedDimer lost the mode immediately after compute()
      auto dimer = std::dynamic_pointer_cast<ImprovedDimer>(minModeMethod);
      if (dimer && !dimer->rotationDidConverge) {
        // Check if it restored to a valid negative curvature state
        if (dimer->getEigenvalue() < 0.0) {
          // Dimer restored to best state - update eigenvector and continue
          // but signal that we should probably stop soon
          eigenvector = minModeMethod->getEigenvector();
          SPDLOG_DEBUG(
              "[MinMode] Dimer restored to best state with C_tau={:.4f}",
              dimer->getEigenvalue());
          throw eonc::DimerModeRestoredException();
        } else {
          // Truly lost - no valid state to restore to
          throw eonc::DimerModeLostException();
        }
      }
      iteration++;
    }

    eigenvector = minModeMethod->getEigenvector();
    double eigenvalue = minModeMethod->getEigenvalue();

    proj =
        (force.array() * eigenvector.array()).sum() * eigenvector.normalized();

    if (0 < eigenvalue) {
      if (params.saddle_search_options.perp_force_ratio > 0.0) {
        // reverse force parallel to eigenvector, and reduce perpendicular force
        double const d = params.saddle_search_options.perp_force_ratio;
        force = d * force - (1. + d) * proj;

        // zero out the smallest forces to keep displacement confined
      } else if (params.saddle_search_options.confine_positive.enabled) {
        if (params.saddle_search_options.confine_positive.bowl_breakout) {
          AtomMatrix forceTemp = matter->getForces();
          double *indices_max;
          indices_max = new double[params.saddle_search_options.confine_positive
                                       .bowl_active];

          // determine the force for the x largest component
          double f_max;
          int i_max;
          for (int j = 0;
               j < params.saddle_search_options.confine_positive.bowl_active;
               j++) {
            f_max = forceTemp.row(0).norm();
            i_max = 0;
            for (int i = 0; i < matter->numberOfAtoms(); i++) {
              if (f_max < forceTemp.row(i).norm()) {
                f_max = forceTemp.row(i).norm();
                i_max = i;
              }
            }
            forceTemp(3 * i_max + 0) = 0;
            forceTemp(3 * i_max + 1) = 0;
            forceTemp(3 * i_max + 2) = 0;
            indices_max[j] = i_max;
          }
          for (int i = 0; i < matter->numberOfAtoms(); i++) {
            forceTemp(3 * i + 0) = 0.0;
            forceTemp(3 * i + 1) = 0.0;
            forceTemp(3 * i + 2) = 0.0;
          }
          // only set the projected forces corresponding to the atoms subject to
          // the largest forces
          for (int j = 0;
               j < params.saddle_search_options.confine_positive.bowl_active;
               j++) {
            size_t tind = static_cast<size_t>(indices_max[j]);
            forceTemp(3 * tind + 0) = -proj(3 * tind + 0);
            forceTemp(3 * tind + 1) = -proj(3 * tind + 1);
            forceTemp(3 * tind + 2) = -proj(3 * tind + 2);
          }
          force = forceTemp;
          delete[] indices_max;
        } else {
          int sufficientForce = 0;
          double minForce =
              params.saddle_search_options.confine_positive.min_force;
          while (sufficientForce <
                 params.saddle_search_options.confine_positive.min_active) {
            sufficientForce = 0;
            force = matter->getForces();
            for (int i = 0; i < 3 * matter->numberOfAtoms(); i++) {
              if (fabs(force(i)) < minForce)
                force(i) = 0;
              else {
                sufficientForce = sufficientForce + 1;
                force(i) =
                    -params.saddle_search_options.confine_positive.boost *
                    proj(i);
              }
            }
            minForce *=
                params.saddle_search_options.confine_positive.scale_ratio;
          }
        }
      } else {
        // follow eigenmode
        force = -proj;
      }
    } else {
      // reversing force parallel to eigenmode
      force += -2. * proj;
    }

    VectorXd forceV = VectorXd::Map(force.data(), 3 * matter->numberOfAtoms());
    return -forceV;
  }
  double getEnergy() { return matter->getPotentialEnergy(); }
  void setPositions(VectorXd x) { matter->setPositionsV(x); }
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
      SPDLOG_DEBUG("[MinModeSaddleSearch] unknown opt_convergence_metric: {}",
                   params.optimizer_options.convergence_metric);
      std::exit(1);
    }
  }

  VectorXd difference(VectorXd a, VectorXd b) { return matter->pbcV(a - b); }
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
  status = STATUS_GOOD;
  iteration = 0;
  log = spdlog::get("combi");

  if (params.saddle_search_options.minmode_method ==
      LowestEigenmode::MINMODE_DIMER) {
    if (params.dimer_options.improved) {
      minModeMethod = std::make_shared<ImprovedDimer>(matter, params, pot);
      if (true) { // TODO(rg): convert to param
        auto dimer = std::static_pointer_cast<ImprovedDimer>(minModeMethod);

        // Convert the input AtomMatrix 'mode' to VectorXd
        VectorXd refVec =
            VectorXd::Map(mode.data(), 3 * matter->numberOfAtoms());
        refVec = refVec.array() * matter->getFreeV().array();

        dimer->setReferenceMode(refVec);
      }
    } else {
      minModeMethod = std::make_shared<Dimer>(matter, params, pot);
    }
  } else if (params.saddle_search_options.minmode_method ==
             LowestEigenmode::MINMODE_LANCZOS) {
    minModeMethod = std::make_shared<Lanczos>(matter, params, pot);
  }
#ifdef WITH_GPRD
  else if (params.saddle_search_options.minmode_method ==
           LowestEigenmode::MINMODE_GPRDIMER) {
    minModeMethod = std::make_shared<AtomicGPDimer>(matter, params, pot);
  }
#endif
}

int MinModeSaddleSearch::run() {
  SPDLOG_LOGGER_DEBUG(
      log, "Saddle point search started from reactant with energy {} eV.",
      reactantEnergy);

  int optStatus;
  int firstIteration = 1;
  const char *forceLabel =
      params.optimizer_options.convergence_metric_label.c_str();

  if (params.saddle_search_options.minmode_method ==
      LowestEigenmode::MINMODE_GPRDIMER) {
    SPDLOG_LOGGER_DEBUG(
        log, "================= Using the GP Dimer Library =================");
    minModeMethod->compute(matter, mode);
    if (minModeMethod->getEigenvalue() > 0) {
      printf("%f\n", minModeMethod->getEigenvalue());
      return STATUS_NONNEGATIVE_ABORT;
    }
    if (getEigenvalue() > 0.0 && status == STATUS_GOOD) {
      SPDLOG_LOGGER_DEBUG(log, "[MinModeSaddleSearch] eigenvalue not negative");
      status = STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE;
    }
    if (fabs(minModeMethod->getEigenvalue()) <
        params.saddle_search_options.zero_mode_abort_curvature) {
      printf("%f\n", minModeMethod->getEigenvalue());
      status = STATUS_ZEROMODE_ABORT;
    }
    // These exist only for the gprdimer
    iteration = minModeMethod->totalIterations;
    forcecalls = minModeMethod->totalForceCalls;
  } else {

    if (params.saddle_search_options.minmode_method ==
        LowestEigenmode::MINMODE_DIMER) {
      SPDLOG_LOGGER_INFO(log,
                         "[Dimer]  {:9s}   {:9s}   {:10s}   {:18s}   {:9s}   "
                         "{:7s}   {:6s}   {:4s}   {:5s}\n",
                         "Step", "Step Size", "Delta E", forceLabel,
                         "Curvature", "Torque", "Angle", "Rots", "Align");
    } else if (params.saddle_search_options.minmode_method ==
               LowestEigenmode::MINMODE_LANCZOS) {
      SPDLOG_LOGGER_INFO(
          log,
          "[Lanczos]  {:9s} {:9s} {:10s} {:18s} {:9s} {:10s} {:7s} {:5s}\n",
          "Step", "Step Size", "Delta E", forceLabel, "Curvature", "Rel Change",
          "Angle", "Iters");
    } else if (params.saddle_search_options.minmode_method ==
               LowestEigenmode::MINMODE_GPRDIMER) {
      SPDLOG_LOGGER_INFO(log,
                         "[GPRDimer]  {:9s}   {:9s}   {:10s}   {:18s}   {:9s}  "
                         " {:7s}   {:6s}   {:4s}\n",
                         "Step", "Step Size", "Delta E", forceLabel,
                         "Curvature", "Torque", "Angle", "Rots");
    }

    ostringstream climb;
    climb << "climb";
    if (params.debug_options.write_movies) {
      matter->matter2con(climb.str(), false);
    }

    AtomMatrix initialPosition = matter->getPositions();

    auto objf = std::make_shared<MinModeObjectiveFunction>(
        matter, minModeMethod, mode, params);
    // objf->getGradient();
    if (params.saddle_search_options.nonnegative_displacement_abort) {
      objf->getGradient();
      if (minModeMethod->getEigenvalue() > 0) {
        printf("%f\n", minModeMethod->getEigenvalue());
        return STATUS_NONNEGATIVE_ABORT;
      }
    }

    auto optim =
        helpers::create::mkOptim(objf, params.optimizer_options.method, params);

    while (!objf->isConverged() || iteration == 0) {

      if (!firstIteration) {

        // Abort if negative mode becomes delocalized
        if (params.saddle_search_options.nonlocal_count_abort != 0) {
          long nm = numAtomsMoved(
              initialPosition - matter->getPositions(),
              params.saddle_search_options.nonlocal_distance_abort);
          if (nm >= params.saddle_search_options.nonlocal_count_abort) {
            status = STATUS_NONLOCAL_ABORT;
            break;
          }
        }

        // Abort if negative mode becomes zero
        // cout << "curvature: "<<fabs(minModeMethod->getEigenvalue())<<"\n";
        if (fabs(minModeMethod->getEigenvalue()) <
            params.saddle_search_options.zero_mode_abort_curvature) {
          printf("%f\n", minModeMethod->getEigenvalue());
          status = STATUS_ZEROMODE_ABORT;
          break;
        }
      }
      firstIteration = 0;

      if (iteration >= params.saddle_search_options.max_iterations) {
        status = STATUS_BAD_MAX_ITERATIONS;
        break;
      }

      AtomMatrix pos = matter->getPositions();

      try {
        if (params.saddle_search_options.confine_positive.bowl_breakout) {
          // use negative step to communicate that the system is the negative
          // region and a max step should be performed
          if ((minModeMethod->getEigenvalue() > 0) and
              (params.optimizer_options.method == OptType::CG)) {
            optStatus = optim->step(-params.optimizer_options.max_move);
          } else {
            optStatus = optim->step(params.optimizer_options.max_move);
          }
        } else {
          optStatus = optim->step(params.optimizer_options.max_move);
        }
      } catch (const eonc::DimerModeRestoredException &e) {
        // Dimer lost mode but restored to valid negative curvature state
        // Check if we're now converged
        SPDLOG_LOGGER_DEBUG(
            log, "Dimer restored to best state.  Checking convergence...");

        // Force might have changed - recompute convergence
        if (objf->isConverged()) {
          SPDLOG_LOGGER_DEBUG(log, "Converged after dimer restoration.");
          status = STATUS_GOOD;
        } else {
          // Not converged, but we have a valid state - report as partial
          // success
          status = STATUS_DIMER_RESTORED_BEST;
        }
        break;
      } catch (const eonc::DimerModeLostException &e) {
        // Truly lost the mode with no valid state
        SPDLOG_LOGGER_WARN(log, "Dimer lost mode completely. Aborting.");
        status = STATUS_DIMER_LOST_MODE;
        break;
      }

      if (optStatus < 0) {
        status = STATUS_OPTIMIZER_ERROR;
        break;
      }

      double de = objf->getEnergy() - reactantEnergy;
      // should be the total displacement of the system not just a single atom
      // double stepSize =
      // helper_functions::maxAtomMotion(matter->pbc(matter->getPositions() -
      // pos));

      double stepSize;

      // Melander, Laasonen, Jonsson, JCTC, 11(3), 1055–1062, 2015
      // http://doi.org/10.1021/ct501155k
      if (params.saddle_search_options.remove_rotation) {
        rotationRemove(pos, matter);
      }
      stepSize = (matter->pbc(matter->getPositions() - pos)).norm();

      iteration++;

      if (params.saddle_search_options.minmode_method ==
          LowestEigenmode::MINMODE_DIMER) {
        SPDLOG_LOGGER_DEBUG(
            log,
            "[Dimer]  {:9}   {:9.7f}   {:10.4f}   {:18.5e}   {:9.4f}   {:7.3f} "
            "  {:6.3f}   {:4}\n",
            iteration, stepSize, matter->getPotentialEnergy() - reactantEnergy,
            objf->getConvergence(), minModeMethod->getEigenvalue(),
            minModeMethod->statsTorque, minModeMethod->statsAngle,
            minModeMethod->statsRotations);
      } else if (params.saddle_search_options.minmode_method ==
                 LowestEigenmode::MINMODE_LANCZOS) {
        SPDLOG_LOGGER_DEBUG(
            log,
            "[Lanczos]  {:9} {:9.6f} {:10.4f} {:18.5e} {:9.4f} {:10.6f} "
            "{:7.3f} {:5}\n",
            iteration, stepSize, matter->getPotentialEnergy() - reactantEnergy,
            objf->getConvergence(), minModeMethod->getEigenvalue(),
            minModeMethod->statsTorque, minModeMethod->statsAngle,
            minModeMethod->statsRotations);
      } else if (params.saddle_search_options.minmode_method ==
                 LowestEigenmode::MINMODE_GPRDIMER) {
        SPDLOG_LOGGER_DEBUG(
            log,
            "[Dimer]  {:9}   {:9.7f}   {:10.4f}   {:18.5e}   {:9.4f}   {:7.3f} "
            "  {:6.3f}   {:4}\n",
            iteration, stepSize, matter->getPotentialEnergy() - reactantEnergy,
            objf->getConvergence(), minModeMethod->getEigenvalue(),
            minModeMethod->statsTorque, minModeMethod->statsAngle,
            minModeMethod->statsRotations);
      } else {
        log = spdlog::get("_traceback");
        SPDLOG_LOGGER_CRITICAL(
            log, "[MinModeSaddleSearch] Unknown min_mode_method: {}",
            params.saddle_search_options.minmode_method);
        std::exit(1);
      }

      if (params.debug_options.write_movies) {
        matter->matter2con(climb.str(), true);
      }

      if (params.main_options.checkpoint) {
        matter->matter2con("displacement_cp.con", false);
        FILE *fileMode = fopen("mode_cp.dat", "wb");
        helper_functions::saveMode(fileMode, matter,
                                   minModeMethod->getEigenvector());
        fclose(fileMode);
      }

      if (de > params.saddle_search_options.max_energy) {
        status = STATUS_BAD_HIGH_ENERGY;
        break;
      }

      if (params.saddle_search_options.minmode_method ==
          LowestEigenmode::MINMODE_DIMER) {
        // We need to cast the pointer to access the ImprovedDimer-specific flag
        auto dimer = std::dynamic_pointer_cast<ImprovedDimer>(minModeMethod);

        if (dimer && !dimer->rotationDidConverge) {
          // This shouldn't happen if we caught the exceptions above,
          // but keep as safety check
          if (dimer->getEigenvalue() < 0.0) {
            SPDLOG_LOGGER_DEBUG(log,
                                "Dimer restored to valid state.  C_tau={:.4f}",
                                dimer->getEigenvalue());
            status = STATUS_DIMER_RESTORED_BEST;
          } else {
            status = STATUS_DIMER_LOST_MODE;
          }
          break;
        }
      }
    }

    if (iteration == 0)
      minModeMethod->compute(matter, mode);

    if (getEigenvalue() > 0.0 && status == STATUS_GOOD) {
      SPDLOG_LOGGER_DEBUG(log, "[MinModeSaddleSearch] eigenvalue not negative");
      status = STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE;
    }
  }

  return status;
}

double MinModeSaddleSearch::getEigenvalue() {
  return minModeMethod->getEigenvalue();
}

AtomMatrix MinModeSaddleSearch::getEigenvector() {
  return minModeMethod->getEigenvector();
}

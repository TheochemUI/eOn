#include "MinModeSaddleSearch.h"
#include "ConjugateGradients.h"
#include "Dimer.h"
#include "HelperFunctions.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "LowestEigenmode.h"
#include "SaddleSearchJob.h"
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
                           std::shared_ptr<Parameters> paramsPassed)
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
      iteration++;
    }

    eigenvector = minModeMethod->getEigenvector();
    double eigenvalue = minModeMethod->getEigenvalue();

    proj =
        (force.array() * eigenvector.array()).sum() * eigenvector.normalized();

    if (0 < eigenvalue) {
      if (params->saddlePerpForceRatio > 0.0) {
        // reverse force parallel to eigenvector, and reduce perpendicular force
        double const d = params->saddlePerpForceRatio;
        force = d * force - (1. + d) * proj;

        // zero out the smallest forces to keep displacement confined
      } else if (params->saddleConfinePositive) {
        if (params->saddleBowlBreakout) {
          AtomMatrix forceTemp = matter->getForces();
          double *indices_max;
          indices_max = new double[params->saddleBowlActive];

          // determine the force for the x largest component
          double f_max;
          int i_max;
          for (int j = 0; j < params->saddleBowlActive; j++) {
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
          for (int j = 0; j < params->saddleBowlActive; j++) {
            size_t tind = static_cast<size_t>(indices_max[j]);
            forceTemp(3 * tind + 0) = -proj(3 * tind + 0);
            forceTemp(3 * tind + 1) = -proj(3 * tind + 1);
            forceTemp(3 * tind + 2) = -proj(3 * tind + 2);
          }
          force = forceTemp;
          delete[] indices_max;
        } else {
          int sufficientForce = 0;
          double minForce = params->saddleConfinePositiveMinForce;
          while (sufficientForce < params->saddleConfinePositiveMinActive) {
            sufficientForce = 0;
            force = matter->getForces();
            for (int i = 0; i < 3 * matter->numberOfAtoms(); i++) {
              if (fabs(force(i)) < minForce)
                force(i) = 0;
              else {
                sufficientForce = sufficientForce + 1;
                force(i) = -params->saddleConfinePositiveBoost * proj(i);
              }
            }
            minForce *= params->saddleConfinePositiveScaleRatio;
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
  bool isConverged() { return getConvergence() < params->saddleConvergedForce; }

  double getConvergence() {
    if (params->optConvergenceMetric == "norm") {
      return matter->getForcesFreeV().norm();
    } else if (params->optConvergenceMetric == "max_atom") {
      return matter->maxForce();
    } else if (params->optConvergenceMetric == "max_component") {
      return matter->getForces().maxCoeff();
    } else {
      SPDLOG_DEBUG("[MinModeSaddleSearch] unknown opt_convergence_metric: {}",
                   params->optConvergenceMetric);
      std::exit(1);
    }
  }

  VectorXd difference(VectorXd a, VectorXd b) { return matter->pbcV(a - b); }
};

MinModeSaddleSearch::MinModeSaddleSearch(
    std::shared_ptr<Matter> matterPassed, AtomMatrix modePassed,
    double reactantEnergyPassed, std::shared_ptr<Parameters> parametersPassed,
    std::shared_ptr<Potential> potPassed)
    : SaddleSearchMethod(potPassed, parametersPassed),
      matter{matterPassed} {
  reactantEnergy = reactantEnergyPassed;
  mode = modePassed;
  status = STATUS_GOOD;
  iteration = 0;
  log = spdlog::get("combi");

  if (params->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
    if (params->dimerImproved) {
      minModeMethod = std::make_shared<ImprovedDimer>(matter, params, pot);
    } else {
      minModeMethod = std::make_shared<Dimer>(matter, params, pot);
    }
  } else if (params->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
    minModeMethod = std::make_shared<Lanczos>(matter, params, pot);
  }
#ifdef WITH_GPRD
  else if (params->saddleMinmodeMethod == LowestEigenmode::MINMODE_GPRDIMER) {
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
  const char *forceLabel = params->optConvergenceMetricLabel.c_str();

  if (params->saddleMinmodeMethod == LowestEigenmode::MINMODE_GPRDIMER) {
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
        params->saddleZeroModeAbortCurvature) {
      printf("%f\n", minModeMethod->getEigenvalue());
      status = STATUS_ZEROMODE_ABORT;
    }
    // These exist only for the gprdimer
    iteration = minModeMethod->totalIterations;
    forcecalls = minModeMethod->totalForceCalls;
  } else {

    if (params->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
      SPDLOG_LOGGER_INFO(log,
                         "[Dimer]  {:9s}   {:9s}   {:10s}   {:18s}   {:9s}   "
                         "{:7s}   {:6s}   {:4s}\n",
                         "Step", "Step Size", "Delta E", forceLabel,
                         "Curvature", "Torque", "Angle", "Rots");
    } else if (params->saddleMinmodeMethod ==
               LowestEigenmode::MINMODE_LANCZOS) {
      SPDLOG_LOGGER_INFO(
          log,
          "[Lanczos]  {:9s} {:9s} {:10s} {:18s} {:9s} {:10s} {:7s} {:5s}\n",
          "Step", "Step Size", "Delta E", forceLabel, "Curvature", "Rel Change",
          "Angle", "Iters");
    } else if (params->saddleMinmodeMethod ==
               LowestEigenmode::MINMODE_GPRDIMER) {
      SPDLOG_LOGGER_INFO(log,
                         "[GPRDimer]  {:9s}   {:9s}   {:10s}   {:18s}   {:9s}  "
                         " {:7s}   {:6s}   {:4s}\n",
                         "Step", "Step Size", "Delta E", forceLabel,
                         "Curvature", "Torque", "Angle", "Rots");
    }

    ostringstream climb;
    climb << "climb";
    if (params->writeMovies) {
      matter->matter2con(climb.str(), false);
    }

    AtomMatrix initialPosition = matter->getPositions();

    auto objf = std::make_shared<MinModeObjectiveFunction>(
        matter, minModeMethod, mode, params);
    // objf->getGradient();
    if (params->saddleNonnegativeDisplacementAbort) {
      objf->getGradient();
      if (minModeMethod->getEigenvalue() > 0) {
        printf("%f\n", minModeMethod->getEigenvalue());
        return STATUS_NONNEGATIVE_ABORT;
      }
    }

    auto optim = helpers::create::mkOptim(objf, params->optMethod, params);

    while (!objf->isConverged() || iteration == 0) {

      if (!firstIteration) {

        // Abort if negative mode becomes delocalized
        if (params->saddleNonlocalCountAbort != 0) {
          long nm = numAtomsMoved(initialPosition - matter->getPositions(),
                                  params->saddleNonlocalDistanceAbort);
          if (nm >= params->saddleNonlocalCountAbort) {
            status = STATUS_NONLOCAL_ABORT;
            break;
          }
        }

        // Abort if negative mode becomes zero
        // cout << "curvature: "<<fabs(minModeMethod->getEigenvalue())<<"\n";
        if (fabs(minModeMethod->getEigenvalue()) <
            params->saddleZeroModeAbortCurvature) {
          printf("%f\n", minModeMethod->getEigenvalue());
          status = STATUS_ZEROMODE_ABORT;
          break;
        }
      }
      firstIteration = 0;

      if (iteration >= params->saddleMaxIterations) {
        status = STATUS_BAD_MAX_ITERATIONS;
        break;
      }

      AtomMatrix pos = matter->getPositions();

      if (params->saddleBowlBreakout) {
        // use negative step to communicate that the system is the negative
        // region and a max step should be performed
        if ((minModeMethod->getEigenvalue() > 0) and
            (params->optMethod == OptType::CG)) {
          optStatus = optim->step(-params->optMaxMove);
        } else {
          optStatus = optim->step(params->optMaxMove);
        }
      } else {
        optStatus = optim->step(params->optMaxMove);
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
      if (params->saddleRemoveRotation) {
        rotationRemove(pos, matter);
      }
      stepSize = (matter->pbc(matter->getPositions() - pos)).norm();

      iteration++;

      if (params->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
        SPDLOG_LOGGER_DEBUG(
            log,
            "[Dimer]  {:9}   {:9.7f}   {:10.4f}   {:18.5e}   {:9.4f}   {:7.3f} "
            "  {:6.3f}   {:4}\n",
            iteration, stepSize, matter->getPotentialEnergy() - reactantEnergy,
            objf->getConvergence(), minModeMethod->getEigenvalue(),
            minModeMethod->statsTorque, minModeMethod->statsAngle,
            minModeMethod->statsRotations);
      } else if (params->saddleMinmodeMethod ==
                 LowestEigenmode::MINMODE_LANCZOS) {
        SPDLOG_LOGGER_DEBUG(
            log,
            "[Lanczos]  {:9} {:9.6f} {:10.4f} {:18.5e} {:9.4f} {:10.6f} "
            "{:7.3f} {:5}\n",
            iteration, stepSize, matter->getPotentialEnergy() - reactantEnergy,
            objf->getConvergence(), minModeMethod->getEigenvalue(),
            minModeMethod->statsTorque, minModeMethod->statsAngle,
            minModeMethod->statsRotations);
      } else if (params->saddleMinmodeMethod ==
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
            params->saddleMinmodeMethod);
        std::exit(1);
      }

      if (params->writeMovies) {
        matter->matter2con(climb.str(), true);
      }

      if (de > params->saddleMaxEnergy) {
        status = STATUS_BAD_HIGH_ENERGY;
        break;
      }

      if (params->checkpoint) {
        matter->matter2con("displacement_cp.con", false);
        FILE *fileMode = fopen("mode_cp.dat", "wb");
        helper_functions::saveMode(fileMode, matter,
                                   minModeMethod->getEigenvector());
        fclose(fileMode);
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

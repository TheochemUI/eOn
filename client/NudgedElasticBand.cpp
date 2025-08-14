#include "NudgedElasticBand.h"
#include "BaseStructures.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "MinModeSaddleSearch.h"
#include "Optimizer.h"
#include "magic_enum/magic_enum.hpp"

#include <filesystem>
#include <spdlog/spdlog.h>

using namespace helper_functions;
namespace fs = std::filesystem;

namespace helper_functions::neb_paths {
std::vector<Matter> linearPath(const Matter &initImg, const Matter &finalImg,
                               const size_t nimgs) {
  std::vector<Matter> all_images_on_path(nimgs + 2, initImg);
  all_images_on_path.front() = Matter(initImg);
  all_images_on_path.back() = Matter(finalImg);
  AtomMatrix posInitial = all_images_on_path.front().getPositions();
  AtomMatrix posFinal = all_images_on_path.back().getPositions();
  AtomMatrix imageSep = initImg.pbc(posFinal - posInitial) / (nimgs + 1);
  // Only the ones which are not the front and back
  for (auto it{std::next(all_images_on_path.begin())};
       it != std::prev(all_images_on_path.end()); ++it) {
    *it = Matter(initImg);
    (*it).setPositions(posInitial +
                       imageSep *
                           int(std::distance(all_images_on_path.begin(), it)));
  }
  return all_images_on_path;
}

std::vector<Matter> filePathInit(const std::vector<fs::path> &fsrcs,
                                 const Matter &refImg, const size_t nimgs) {
  std::vector<Matter> all_images_on_path;
  if (nimgs + 2 != fsrcs.size()) {
    throw std::runtime_error("Error in filePathInit: Expected " +
                             std::to_string(nimgs + 2) + " files, but got " +
                             std::to_string(fsrcs.size()) + ".");
  }
  all_images_on_path.reserve(nimgs + 2);
  // For all images
  for (const auto &filePath : fsrcs) {
    Matter img(refImg);
    img.con2matter(filePath.string());
    all_images_on_path.push_back(img);
  }
  return all_images_on_path;
}

std::vector<fs::path> readFilePaths(const std::string &listFilePath) {
  std::vector<fs::path> paths;
  std::ifstream inputFile(listFilePath);

  if (!inputFile.is_open()) {
    throw std::runtime_error("Error: Could not open path list file: " +
                             listFilePath);
  }

  std::string line;
  while (std::getline(inputFile, line)) {
    // Skip any empty lines in the input file
    if (!line.empty()) {
      paths.emplace_back(line);
    }
  }

  return paths;
}
} // namespace helper_functions::neb_paths

// NEBObjectiveFunction definitions
VectorXd NEBObjectiveFunction::getGradient(bool fdstep) {
  VectorXd forceV;
  forceV.resize(3 * neb->atoms * neb->numImages);
  if (neb->movedAfterForceCall)
    neb->updateForces();
  for (long i = 1; i <= neb->numImages; i++) {
    forceV.segment(3 * neb->atoms * (i - 1), 3 * neb->atoms) =
        VectorXd::Map(neb->projectedForce[i]->data(), 3 * neb->atoms);
  }
  return -forceV;
}

double NEBObjectiveFunction::getEnergy() {
  double Energy{0};
  for (long i = 1; i <= neb->numImages; i++) {
    Energy += neb->path[i]->getPotentialEnergy();
  }
  return Energy;
}

void NEBObjectiveFunction::setPositions(VectorXd x) {
  neb->movedAfterForceCall = true;
  for (long i = 1; i <= neb->numImages; i++) {
    neb->path[i]->setPositions(MatrixXd::Map(
        x.segment(3 * neb->atoms * (i - 1), 3 * neb->atoms).data(), neb->atoms,
        3));
  }
}

VectorXd NEBObjectiveFunction::getPositions() {
  VectorXd posV;
  posV.resize(3 * neb->atoms * neb->numImages);
  for (long i = 1; i <= neb->numImages; i++) {
    posV.segment(3 * neb->atoms * (i - 1), 3 * neb->atoms) =
        VectorXd::Map(neb->path[i]->getPositions().data(), 3 * neb->atoms);
  }
  return posV;
}

int NEBObjectiveFunction::degreesOfFreedom() {
  return 3 * neb->numImages * neb->atoms;
}

bool NEBObjectiveFunction::isUncertain() {
  double maxMaxUnc = std::numeric_limits<double>::lowest();
  double currentMaxUnc{0};
  for (long idx = 0; idx <= neb->numImages; idx++) {
    currentMaxUnc = neb->path[idx]->getEnergyVariance();
    if (currentMaxUnc > maxMaxUnc) {
      maxMaxUnc = currentMaxUnc;
    }
  }
  bool unc_conv{maxMaxUnc > params->gp_uncertainity};
  if (unc_conv) {
    this->status = NudgedElasticBand::NEBStatus::MAX_UNCERTAINITY;
  }
  return unc_conv;
}

bool NEBObjectiveFunction::isConverged() {
  bool force_conv = getConvergence() < params->nebConvergedForce;
  return force_conv;
}

double NEBObjectiveFunction::getConvergence() {
  return neb->convergenceForce();
}

VectorXd NEBObjectiveFunction::difference(VectorXd a, VectorXd b) {
  VectorXd pbcDiff(3 * neb->numImages * neb->atoms);
  for (int i = 1; i <= neb->numImages; i++) {
    int n = (i - 1) * 3 * neb->atoms;
    int m = 3 * neb->atoms;
    pbcDiff.segment(n, m) =
        neb->path[i]->pbcV(a.segment(n, m) - b.segment(n, m));
  }
  return pbcDiff;
}

// Nudged Elastic Band definitions
// First constructor: Now delegates to the second constructor
NudgedElasticBand::NudgedElasticBand(
    std::shared_ptr<Matter> initialPassed, std::shared_ptr<Matter> finalPassed,
    std::shared_ptr<Parameters> parametersPassed,
    std::shared_ptr<Potential> potPassed)
    : NudgedElasticBand(
          [&, initialPassed]() {
            if (!parametersPassed->nebIpath.empty()) {
              std::vector<fs::path> file_paths =
                  helper_functions::neb_paths::readFilePaths(
                      parametersPassed->nebIpath);
              return helper_functions::neb_paths::filePathInit(
                  file_paths, *initialPassed, parametersPassed->nebImages);
            } else {
              return helper_functions::neb_paths::linearPath(
                  *initialPassed, *finalPassed, parametersPassed->nebImages);
            }
          }(),
          parametersPassed, potPassed) {}

// Second constructor: Contains all the common setup code
NudgedElasticBand::NudgedElasticBand(
    std::vector<Matter> initPath, std::shared_ptr<Parameters> parametersPassed,
    std::shared_ptr<Potential> potPassed)
    : params{parametersPassed},
      pot{potPassed} {

  log = spdlog::get("combi");
  this->status = NEBStatus::INIT;
  numImages = params->nebImages;
  atoms = initPath.front().numberOfAtoms();

  // Common initialization logic
  path.resize(numImages + 2);
  tangent.resize(numImages + 2);
  projectedForce.resize(numImages + 2);
  extremumPosition.resize(2 * (numImages + 1));
  extremumEnergy.resize(2 * (numImages + 1));
  extremumCurvature.resize(2 * (numImages + 1));
  numExtrema = 0;

  for (long i = 0; i <= numImages + 1; i++) {
    path[i] = std::make_shared<Matter>(pot, params);
    *path[i] = initPath[i];
    tangent[i] = std::make_shared<AtomMatrix>();
    tangent[i]->resize(atoms, 3);
    projectedForce[i] = std::make_shared<AtomMatrix>();
    projectedForce[i]->resize(atoms, 3);
  }

  // Common final setup
  movedAfterForceCall = true;
  path[0]->getPotentialEnergy();
  path[numImages + 1]->getPotentialEnergy();
  climbingImage = 0;

  // Setup springs
  k_u = params->nebKSPMax;
  k_l = params->nebKSPMin;
  if (params->nebEnergyWeighted) {
    ksp = k_l;
  } else {
    ksp = params->nebSpring;
  }

  // Optional debugging setup
  if (params->estNEBeig) {
    eigenmode_solvers.resize(numImages + 2);
    for (long i = 0; i <= numImages + 1; i++) {
      if (params->nebMMF == LowestEigenmode::MINMODE_DIMER) {
        eigenmode_solvers[i] =
            std::make_shared<ImprovedDimer>(path[i], parametersPassed, pot);
      } else if (params->nebMMF == LowestEigenmode::MINMODE_LANCZOS) {
        eigenmode_solvers[i] =
            std::make_shared<Lanczos>(path[i], parametersPassed, pot);
      } else {
        log = spdlog::get("_traceback");
        SPDLOG_LOGGER_CRITICAL(log, "[Debug] unknown neb_mmf_estimator: {}",
                               params->nebMMF);
        std::exit(1);
      }
    }
  }
}

NudgedElasticBand::NEBStatus NudgedElasticBand::compute(void) {
  long iteration = 0;
  this->status = NEBStatus::RUNNING;

  SPDLOG_LOGGER_DEBUG(log, "Nudged elastic band calculation started.");
  if (params->nebEnergyWeighted) {
    E_ref = std::min(path[0]->getPotentialEnergy(),
                     path[numImages + 1]->getPotentialEnergy());
  }
  updateForces();

  auto objf = std::make_shared<NEBObjectiveFunction>(this, params);

  bool switched{false};
  auto optim = helpers::create::mkOptim(objf, params->optMethod, params);
  std::unique_ptr<Optimizer> refine_optim{nullptr};
  if (params->refineOptMethod != OptType::None) {
    refine_optim =
        helpers::create::mkOptim(objf, params->refineOptMethod, params);
  }
  SPDLOG_DEBUG("{:>10s} {:>12s} {:>14s} {:>11s} {:>12s}", "iteration",
               "step size", params->optConvergenceMetricLabel, "max image",
               "max energy");
  SPDLOG_DEBUG(
      "---------------------------------------------------------------\n");

  while (this->status != NEBStatus::GOOD) {
    if (params->writeMovies) {
      bool append = true;
      if (iteration == 0) {
        append = false;
      }
      path[maxEnergyImage]->matter2con("neb_maximage.con", append);
      std::string nebFilename(fmt::format("neb_path_{:03d}.con", iteration));
      FILE *fileNEBPath = fopen(nebFilename.c_str(), "wb");
      for (long idx = 0; idx <= numImages + 1; idx++) {
        path[idx]->matter2con(fileNEBPath);
      }
      fclose(fileNEBPath);
      printImageData(true, iteration);
    }
    VectorXd pos = objf->getPositions();
    double convForce{convergenceForce()};
    bool originalCIflag = params->nebClimbingImageMethod;
    if (convForce >= params->nebciAfter) {
      params->nebClimbingImageMethod = false;
    }
    if (iteration) {
      // so that we print forces before taking an optimizer step
      if (climbingImage > 0 && climbingImage <= numImages &&
          params->nebciWithMMF && params->nebClimbingImageMethod &&
          convForce < params->nebciMMFAfter && iteration > 1) {

        auto tempMinModeSearch = std::make_shared<MinModeSaddleSearch>(
            path[climbingImage], *tangent[climbingImage],
            path[climbingImage]->getPotentialEnergy(), params, pot);

        // Save the original value of saddleMaxIterations
        auto originalSaddleMaxIterations = params->saddleMaxIterations;
        params->saddleMaxIterations = params->nebciMMFnSteps;
        int minModeStatus = tempMinModeSearch->run();

        // Restore the original value of saddleMaxIterations
        params->saddleMaxIterations = originalSaddleMaxIterations;
        if (minModeStatus != MinModeSaddleSearch::STATUS_GOOD &&
            minModeStatus != MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS) {
          SPDLOG_WARN(
              "MinModeSaddleSearch for climbing image returned status: {}",
              minModeStatus);
        }

        // Now, the climbing image (path[climbingImage]) has been optimized for
        // a few steps by the MinModeSaddleSearch.
        movedAfterForceCall = true;
        updateForces();
      }

      if (iteration >= params->nebMaxIterations) {
        status = NEBStatus::BAD_MAX_ITERATIONS;
        break;
      }
      if (refine_optim) {
        if (optim && convForce > params->refineThreshold) {
          optim->step(params->optMaxMove);
        } else {
          if (!switched) {
            switched = true;
            SPDLOG_DEBUG("Switched to {}", magic_enum::enum_name<OptType>(
                                               params->refineOptMethod));
          }
          refine_optim->step(params->optMaxMove);
        }
      } else {
        optim->step(params->optMaxMove);
      }
    }
    iteration++;
    params->nebClimbingImageMethod = originalCIflag;

    double dE = path[maxEnergyImage]->getPotentialEnergy() -
                path[0]->getPotentialEnergy();
    double stepSize = helper_functions::maxAtomMotionV(
        path[0]->pbcV(objf->getPositions() - pos));
    SPDLOG_LOGGER_DEBUG(log, "{:>10} {:>12.4e} {:>14.4e} {:>11} {:>12.4}",
                        iteration, stepSize, convergenceForce(), maxEnergyImage,
                        dE);

    if (pot->getType() == PotType::CatLearn) {
      if (objf->isUncertain()) {
        SPDLOG_LOGGER_DEBUG(log, "NEB failed due to high uncertainity");
        status = NEBStatus::MAX_UNCERTAINITY;
        break;
      } else if (objf->isConverged()) {
        SPDLOG_LOGGER_DEBUG(log, "NEB converged\n");
        status = NEBStatus::GOOD;
        break;
      }

    } else {
      if (objf->isConverged()) {
        SPDLOG_LOGGER_DEBUG(log, "NEB converged\n");
        status = NEBStatus::GOOD;
        break;
      }
    }
  }
  return status;
}

// generate the force value that is compared to the convergence criterion
double NudgedElasticBand::convergenceForce(void) {
  if (movedAfterForceCall)
    updateForces();
  double fmax = 0;

  for (long i = 1; i <= numImages; i++) {

    if (params->nebClimbingImageConvergedOnly == true &&
        params->nebClimbingImageMethod && climbingImage != 0) {
      i = climbingImage;
    }
    if (params->optConvergenceMetric == "norm") {
      fmax = max(fmax, projectedForce[i]->norm());
    } else if (params->optConvergenceMetric == "max_atom") {
      for (int j = 0; j < path[0]->numberOfAtoms(); j++) {
        if (path[0]->getFixed(j))
          continue;
        fmax = max(fmax, projectedForce[i]->row(j).norm());
      }
    } else if (params->optConvergenceMetric == "max_component") {
      fmax = max(fmax, projectedForce[i]->maxCoeff());
    } else {
      log = spdlog::get("_traceback");
      SPDLOG_LOGGER_CRITICAL(
          log, "[Nudged Elastic Band] unknown opt_convergence_metric: {}",
          params->optConvergenceMetric);
      std::exit(1);
    }
    if (params->nebClimbingImageConvergedOnly == true &&
        params->nebClimbingImageMethod && climbingImage != 0) {
      break;
    }
  }
  return fmax;
}

// Update the forces, do the projections, and add spring forces
void NudgedElasticBand::updateForces(void) {
  // variables for tangent
  double maxDiffEnergy, minDiffEnergy;
  double energyDiffPrev, energyDiffNext;
  double energy, energyPrev, energyNext;
  // bool higherEnergyPrev, higherEnergyNext;

  // variables for climbing image
  double maxEnergy;

  // variables for the energy weighted springs
  k_l = params->nebKSPMin;
  k_u = params->nebKSPMax;
  std::vector<double> springConstants(numImages + 1, k_l);

  // variables for force projections
  AtomMatrix force(atoms, 3), forcePerp(atoms, 3), forcePar(atoms, 3);
  AtomMatrix forceSpringPar(atoms, 3), forceSpring(atoms, 3),
      forceSpringPerp(atoms, 3);
  AtomMatrix forceDNEB(atoms, 3);
  AtomMatrix pos(atoms, 3), posNext(atoms, 3), posPrev(atoms, 3),
      posDiffNext(atoms, 3), posDiffPrev(atoms, 3);
  double distNext, distPrev;

  // update the forces on the numImages and find the highest energy image
  maxEnergy = path[0]->getPotentialEnergy();
  maxEnergyImage = 0;
  for (long i = 1; i <= numImages + 1; i++) {
    path[i]->getForces();
    if (path[i]->getPotentialEnergy() > maxEnergy) {
      maxEnergy = path[i]->getPotentialEnergy();
      maxEnergyImage = i;
    }
  }

  // Energy weighted springs, calculated here since all the springs are used
  // internally
  if (params->nebEnergyWeighted) {
    for (int idx = 1; idx <= numImages + 1; idx++) {
      double Ei = std::max(path[idx]->getPotentialEnergy(),
                           path[idx - 1]->getPotentialEnergy());
      if (Ei > E_ref) {
        double alpha_i = (maxEnergy - Ei) / (maxEnergy - E_ref);
        springConstants[idx - 1] =
            (1 - alpha_i) * k_u + alpha_i * k_l; // Equation (3) and (4)
      } // else always k_l
    }
  }

  for (long i = 1; i <= numImages; i++) {
    // set local variables
    force = path[i]->getForces();
    pos = path[i]->getPositions();
    posPrev = path[i - 1]->getPositions();
    posNext = path[i + 1]->getPositions();
    energy = path[i]->getPotentialEnergy();
    energyPrev = path[i - 1]->getPotentialEnergy();
    energyNext = path[i + 1]->getPotentialEnergy();
    posDiffNext = path[i]->pbc(posNext - pos); // R[i+1] - R[i]
    posDiffPrev = path[i]->pbc(pos - posPrev); // R[i] - R[i-1]
    distNext = posDiffNext.norm();             // Distance to next image
    distPrev = posDiffPrev.norm();             // Distance to previous image

    // determine the tangent
    if (params->nebOldTangent) {
      // old tangent
      *tangent[i] = posDiffNext;
    } else {
      // new improved tangent
      // higherEnergyPrev = energyPrev > energyNext;
      // higherEnergyNext = energyNext > energyPrev;

      if (energyNext > energy && energy > energyPrev) {
        *tangent[i] = posDiffNext;
      } else if (energy > energyNext && energyPrev > energy) {
        *tangent[i] = posDiffPrev;
      } else {
        // we are at an extremum
        energyDiffPrev = energyPrev - energy;
        energyDiffNext = energyNext - energy;

        // calculate the energy difference to neighboring numImages
        minDiffEnergy = min(abs(energyDiffPrev), abs(energyDiffNext));
        maxDiffEnergy = max(abs(energyDiffPrev), abs(energyDiffNext));

        // use these energy differences to weight the tangent
        if (energyDiffPrev > energyDiffNext) {
          *tangent[i] = posDiffNext * minDiffEnergy;
          *tangent[i] += posDiffPrev * maxDiffEnergy;
        } else {
          *tangent[i] = posDiffNext * maxDiffEnergy;
          *tangent[i] += posDiffPrev * minDiffEnergy;
        }
      }
    }
    tangent[i]->normalize();

    // project the forces and add springs
    force = path[i]->getForces();

    // calculate the force perpendicular to the tangent
    forcePerp =
        force - (force.array() * (*tangent[i]).array()).sum() * *tangent[i];
    if (params->nebEnergyWeighted) {
      // Spring for segment (i) -> (i+1)
      double kspNext = springConstants[i];
      // Spring for segment (i-1) -> (i)
      double kspPrev = springConstants[i - 1];
      forceSpringPar =
          ((kspNext * distNext) - (kspPrev * distPrev)) * *tangent[i];
    } else {
      this->ksp = params->nebSpring;
      forceSpringPar = this->ksp * (distNext - distPrev) * *tangent[i];
      forceSpring = this->ksp * path[i]->pbc((posNext - pos) - (pos - posPrev));
    }

    if (params->nebDoublyNudged) {
      if (not params->nebEnergyWeighted) {
        forceSpringPerp =
            forceSpring -
            (forceSpring.array() * (*tangent[i]).array()).sum() * *tangent[i];
        forceDNEB =
            forceSpringPerp -
            (forceSpringPerp.array() * forcePerp.normalized().array()).sum() *
                forcePerp.normalized();
        if (params->nebDoublyNudgedSwitching) {
          double switching;
          switching =
              2.0 / M_PI *
              atan(pow(forcePerp.norm(), 2) / pow(forceSpringPerp.norm(), 2));
          forceDNEB *= switching;
        }
      } else {
        SPDLOG_WARN("Not using doubly nudged since energy_weighted is set");
      }
    } else {
      forceDNEB.setZero();
    }

    if (params->nebClimbingImageMethod && i == maxEnergyImage) {
      // we are at the climbing image
      climbingImage = maxEnergyImage;
      *projectedForce[i] =
          force -
          (2.0 * (force.array() * (*tangent[i]).array()).sum() * *tangent[i]) +
          forceDNEB;
    } else // all non-climbing numImages
    {
      // sum the spring and potential forces for the neb force
      if (params->nebElasticBand) {
        if (not params->nebEnergyWeighted) {
          *projectedForce[i] = forceSpring + force;
        } else {
          SPDLOG_WARN("Not using elastic_band  since energy_weighted is set");
        }
      } else {
        *projectedForce[i] = forceSpringPar + forcePerp + forceDNEB;
      }
      //*projectedForce[i] = forceSpring + forcePerp;

      // if (params->nebFullSpring) {

      movedAfterForceCall = false; // so that we don't repeat a force call
    }

    // zero net translational force
    if (path[i]->numberOfFreeAtoms() == path[i]->numberOfAtoms()) {
      for (int j = 0; j <= 2; j++) {
        double translationMag = projectedForce[i]->col(j).sum();
        int natoms = projectedForce[i]->col(j).size();
        projectedForce[i]->col(j).array() -= translationMag / ((double)natoms);
      }
    }
  }

  return;
}

// Print NEB image data
void NudgedElasticBand::printImageData(bool writeToFile, size_t idx) {
  double dist, distTotal = 0;
  AtomMatrix tangentStart =
      path[0]->pbc(path[1]->getPositions() - path[0]->getPositions());
  AtomMatrix tangentEnd = path[numImages]->pbc(
      path[numImages + 1]->getPositions() - path[numImages]->getPositions());
  AtomMatrix tang;
  std::string header, fmttr;
  if (params->estNEBeig) {
    header = fmt::format("{:>3s} {:>12s} {:>12s} {:>12s} {:>12s}", "img",
                         "rxn_coord", "energy", "f_para", "eigval");
    fmttr = "{:>3} {:>12.6f} {:>12.6f} {:>12.6f} {:>12.6f}";
  } else {
    header = fmt::format("{:>3s} {:>12s} {:>12s} {:>12s}", "img", "rxn_coord",
                         "energy", "f_para");
    fmttr = "{:>3} {:>12.6f} {:>12.6f} {:>12.6f}";
  }

  // Endpoint tangents must be normalized for a correct projection
  tangentStart.normalize();
  tangentEnd.normalize();

  std::shared_ptr<spdlog::logger> fileLogger;
  if (writeToFile) {
    if (spdlog::get("file_logger")) {
      spdlog::drop("file_logger");
    }
    auto neb_dat_fs = fmt::format("neb_{:03}.dat", idx);
    if (fs::exists(neb_dat_fs)) {
      fs::remove(neb_dat_fs);
    }
    fileLogger = spdlog::basic_logger_mt("file_logger", neb_dat_fs);
    fileLogger->set_pattern("%v");
    fileLogger->info(header);
  }
  const double energy_reactant = path[0]->getPotentialEnergy();

  for (long i = 0; i <= numImages + 1; i++) {
    if (i == 0) {
      tang = tangentStart;
    } else if (i == numImages + 1) {
      tang = tangentEnd;
    } else {
      tang = *tangent[i];
    }

    if (i > 0) {
      dist = path[i]->distanceTo(*path[i - 1]);
      distTotal += dist;
    }

    // Calculate standard values for logging
    double relative_energy = path[i]->getPotentialEnergy() - energy_reactant;
    double parallel_force = (path[i]->getForces().array() * tang.array()).sum();

    // Conditionally run an MMF and log the eigenvalue
    if (params->estNEBeig) {
      eigenmode_solvers[i]->compute(path[i], tang);
      double lowest_eigenvalue = eigenmode_solvers[i]->getEigenvalue();
      if (fileLogger) {
        fileLogger->info(fmttr, i, distTotal, relative_energy, parallel_force,
                         lowest_eigenvalue);
      } else {
        SPDLOG_LOGGER_DEBUG(log, fmttr, i, distTotal, relative_energy,
                            parallel_force, lowest_eigenvalue);
      }
    } else { // Standard output without the eigenvalue
      if (fileLogger) {
        fileLogger->info(fmttr, i, distTotal, relative_energy, parallel_force);
      } else {
        SPDLOG_LOGGER_DEBUG(log, fmttr, i, distTotal, relative_energy,
                            parallel_force);
      }
    }
  }

  if (fileLogger) {
    spdlog::drop("file_logger");
  }
}

// Estimate the barrier using a cubic spline
void NudgedElasticBand::findExtrema(void) {
  // calculate the cubic parameters for each interval (a,b,c,d)

  AtomMatrix tangentEndpoint;
  double a[numImages + 1], b[numImages + 1], c[numImages + 1], d[numImages + 1];
  double F1, F2, U1, U2, dist;

  for (long i = 0; i <= numImages; i++) {
    dist = path[i]->distanceTo(*path[i + 1]);
    if (i == 0) {
      tangentEndpoint =
          path[i]->pbc(path[1]->getPositions() - path[0]->getPositions());
      tangentEndpoint.normalize();
      F1 =
          (path[i]->getForces().array() * tangentEndpoint.array()).sum() * dist;
    } else {
      F1 = (path[i]->getForces().array() * (*tangent[i]).array()).sum() * dist;
    }
    if (i == numImages) {
      tangentEndpoint = path[i + 1]->pbc(path[numImages + 1]->getPositions() -
                                         path[numImages]->getPositions());
      tangentEndpoint.normalize();
      F2 = (path[i + 1]->getForces().array() * tangentEndpoint.array()).sum() *
           dist;
    } else {
      F2 =
          (path[i + 1]->getForces().array() * (*tangent[i + 1]).array()).sum() *
          dist;
    }
    U1 = path[i]->getPotentialEnergy();
    U2 = path[i + 1]->getPotentialEnergy();
    a[i] = U1;
    b[i] = -F1;
    c[i] = 3. * (U2 - U1) + 2. * F1 + F2;
    d[i] = -2. * (U2 - U1) - (F1 + F2);
  }

  // finding extrema along the MEP

  //    long numExtrema = 0;
  //    double extremaEnergy[2*(numImages+1)]; // the maximum number of
  //    extrema double extremaPosition[2*(numImages+1)];
  double discriminant, f;

  for (long i = 0; i <= numImages; i++) {
    discriminant = pow(c[i], 2) - 3. * b[i] * d[i];
    if (discriminant >= 0) {
      f = -1;

      // quadratic case
      if ((d[i] == 0) && (c[i] != 0)) {
        f = (-b[i] / (2. * c[i]));
      }
      // cubic case 1
      else if (d[i] != 0) {
        f = -(c[i] + sqrt(discriminant)) / (3. * d[i]);
      }
      if ((f >= 0) && (f <= 1)) {
        extremumPosition[numExtrema] = i + f;
        extremumEnergy[numExtrema] =
            d[i] * pow(f, 3) + c[i] * pow(f, 2) + b[i] * f + a[i];
        extremumCurvature[numExtrema] = 6.0 * d[i] * f + 2 * c[i];
        numExtrema++;
      }
      // cubic case 2
      if (d[i] != 0) {
        f = (-(c[i] - sqrt(discriminant)) / (3. * d[i]));
      }
      if ((f >= 0) && (f <= 1)) {
        extremumPosition[numExtrema] = i + f;
        extremumEnergy[numExtrema] =
            d[i] * pow(f, 3) + c[i] * pow(f, 2) + b[i] * f + a[i];
        extremumCurvature[numExtrema] = 6 * d[i] * f + 2 * c[i];
        numExtrema++;
      }
    }
  }

  SPDLOG_LOGGER_DEBUG(log, "Found {} extrema", numExtrema);
  SPDLOG_LOGGER_DEBUG(log, "Energy reference: {}",
                      path[0]->getPotentialEnergy());
  for (long i = 0; i < numExtrema; i++) {
    SPDLOG_LOGGER_DEBUG(
        log, "extrema #{} at image position {} with energy {} and curvature {}",
        i + 1, extremumPosition[i],
        extremumEnergy[i] - path[0]->getPotentialEnergy(),
        extremumCurvature[i]);
  }
}

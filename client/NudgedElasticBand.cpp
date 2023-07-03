#include "NudgedElasticBand.h"
#include "Log.h"
#include "Optimizer.h"

using namespace helper_functions;

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
    (*it).setPositions(
        posInitial +
        imageSep * double(std::distance(all_images_on_path.begin(), it)));
  }
  return all_images_on_path;
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
  double Energy = 0;
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

bool NEBObjectiveFunction::isConverged() {
  return getConvergence() < params->nebConvergedForce;
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
NudgedElasticBand::NudgedElasticBand(
    std::shared_ptr<Matter> initialPassed, std::shared_ptr<Matter> finalPassed,
    std::shared_ptr<Parameters> parametersPassed,
    std::shared_ptr<Potential> potPassed)
    : params{parametersPassed}, pot{potPassed} {
  numImages = params->nebImages;
  atoms = initialPassed->numberOfAtoms();
  auto linear_path = helper_functions::neb_paths::linearPath(
      *initialPassed, *finalPassed, params->nebImages);
  path.resize(numImages + 2);
  tangent.resize(numImages + 2);
  projectedForce.resize(numImages + 2);
  extremumPosition.resize(2 * (numImages + 1));
  extremumEnergy.resize(2 * (numImages + 1));
  extremumCurvature.resize(2 * (numImages + 1));
  numExtrema = 0;

  log("\nNEB: initialize\n");
  for (long i = 0; i <= numImages + 1; i++) {
    path[i] = std::make_shared<Matter>(pot, params);
    *path[i] = linear_path[i];
    tangent[i] = std::make_shared<AtomMatrix>();
    tangent[i]->resize(atoms, 3);
    projectedForce[i] = std::make_shared<AtomMatrix>();
    projectedForce[i]->resize(atoms, 3);
  }
  *path[numImages + 1] = *finalPassed; // final image

  movedAfterForceCall = true;

  // Make sure that the endpoints know their energy
  path[0]->getPotentialEnergy();
  path[numImages + 1]->getPotentialEnergy();
  climbingImage = 0;

  return;
}

NudgedElasticBand::NEBStatus NudgedElasticBand::compute(void) {
  NudgedElasticBand::NEBStatus status = NEBStatus::STATUS_INIT;
  long iteration = 0;

  log("Nudged elastic band calculation started.\n");

  updateForces();

  NEBObjectiveFunction objf(this, params);

  Optimizer *optimizer = Optimizer::getOptimizer(&objf, params.get());

  const char *forceLabel = params->optConvergenceMetricLabel.c_str();
  log("%10s %12s %14s %11s %12s\n", "iteration", "step size", forceLabel,
      "max image", "max energy");
  log("---------------------------------------------------------------\n");

  char fmt[] = "%10li %12.4e %14.4e %11li %12.4f\n";
  char fmtTiny[] = "%10li %12.4e %14.4e %11li %12.4e\n";

  while (!objf.isConverged()) {
    if (params->writeMovies) {
      bool append = true;
      if (iteration == 0)
        append = false;
      path[maxEnergyImage]->matter2con("neb_maximage.con", append);
    }
    VectorXd pos = objf.getPositions();
    if (iteration) { // so that we print forces before taking an optimizer step
      if (iteration >= params->nebMaxIterations) {
        status = NEBStatus::STATUS_BAD_MAX_ITERATIONS;
        break;
      }
      optimizer->step(params->optMaxMove);
    }
    iteration++;

    double dE = path[maxEnergyImage]->getPotentialEnergy() -
                path[0]->getPotentialEnergy();
    double stepSize = helper_functions::maxAtomMotionV(
        path[0]->pbcV(objf.getPositions() - pos));
    if (dE > 0.01) {
      log(fmt, iteration, stepSize, convergenceForce(), maxEnergyImage, dE);
    } else {
      log(fmtTiny, iteration, stepSize, convergenceForce(), maxEnergyImage, dE);
    }
  }

  if (objf.isConverged()) {
    status = NEBStatus::STATUS_GOOD;
    log("NEB converged\n");
  }

  printImageData();
  findExtrema();

  delete optimizer;
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
      log("[Nudged Elastic Band] unknown opt_convergence_metric: %s\n",
          params->optConvergenceMetric.c_str());
      exit(1);
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

  // variables for force projections
  AtomMatrix force(atoms, 3), forcePerp(atoms, 3), forcePar(atoms, 3);
  AtomMatrix forceSpringPar(atoms, 3), forceSpring(atoms, 3),
      forceSpringPerp(atoms, 3);
  AtomMatrix forceDNEB(atoms, 3);
  AtomMatrix pos(atoms, 3), posNext(atoms, 3), posPrev(atoms, 3);
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

  for (long i = 1; i <= numImages; i++) {
    // set local variables
    force = path[i]->getForces();
    pos = path[i]->getPositions();
    posPrev = path[i - 1]->getPositions();
    posNext = path[i + 1]->getPositions();
    energy = path[i]->getPotentialEnergy();
    energyPrev = path[i - 1]->getPotentialEnergy();
    energyNext = path[i + 1]->getPotentialEnergy();

    // determine the tangent
    if (params->nebOldTangent) {
      // old tangent
      *tangent[i] = path[i]->pbc(posNext - posPrev);
    } else {
      // new improved tangent
      // higherEnergyPrev = energyPrev > energyNext;
      // higherEnergyNext = energyNext > energyPrev;

      if (energyNext > energy && energy > energyPrev) {
        *tangent[i] = path[i]->pbc(posNext - pos);
      } else if (energy > energyNext && energyPrev > energy) {
        *tangent[i] = path[i]->pbc(pos - posPrev);
      } else {
        // we are at an extremum
        energyDiffPrev = energyPrev - energy;
        energyDiffNext = energyNext - energy;

        // calculate the energy difference to neighboring numImages
        minDiffEnergy = min(abs(energyDiffPrev), abs(energyDiffNext));
        maxDiffEnergy = max(abs(energyDiffPrev), abs(energyDiffNext));

        // use these energy differences to weight the tangent
        if (energyDiffPrev > energyDiffNext) {
          *tangent[i] = path[i]->pbc(posNext - pos) * minDiffEnergy;
          *tangent[i] += path[i]->pbc(pos - posPrev) * maxDiffEnergy;
        } else {
          *tangent[i] = path[i]->pbc(posNext - pos) * maxDiffEnergy;
          *tangent[i] += path[i]->pbc(pos - posPrev) * minDiffEnergy;
        }
      }
    }
    tangent[i]->normalize();

    // project the forces and add springs
    force = path[i]->getForces();

    // calculate the force perpendicular to the tangent
    forcePerp =
        force - (force.array() * (*tangent[i]).array()).sum() * *tangent[i];
    forceSpring =
        params->nebSpring * path[i]->pbc((posNext - pos) - (pos - posPrev));

    // calculate the spring force
    distPrev = path[i]->pbc(posPrev - pos).squaredNorm();
    distNext = path[i]->pbc(posNext - pos).squaredNorm();
    forceSpringPar = params->nebSpring * (distNext - distPrev) * *tangent[i];

    if (params->nebDoublyNudged) {
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
        *projectedForce[i] = forceSpring + force;
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
void NudgedElasticBand::printImageData(bool writeToFile) {
  double dist, distTotal = 0;
  AtomMatrix tangentStart =
      path[0]->pbc(path[1]->getPositions() - path[0]->getPositions());
  AtomMatrix tangentEnd = path[numImages]->pbc(
      path[numImages + 1]->getPositions() - path[numImages]->getPositions());
  AtomMatrix tang;

  log("Image data (as in neb.dat)\n");

  FILE *fh = NULL;
  if (writeToFile) {
    fh = fopen("neb.dat", "w");
  }

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
    if (fh == NULL) {
      log("%3li %12.6f %12.6f %12.6f\n", i, distTotal,
          path[i]->getPotentialEnergy() - path[0]->getPotentialEnergy(),
          (path[i]->getForces().array() * tang.array()).sum());
    } else {
      fprintf(fh, "%3li %12.6f %12.6f %12.6f\n", i, distTotal,
              path[i]->getPotentialEnergy() - path[0]->getPotentialEnergy(),
              (path[i]->getForces().array() * tang.array()).sum());
    }
  }
  if (writeToFile) {
    fclose(fh);
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
  //    double extremaEnergy[2*(numImages+1)]; // the maximum number of extrema
  //    double extremaPosition[2*(numImages+1)];
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

  log("\nFound %li extrema\n", numExtrema);
  log("Energy reference: %f\n", path[0]->getPotentialEnergy());
  for (long i = 0; i < numExtrema; i++) {
    log("extrema #%li at image position %f with energy %f and curvature %f\n",
        i + 1, extremumPosition[i],
        extremumEnergy[i] - path[0]->getPotentialEnergy(),
        extremumCurvature[i]);
  }
}

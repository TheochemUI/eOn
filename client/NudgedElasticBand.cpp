#include "NudgedElasticBand.h"

#include "Log.h"
#include "Optimizer.h"

#include <algorithm>
#include <filesystem>
#include <fstream>

#include <fmt/core.h>
#include <fmt/printf.h>

using namespace helper_functions;
namespace fs = std::filesystem;

int NudgedElasticBand::compute() {
    nebStatus status{nebStatus::STATUS_INIT};
    size_t iteration = 0;

    log("Nudged elastic band calculation started.\n");

    updateForces();

    NEBObjectiveFunction objf(this, parameters);

    Optimizer *optimizer = Optimizer::getOptimizer(&objf, parameters);

    const char *forceLabel = parameters->optConvergenceMetricLabel.c_str();
    log("%10s %12s %14s %11s %12s\n",
        "iteration",
        "step size",
        forceLabel,
        "max image",
        "max energy");
    log("---------------------------------------------------------------\n");

    char fmt[] = "%10li %12.4e %14.4e %11li %12.4f\n";
    char fmtTiny[] = "%10li %12.4e %14.4e %11li %12.4e\n";

    while (!objf.isConverged()) {
        if (parameters->writeMovies) {
            bool append = true;
            if (iteration == 0)
                append = false;
            neb_images[maxEnergyImage].matter2con("neb_maximage.con", append);
        }
        VectorXd pos = objf.getPositions();
        if (iteration) { // so that we print forces before taking an optimizer step
            if (iteration >= parameters->nebMaxIterations) {
                status = nebStatus::STATUS_BAD_MAX_ITERATIONS;
                break;
            }
            optimizer->step(parameters->optMaxMove);
        }
        iteration++;

        double dE
            = neb_images[maxEnergyImage].getPotentialEnergy() - neb_images[0].getPotentialEnergy();
        double stepSize
            = helper_functions::maxAtomMotionV(neb_images[0].pbcV(objf.getPositions() - pos));
        if (dE > 0.01) {
            log(fmt, iteration, stepSize, convergenceForce(), maxEnergyImage, dE);
        } else {
            log(fmtTiny, iteration, stepSize, convergenceForce(), maxEnergyImage, dE);
        }
    }

    if (objf.isConverged()) {
        status = nebStatus::STATUS_GOOD;
        log("NEB converged\n");
    }

    printImageData();
    findExtrema();

    delete optimizer;
    return static_cast<int>(status);
}

// generate the force value that is compared to the convergence criterion
double NudgedElasticBand::convergenceForce(void) {
    if (movedAfterForceCall)
        updateForces();
    double fmax = 0;

    for (size_t i{1}; i <= nimages; i++) {

        if (parameters->nebClimbingImageConvergedOnly == true && parameters->nebClimbingImageMethod
            && climbingImage != 0) {
            i = climbingImage;
        }
        if (parameters->optConvergenceMetric == "norm") {
            fmax = max(fmax, projectedForce[i].norm());
        } else if (parameters->optConvergenceMetric == "max_atom") {
            for (int j = 0; j < neb_images[0].numberOfAtoms(); j++) {
                if (neb_images[0].getFixed(j))
                    continue;
                fmax = max(fmax, projectedForce[i].row(j).norm());
            }
        } else if (parameters->optConvergenceMetric == "max_component") {
            fmax = max(fmax, projectedForce[i].maxCoeff());
        } else {
            log("[Nudged Elastic Band] unknown opt_convergence_metric: %s\n",
                parameters->optConvergenceMetric.c_str());
            exit(1);
        }
        if (parameters->nebClimbingImageConvergedOnly == true && parameters->nebClimbingImageMethod
            && climbingImage != 0) {
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
    AtomMatrix forceSpringPar(atoms, 3), forceSpring(atoms, 3), forceSpringPerp(atoms, 3);
    AtomMatrix forceDNEB(atoms, 3);
    AtomMatrix pos(atoms, 3), posNext(atoms, 3), posPrev(atoms, 3);
    double distNext, distPrev;

    // update the forces on the nimages and find the highest energy image
    auto comparer
        = [](Matter &m1, Matter &m2) { return m1.getPotentialEnergy() < m2.getPotentialEnergy(); };

    // maxEnergy = std::max(neb_images.begin(), neb_images.end(), comparer);
    maxEnergyImage = std::distance(
        neb_images.begin(), std::max_element(neb_images.begin(), neb_images.end(), comparer));
    maxEnergy = neb_images[maxEnergyImage].getPotentialEnergy();

    // maxEnergy = neb_images[0].getPotentialEnergy();
    // maxEnergyImage = 0;
    // for (size_t i = 1; i <= nimages + 1; i++) {
    //     neb_images[i].getForces();
    //     if (neb_images[i].getPotentialEnergy() > maxEnergy) {
    //         maxEnergy = neb_images[i].getPotentialEnergy();
    //         maxEnergyImage = i;
    //     }
    // }

    for (size_t i{1}; i <= nimages; i++) {
        // set local variables
        force = neb_images[i].getForces();
        pos = neb_images[i].getPositions();
        posPrev = neb_images[i - 1].getPositions();
        posNext = neb_images[i + 1].getPositions();
        energy = neb_images[i].getPotentialEnergy();
        energyPrev = neb_images[i - 1].getPotentialEnergy();
        energyNext = neb_images[i + 1].getPotentialEnergy();

        // determine the tangent
        if (parameters->nebOldTangent) {
            // old tangent
            neb_tangents[i] = neb_images[i].pbc(posNext - posPrev);
        } else {
            // new improved tangent
            // higherEnergyPrev = energyPrev > energyNext;
            // higherEnergyNext = energyNext > energyPrev;

            if (energyNext > energy && energy > energyPrev) {
                neb_tangents[i] = neb_images[i].pbc(posNext - pos);
            } else if (energy > energyNext && energyPrev > energy) {
                neb_tangents[i] = neb_images[i].pbc(pos - posPrev);
            } else {
                // we are at an extremum
                energyDiffPrev = energyPrev - energy;
                energyDiffNext = energyNext - energy;

                // calculate the energy difference to neighboring nimages
                minDiffEnergy = min(abs(energyDiffPrev), abs(energyDiffNext));
                maxDiffEnergy = max(abs(energyDiffPrev), abs(energyDiffNext));

                // use these energy differences to weight the tangent
                if (energyDiffPrev > energyDiffNext) {
                    neb_tangents[i] = neb_images[i].pbc(posNext - pos) * minDiffEnergy;
                    neb_tangents[i] += neb_images[i].pbc(pos - posPrev) * maxDiffEnergy;
                } else {
                    neb_tangents[i] = neb_images[i].pbc(posNext - pos) * maxDiffEnergy;
                    neb_tangents[i] += neb_images[i].pbc(pos - posPrev) * minDiffEnergy;
                }
            }
        }
        neb_tangents[i].normalize();

        // project the forces and add springs
        force = neb_images[i].getForces();

        // calculate the force perpendicular to the tangent
        forcePerp = force - (force.array() * (neb_tangents[i]).array()).sum() * neb_tangents[i];
        forceSpring = parameters->nebSpring * neb_images[i].pbc((posNext - pos) - (pos - posPrev));

        // calculate the spring force
        distPrev = neb_images[i].pbc(posPrev - pos).squaredNorm();
        distNext = neb_images[i].pbc(posNext - pos).squaredNorm();
        forceSpringPar = parameters->nebSpring * (distNext - distPrev) * neb_tangents[i];

        if (parameters->nebDoublyNudged) {
            forceSpringPerp
                = forceSpring
                  - (forceSpring.array() * (neb_tangents[i]).array()).sum() * neb_tangents[i];
            forceDNEB = forceSpringPerp
                        - (forceSpringPerp.array() * forcePerp.normalized().array()).sum()
                              * forcePerp.normalized();
            if (parameters->nebDoublyNudgedSwitching) {
                long double switching;
                switching
                    = 2.0 / M_PI * atan(pow(forcePerp.norm(), 2) / pow(forceSpringPerp.norm(), 2));
                forceDNEB *= switching;
            }
        } else {
            forceDNEB.setZero();
        }

        if (parameters->nebClimbingImageMethod && i == maxEnergyImage) {
            // we are at the climbing image
            climbingImage = maxEnergyImage;
            projectedForce[i]
                = force
                  - (2.0 * (force.array() * (neb_tangents[i]).array()).sum() * neb_tangents[i])
                  + forceDNEB;
        } else // all non-climbing nimages
        {
            // sum the spring and potential forces for the neb force
            if (parameters->nebElasticBand) {
                projectedForce[i] = forceSpring + force;
            } else {
                projectedForce[i] = forceSpringPar + forcePerp + forceDNEB;
            }
            //*projectedForce[i] = forceSpring + forcePerp;

            // if (parameters->nebFullSpring) {

            movedAfterForceCall = false; // so that we don't repeat a force call
        }

        // zero net translational force
        if (neb_images[i].numberOfFreeAtoms() == neb_images[i].numberOfAtoms()) {
            for (int j = 0; j <= 2; j++) {
                double translationMag = projectedForce[i].col(j).sum();
                int natoms = projectedForce[i].col(j).size();
                projectedForce[i].col(j).array() -= translationMag / ((double) natoms);
            }
        }
    }

    return;
}

// Print NEB image data
void NudgedElasticBand::printImageData(bool writeToFile) {
    double dist, distTotal = 0;
    AtomMatrix tangentStart
        = neb_images[0].pbc(neb_images[1].getPositions() - neb_images[0].getPositions());
    AtomMatrix tangentEnd = neb_images[nimages].pbc(neb_images[nimages + 1].getPositions()
                                                    - neb_images[nimages].getPositions());
    AtomMatrix tang;


    fs::path filePath{"neb.dat"s};
    std::ofstream filedat{filePath};
    if (writeToFile) {
        filedat.open(filePath, std::ofstream::out);
    } else {
         log("Image data (as in neb.dat)\n");
    }

    for (long i = 0; i <= nimages + 1; i++) {
        if (i == 0) {
            tang = tangentStart;
        } else if (i == nimages + 1) {
            tang = tangentEnd;
        } else {
            tang = neb_tangents[i];
        }
        if (i > 0) {
            dist = neb_images[i].distanceTo(neb_images[i - 1]);
            distTotal += dist;
        }
        if (filedat.is_open() && writeToFile) {
            filedat << fmt::sprintf("%3li %12.6f %12.6f %12.6f\n",
                                    i,
                                    distTotal,
                                    neb_images[i].getPotentialEnergy()
                                        - neb_images[0].getPotentialEnergy(),
                                    (neb_images[i].getForces().array() * tang.array()).sum());
        } else {
            log("%3li %12.6f %12.6f %12.6f\n",
                i,
                distTotal,
                neb_images[i].getPotentialEnergy() - neb_images[0].getPotentialEnergy(),
                (neb_images[i].getForces().array() * tang.array()).sum());
        }
    }
    if (writeToFile) {
        filedat.close();
    }
}

// Estimate the barrier using a cubic spline
void NudgedElasticBand::findExtrema(void) {
    // calculate the cubic parameters for each interval (a,b,c,d)

    AtomMatrix tangentEndpoint;
    std::vector<double> a(nimages + 1, 0), b(nimages + 1, 0), c(nimages + 1, 0), d(nimages + 1, 0);
    double F1, F2, U1, U2, dist;

    for (long i = 0; i <= nimages; i++) {
        dist = neb_images[i].distanceTo(neb_images[i + 1]);
        if (i == 0) {
            tangentEndpoint
                = neb_images[i].pbc(neb_images[1].getPositions() - neb_images[0].getPositions());
            tangentEndpoint.normalize();
            F1 = (neb_images[i].getForces().array() * tangentEndpoint.array()).sum() * dist;
        } else {
            F1 = (neb_images[i].getForces().array() * (neb_tangents[i]).array()).sum() * dist;
        }
        if (i == nimages) {
            tangentEndpoint = neb_images[i + 1].pbc(neb_images[nimages + 1].getPositions()
                                                    - neb_images[nimages].getPositions());
            tangentEndpoint.normalize();
            F2 = (neb_images[i + 1].getForces().array() * tangentEndpoint.array()).sum() * dist;
        } else {
            F2 = (neb_images[i + 1].getForces().array() * (neb_tangents[i + 1]).array()).sum()
                 * dist;
        }
        U1 = neb_images[i].getPotentialEnergy();
        U2 = neb_images[i + 1].getPotentialEnergy();
        a[i] = U1;
        b[i] = -F1;
        c[i] = 3. * (U2 - U1) + 2. * F1 + F2;
        d[i] = -2. * (U2 - U1) - (F1 + F2);
    }

    // finding extrema along the MEP

    //    long numExtrema = 0;
    //    double extremaEnergy[2*(nimages+1)]; // the maximum number of extrema
    //    double extremaPosition[2*(nimages+1)];
    double discriminant, f;

    for (long i = 0; i < nimages+1; i++) {
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
                extremumPosition.push_back(i+f);
                extremumEnergy.push_back(d[i] * pow(f, 3) + c[i] * pow(f, 2) + b[i] * f + a[i]);
                extremumCurvature.push_back(6.0 * d[i] * f + 2 * c[i]);
            }
            // cubic case 2
            if (d[i] != 0) {
                f = (-(c[i] - sqrt(discriminant)) / (3. * d[i]));
            }
            if ((f >= 0) && (f <= 1)) {
                extremumPosition.push_back(i + f);
                extremumEnergy.push_back(d[i] * pow(f, 3) + c[i] * pow(f, 2) + b[i] * f + a[i]);
                extremumCurvature.push_back(6 * d[i] * f + 2 * c[i]);
            }
        }
    }

    numExtrema = extremumEnergy.size();
    log("\nFound %li extrema\n", numExtrema);
    log("Energy reference: %f\n", neb_images[0].getPotentialEnergy());
    for (size_t i = 0; i < numExtrema; i++) {
        log("extrema #%li at neb_images position %f with energy %f and curvature %f\n",
            i + 1,
            extremumPosition[i],
            extremumEnergy[i] - neb_images[0].getPotentialEnergy(),
            extremumCurvature[i]);
    }
}

// NEBObjective
VectorXd NEBObjectiveFunction::getGradient(bool fdstep) {
    VectorXd forceV;
    forceV.resize(3 * neb->atoms * neb->nimages);
    if (neb->movedAfterForceCall)
        neb->updateForces();
    for (size_t i = 1; i <= neb->nimages; i++) {
        forceV.segment(3 * neb->atoms * (i - 1), 3 * neb->atoms)
            = VectorXd::Map(neb->projectedForce[i].data(), 3 * neb->atoms);
    }
    return -forceV;
}
double NEBObjectiveFunction::getEnergy() {
    double Energy = 0;
    for (long i = 1; i <= neb->nimages; i++) {
        Energy += neb->neb_images[i].getPotentialEnergy();
    }
    return Energy;
}
void NEBObjectiveFunction::setPositions(VectorXd x) {
    neb->movedAfterForceCall = true;
    for (size_t i = 1; i <= neb->nimages; i++) {
        neb->neb_images[i].setPositions(MatrixXd::Map(
            x.segment(3 * neb->atoms * (i - 1), 3 * neb->atoms).data(), neb->atoms, 3));
    }
}
VectorXd NEBObjectiveFunction::getPositions() {
    VectorXd posV;
    posV.resize(3 * neb->atoms * neb->nimages);
    for (long i = 1; i <= neb->nimages; i++) {
        posV.segment(3 * neb->atoms * (i - 1), 3 * neb->atoms)
            = VectorXd::Map(neb->neb_images[i].getPositions().data(), 3 * neb->atoms);
    }
    return posV;
}
int NEBObjectiveFunction::degreesOfFreedom() { return 3 * neb->nimages * neb->atoms; }
bool NEBObjectiveFunction::isConverged() {
    return getConvergence() < parameters->nebConvergedForce;
}
double NEBObjectiveFunction::getConvergence() { return neb->convergenceForce(); }
VectorXd NEBObjectiveFunction::difference(VectorXd a, VectorXd b) {
    VectorXd pbcDiff(3 * neb->nimages * neb->atoms);
    for (size_t i = 1; i <= neb->nimages; i++) {
        int n = (i - 1) * 3 * neb->atoms;
        int m = 3 * neb->atoms;
        pbcDiff.segment(n, m) = neb->neb_images[i].pbcV(a.segment(n, m) - b.segment(n, m));
    }
    return pbcDiff;
}

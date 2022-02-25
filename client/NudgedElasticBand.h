#ifndef NudgedElasticBand_H
#define NudgedElasticBand_H

#include "HelperFunctions.h"
#include "Log.h"
#include "Matter.h"
#include "Parameters.h"

#include <cmath>
#include <math.h>

class Matter;
class Parameters;

// NEB method for determining a minimum energy path between two matter objects
class NudgedElasticBand {

private:
    Parameters *parameters;

public:
    enum class nebStatus {
        STATUS_GOOD = 0,              // 0
        STATUS_INIT = 1,              // 1
        STATUS_BAD_MAX_ITERATIONS = 2 // 2
    };

    NudgedElasticBand(Matter *initialPassed, Matter *finalPassed, Parameters *parametersPassed)
        : // INITIALIZED by DECLARATION ORDER
          atoms{initialPassed->numberOfAtoms()}, climbingImage{0}, numExtrema{0} {
        parameters = parametersPassed;
        nimages = parameters->nebImages; // readin
        for ([[maybe_unused]] size_t i{0}; i < nimages + 2; i++) {
            neb_images.push_back(*initialPassed);
            neb_tangents.push_back(AtomMatrix::Constant(atoms, 3, 0));
            projectedForce.push_back(AtomMatrix::Constant(atoms, 3, 0));
        }
        log("\nNEB: initialize\n");
        neb_images[nimages + 1] = (*finalPassed);
        AtomMatrix posInitial = neb_images.front().getPositions();
        AtomMatrix posFinal = neb_images.back().getPositions();
        AtomMatrix imgSep = neb_images.front().pbc(posFinal - posInitial) / (nimages + 1);
        for (long double idx{0}; auto &&img : neb_images) {
            if ((idx == 0) /*initial*/ || (idx == nimages + 1) /*final*/) {
                continue;
            }
            img.setPositions(posInitial + imgSep * idx);
            ++idx;
        }

        movedAfterForceCall = true;

        // Make sure that the endpoints know their energy
        neb_images.front().getPotentialEnergy();
        neb_images.back().getPotentialEnergy();
    }

    ~NudgedElasticBand(){};

    int compute(void);
    void updateForces(void);
    double convergenceForce(void);
    void findExtrema(void);
    void printImageData(bool writeToFile = false);

    size_t atoms, nimages, climbingImage, numExtrema;
    std::vector<Matter> neb_images; // NEB images
    std::vector<AtomMatrix> neb_tangents;
    std::vector<AtomMatrix> projectedForce;
    bool movedAfterForceCall;
    std::vector<double> extremumEnergy;
    std::vector<double> extremumPosition;
    std::vector<double> extremumCurvature;

    long maxEnergyImage;
};

class NEBObjectiveFunction : public ObjectiveFunction {
public:
    NEBObjectiveFunction(NudgedElasticBand *nebPassed, Parameters *parametersPassed)
        : neb{nebPassed}, parameters{parametersPassed} {}
    ~NEBObjectiveFunction(){};
    double getEnergy();
    double getConvergence();
    bool isConverged();
    int degreesOfFreedom();
    VectorXd getPositions();
    VectorXd getGradient(bool fdstep = false);
    VectorXd difference(VectorXd a, VectorXd b);
    void setPositions(VectorXd x);

private:
    NudgedElasticBand *neb;
    Parameters *parameters;
};
#endif

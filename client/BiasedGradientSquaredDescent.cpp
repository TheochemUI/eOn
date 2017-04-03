
#include "Log.h"
#include "Matter.h"
#include "HelperFunctions.h"
#include "Optimizer.h"
#include "ObjectiveFunction.h"
#include "BiasedGradientSquaredDescent.h"
#include "Lanczos.h"
#include "Dimer.h"
#include "ImprovedDimer.h"

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <string.h>
#include <map>

//static const char LOG_PREFIX[] = "[BGSD]";

class BGSDObjectiveFunction : public ObjectiveFunction
{
    public:
        BGSDObjectiveFunction(Matter *matterPassed,
                              double reactantEnergyPassed,
			      double bgsdAlphaPassed,
                              Parameters *parametersPassed)
        {
            matter = matterPassed;
            parameters = parametersPassed;
	    bgsdAlpha = bgsdAlphaPassed;
            reactantEnergy = reactantEnergyPassed;
        }

        ~BGSDObjectiveFunction(void){};

        double getEnergy()
        {
            VectorXd Vforce = matter->getForcesFreeV();
            double Henergy = 0.5 * Vforce.dot(Vforce) + 0.5 * bgsdAlpha * (matter->getPotentialEnergy() - (reactantEnergy + parameters->beta)) * (matter->getPotentialEnergy() -(reactantEnergy + parameters->beta));
            return Henergy;
        }

        VectorXd getGradient(bool fdstep=false) 
        {   
            VectorXd Vforce = matter->getForcesFreeV();
            double magVforce = Vforce.norm();
            VectorXd normVforce = Vforce/magVforce;
            VectorXd Vpositions = matter->getPositionsFreeV();
            matter -> setPositionsFreeV(matter->getPositionsFreeV() - normVforce * parameters->gradientfinitedifference);
            VectorXd Vforcenew = matter->getForcesFreeV();
            matter -> setPositionsFreeV(Vpositions);
            VectorXd Hforce = magVforce * (Vforcenew - Vforce)/parameters->gradientfinitedifference + bgsdAlpha * (matter->getPotentialEnergy() - (reactantEnergy + parameters->beta)) * Vforce;
            return -Hforce; }

        double getGradientnorm() { VectorXd Hforce = getGradient();
                                   double Hnorm = Hforce.norm();
                                    return Hnorm;}

        void setPositions(VectorXd x) { matter->setPositionsFreeV(x); }
        VectorXd getPositions() { return matter->getPositionsFreeV(); }
        int degreesOfFreedom() { return 3*matter->numberOfFreeAtoms(); }
        bool isConverged() { return isConvergedH() && isConvergedV() ; }
        bool isConvergedH() { return getConvergenceH() < parameters->Hforceconvergence; }
        bool isConvergedV() { return getConvergenceV() < parameters->grad2energyconvergence;} 
        bool isConvergedIP() { return getConvergenceH() < parameters->grad2forceconvergence; }

        double getConvergence() { return getEnergy() && getGradient().norm(); }
        double getConvergenceH() { return getGradient().norm(); }
        double getConvergenceV() { return getEnergy(); }
        VectorXd difference(VectorXd a, VectorXd b) {
            return matter->pbcV(a-b);
        }
    private:
        Matter *matter;
        Parameters *parameters;
        double reactantEnergy;
	double bgsdAlpha;  
};



BiasedGradientSquaredDescent::BiasedGradientSquaredDescent(Matter *matterPassed, double reactantEnergyPassed, Parameters *parametersPassed)
{
    parameters = parametersPassed;
    reactantEnergy = reactantEnergyPassed;
    saddle = matterPassed;
    eigenvector.resize(saddle->numberOfAtoms(), 3);
    eigenvector.setZero();
}

BiasedGradientSquaredDescent::~BiasedGradientSquaredDescent()
{
}

int BiasedGradientSquaredDescent::run()
{
    BGSDObjectiveFunction objf(saddle, reactantEnergy, parameters->alpha, parameters);
    Optimizer *optimizer = Optimizer::getOptimizer(&objf, parameters);
    int iteration = 0;
    printf("starting optimization of H with parameters alpha and beta: %.2f %.2f\n",parameters->alpha,parameters->beta);
    while (!objf.isConvergedH() || iteration == 0) {
        optimizer->step(parameters->optMaxMove);
        printf("iteration %i Henergy, gradientHnorm, and Venergy: %.8f %.8f  %.8f\n", iteration, objf.getEnergy(),objf.getGradientnorm(),saddle->getPotentialEnergy());           
        iteration++;
    }
    BGSDObjectiveFunction objf2(saddle, reactantEnergy, 0.0, parameters);
    Optimizer *optimizer2 = Optimizer::getOptimizer(&objf2, parameters);
    while (!objf2.isConvergedV() || iteration == 0) {
        if (objf2.isConvergedIP()) {break;};
        optimizer2->step(parameters->optMaxMove);
        printf("gradient squared iteration %i Henergy, gradientHnorm, and Venergy: %.8f %.8f  %.8f\n", iteration, objf2.getEnergy(),objf2.getGradientnorm(),saddle->getPotentialEnergy());
        iteration++;
    }

    LowestEigenmode *minModeMethod;
    if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
        if (parameters->dimerImproved) {
            minModeMethod = new ImprovedDimer(saddle, parameters);
         }else{
            minModeMethod = new Dimer(saddle, parameters);
         }
    }else if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
        minModeMethod = new Lanczos(saddle, parameters);
    }

 //   eigenvector.setZero();
 //   eigenvector(384,0) = 1.;
    eigenvector.setRandom();
    for (int i=0;i<saddle->numberOfAtoms();i++){
        for (int j=0;j<3;j++) {
            if (saddle->getFixed(i)){
                eigenvector(i,j) = 0.0;
            };
        }
    }
    eigenvector.normalize();
    minModeMethod->compute(saddle, eigenvector);
    eigenvector = minModeMethod->getEigenvector();
    eigenvalue = minModeMethod->getEigenvalue();
    printf("lowest eigenvalue %.8f\n",eigenvalue);
    if (objf2.isConvergedV()) {return 0;}
    else if (objf2.isConvergedIP()) {return 1;}
    else {return 1;};
}

double BiasedGradientSquaredDescent::getEigenvalue()
{
    return eigenvalue;
}

AtomMatrix BiasedGradientSquaredDescent::getEigenvector()
{
    return eigenvector;
}

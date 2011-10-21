#ifndef OBJECTIVE_FUNCTION_H
#define OBJECTIVE_FUNCTION_H

#include "Matter.h"
#include "LowestEigenmode.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "ExactMinMode.h"

class ObjectiveFunction
{
    public:
        virtual ~ObjectiveFunction() {}
        virtual double getEnergy()=0;
        virtual VectorXd getGradient()=0;
        virtual void setPositions(VectorXd x)=0;
        virtual VectorXd getPositions()=0;
        virtual int degreesOfFreedom()=0;
        virtual bool isConverged()=0;
        virtual double getConvergence()=0;
};

class MatterObjectiveFunction : public ObjectiveFunction
{
    public:
        MatterObjectiveFunction(Matter *matterPassed,
                                Parameters *parametersPassed)
        {
            matter = matterPassed;
            parameters = parametersPassed;
        }
        ~MatterObjectiveFunction(void){};
        double getEnergy() { return matter->getPotentialEnergy(); }
        VectorXd getGradient() { return -matter->getForcesFreeV(); }
        void setPositions(VectorXd x) { matter->setPositionsFreeV(x); }
        VectorXd getPositions() { return matter->getPositionsFreeV(); }
        int degreesOfFreedom() { return 3*matter->numberOfFreeAtoms(); }
        bool isConverged() { return getConvergence() < parameters->optConvergedForce; }
        double getConvergence() { return matter->maxForce(); }
    private:
        Matter *matter;
        Parameters *parameters;
};

class MinModeObjectiveFunction : public ObjectiveFunction
{
    public:
        MinModeObjectiveFunction(Matter *matterPassed, AtomMatrix modePassed, 
                                 Parameters *parametersPassed)
        {
            matter = matterPassed;
            eigenvector = modePassed;
            parameters = parametersPassed;
            if(parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER)
            {
                if(parameters->dimerImproved)
                {
                    minModeMethod = new ImprovedDimer(matter, parameters);
                }
                else
                {
                    minModeMethod = new Dimer(matter, parameters);
                }
                
            }
            else if(parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS)
            {
                minModeMethod = new Lanczos(matter, parameters);
            }
            else if(parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_EXACT)
            {
                minModeMethod = new ExactMinMode(matter, parameters);
            }
        }
        ~MinModeObjectiveFunction(void){};
        double getEnergy() { return matter->getPotentialEnergy(); }
        VectorXd getGradient() 
        { 
            AtomMatrix proj;
            AtomMatrix force = matter->getForces();

            minModeMethod->compute(matter, eigenvector);
            eigenvector = minModeMethod->getEigenvector();
            eigenvalue = minModeMethod->getEigenvalue();

            proj = (force.cwise() * eigenvector).sum() * eigenvector.normalized();

            if (0 < eigenvalue) {
                if (parameters->saddlePerpForceRatio > 0.0) 
                {
                    // reverse force parallel to eigenvector, and reduce perpendicular force
                    double const d = parameters->saddlePerpForceRatio;
                    force = d*force - (1.+d)*proj;
                }
                else 
                {
                    // follow eigenmode
                    force = -proj;
                }
            }else{
                // reversing force parallel to eigenmode
                force += -2.*proj;
            }

            VectorXd forceV = VectorXd::Map(force.data(), 3*matter->numberOfAtoms());
            return -forceV;
        }
        void setPositions(VectorXd x) { matter->setPositionsV(x); }
        VectorXd getPositions() { return matter->getPositionsV(); }
        int degreesOfFreedom() { return 3*matter->numberOfAtoms(); }
        bool isConverged() { return getConvergence() < parameters->optConvergedForce; }
        double getConvergence() { return matter->maxForce(); }
        AtomMatrix eigenvector;
        double eigenvalue;
        LowestEigenmode *minModeMethod;

    private:
        Matter *matter;
        Parameters *parameters;

};

#endif

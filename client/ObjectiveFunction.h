#ifndef OBJECTIVE_FUNCTION_H
#define OBJECTIVE_FUNCTION_H

#include "Matter.h"

class ObjectiveFunction
{
    public:
        virtual ~ObjectiveFunction() {}
        virtual double getEnergy()=0;
        virtual VectorXd getGradient()=0;
        virtual void setPositions(VectorXd x)=0;
        virtual VectorXd getPositions()=0;
        virtual int degreesOfFreedom()=0;
        virtual bool converged()=0;
        virtual double convergence()=0;
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
        bool converged() { return convergence() < parameters->optConvergedForce; }

        double convergence()
        {
            int freeAtoms = matter->numberOfFreeAtoms();
            double maxAtomForce=0;
            AtomMatrix force = matter->getForcesFree();
            for (int i=0; i<freeAtoms; i++) {
                maxAtomForce = max(maxAtomForce, force.row(i).norm());
            }

            return maxAtomForce;
        }
    private:
        Matter *matter;
        Parameters *parameters;
};

#endif

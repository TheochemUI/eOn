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

#endif

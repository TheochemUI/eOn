#ifndef _BIASED_GRADIENT_SQUARED_DESCENT_
#define _BIASED_GRADIENT_SQUARED_DESCENT_

#include "Eigen.h"
#include "Matter.h"
#include "SaddleSearchMethod.h"
#include "MinModeSaddleSearch.h"
#include <vector>

class BiasedGradientSquaredDescent : public SaddleSearchMethod
{
    public:
        BiasedGradientSquaredDescent(Matter *matterPassed, double reactantEnergyPassed, Parameters *parametersPassed);
        ~BiasedGradientSquaredDescent();

        int run(void);
        double getEigenvalue();
        AtomMatrix getEigenvector();

        double eigenvalue;
        AtomMatrix eigenvector;

        Matter *saddle;

        int status;

    private:
        Parameters *parameters;
        double reactantEnergy;
//        double bgsdAlpha;
};

#endif

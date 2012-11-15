#ifndef DYNAMICSSADDLESEARCH_H
#define DYNAMICSSADDLESEARCH_H

#include "Eigen.h"
#include "Matter.h"
#include "SaddleSearchMethod.h"
#include "MinModeSaddleSearch.h"
#include <vector>

class DynamicsSaddleSearch : public SaddleSearchMethod
{
    public:
        DynamicsSaddleSearch(Matter *matterPassed, Parameters *parametersPassed);
        ~DynamicsSaddleSearch();

        int run(void);
        double getEigenvalue();
        AtomMatrix getEigenvector();

        int refineTransition(std::vector<Matter*>, Matter *product);

        double eigenvalue;
        AtomMatrix eigenvector;

        double time;

        Matter *reactant;
        Matter *saddle;
        Matter *product;

        int status;

    private:
        Parameters *parameters;
};

#endif

#ifndef DYNAMICSSADDLESEARCH_H
#define DYNAMICSSADDLESEARCH_H

#include "Eigen.h"
#include "Matter.h"
#include "SaddleSearchMethod.h"
#include "MinModeSaddleSearch.h"

class DynamicsSaddleSearch : public SaddleSearchMethod
{
    public:
        DynamicsSaddleSearch(Matter *matterPassed, Parameters *parametersPassed);
        ~DynamicsSaddleSearch();

        int run(void);
        double getEigenvalue();
        AtomMatrix getEigenvector();

        double eigenvalue;
        AtomMatrix eigenvector;

        Matter *reactant;
        Matter *saddle;
        Matter *product;

        int status;

    private:
        Parameters *parameters;
};

#endif

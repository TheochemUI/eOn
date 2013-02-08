#ifndef BASINHOPPINGSADDLESEARCH_H
#define BASINHOPPINGSADDLESEARCH_H

#include "Eigen.h"
#include "Matter.h"
#include "SaddleSearchMethod.h"
#include "MinModeSaddleSearch.h"
#include <vector>

class BasinHoppingSaddleSearch : public SaddleSearchMethod
{
    public:
        BasinHoppingSaddleSearch(Matter *reactant, Matter *displacement, Parameters *parametersPassed);
        ~BasinHoppingSaddleSearch();

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

#ifndef SADDLESEARCHMETHOD_H
#define SADDLESEARCHMETHOD_H

#include "Eigen.h"

class SaddleSearchMethod {
    public:

        virtual ~SaddleSearchMethod() {};
        virtual int run()=0;
        virtual double getEigenvalue()=0;
        virtual AtomMatrix getEigenvector()=0;

        int status;
};

#endif

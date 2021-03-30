#ifndef LANCZOS_H
#define LANCZOS_H

#include "Eigen.h"
#include "Matter.h"
#include "Parameters.h"
#include "LowestEigenmode.h"

// Lanczos method to find the lowest curvature mode
class Lanczos : public LowestEigenmode
{

    public:
        Lanczos(Matter *matter, Parameters *parameters);
        ~Lanczos();

        void compute(Matter *matter, AtomMatrix direction);
        double getEigenvalue();
        AtomMatrix getEigenvector();

    private:
        Parameters *parameters;
        AtomMatrix lowestEv;
        double lowestEw;
};

#endif



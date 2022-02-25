#ifndef PY_POTENTIAL_H
#define PY_POTENTIAL_H

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"
#include "py_wrapper.hpp"

// Potentials
template <class PotentialBase = Potential>
class PyPotential : public PotentialBase {
public:
    /* Inherit the constructors */
    using PotentialBase::PotentialBase;
    void initialize() override {
        PYBIND11_OVERRIDE_PURE(void,          /* Return type */
                               PotentialBase, /* Parent class */
                               initialize /* Name of function in C++ (must match Python name) */
                               , );       /* Trailing comma for no arguments */
    };
    void force(long nAtoms,
               const double *positions,
               const int *atomicNrs,
               double *forces,
               double *energy,
               const double *box,
               int nImages) override {
        PYBIND11_OVERRIDE_PURE(void,      /* Return type */
                               Potential, /* Parent class */
                               force,     /* Name of function in C++ (must match Python name) */
                               nAtoms,    /* Argument(s) */
                               positions,
                               atomicNrs,
                               forces,
                               energy,
                               box,
                               nImages);
    };
};

#endif /* PY_POTENTIAL_H */

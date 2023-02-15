#ifndef PY_POTENTIAL_H_
#define PY_POTENTIAL_H_

#include <pybind11/pybind11.h>

#include "../../Potential.h"
#include "py_wrapper.hpp"

template <class PotentialBase = Potential>
class PyPotential : public PotentialBase {
public:
    /* Constructors and inherited */
    using PotentialBase::PotentialBase;
    void force(long nAtoms,
               const double *positions,
               const int *atomicNrs,
               double *forces,
               double *energy,
               const double *box) override {
        PYBIND11_OVERRIDE_PURE(void,
                               PotentialBase,
                               force,
                               nAtoms,
                               positions,
                               atomicNrs,
                               forces,
                               energy,
                               box);
    }
};

#endif // PY_POTENTIAL_H_

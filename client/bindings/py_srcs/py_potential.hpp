#ifndef PY_POTENTIAL_H_
#define PY_POTENTIAL_H_

#include <pybind11/pybind11.h>

#include "../../Potential.h"
#include "py_wrapper.hpp"
#include <utility>

template <class PotentialBase = Potential>
class PyPotential : public PotentialBase {
public:
    /* Constructors and inherited */
    using PotentialBase::PotentialBase;
    std::pair<double, AtomMatrix> get_ef(const AtomMatrix positions,
               const VectorXi atomicNrs,
               const Matrix3d box) override {
        // This is needed to prevent substitution errors
        using Return = std::pair<double, AtomMatrix>;
        PYBIND11_OVERRIDE_PURE(Return,
                               PotentialBase,
                               get_ef,
                               positions,
                               atomicNrs,
                               box);
    }
};

#endif // PY_POTENTIAL_H_

#ifndef PY_WRAPPER_H
#define PY_WRAPPER_H

#include <fstream>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../Parameters.h"
#include "../Matter.h"
#include "../Potential.h"
#include "../potentials/Morse/Morse.h"

using namespace std::string_literals; // For ""s
using namespace pybind11::literals;   // For ""_a
namespace py = pybind11;              // Convention

void py_parameters(py::module_ &m);
void py_matter(py::module_ &m);
void py_potential(py::module_ &m);
void py_morse(py::module_ &m);

template <class PotentialBase = Potential>
class PyPotential : public PotentialBase {
public:
    /* Inherit the constructors */
    using PotentialBase::PotentialBase;
    void initialize() override{
        PYBIND11_OVERRIDE_PURE(void,          /* Return type */
                               PotentialBase, /* Parent class */
                               initialize /* Name of function in C++ (must match Python name) */
                               ,); /* Trailing comma for no arguments */
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
template <class MorseBase = Morse> class PyMorse : public PyPotential<MorseBase>{
        public:
        using PyPotential<MorseBase>::PyPotential; // Inherit constructor
        // Override pure virtual with non-pure
        void initialize() override { PYBIND11_OVERRIDE(void, MorseBase, initialize,); }
    void force(long nAtoms,
               const double *positions,
               const int *atomicNrs,
               double *forces,
               double *energy,
               const double *box,
               int nImages) override {
        PYBIND11_OVERRIDE(void,      /* Return type */
                               MorseBase, /* Parent class */
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

#endif /* PY_WRAPPER_H */

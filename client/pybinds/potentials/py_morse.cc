// clang-format off
#include "../py_wrapper.hpp"
#include "../py_potential.hpp"
// Binding
#include "../../potentials/Morse/Morse.h"
// clang-format on

template <class MorseBase = Morse>
class PyMorse : public PyPotential<MorseBase> {
public:
    using PyPotential<MorseBase>::PyPotential; // Inherit constructor
    // Override pure virtual with non-pure
    void initialize() override { PYBIND11_OVERRIDE(void, MorseBase, initialize, ); }
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


void py_morse(py::module_ &m) {
    py::class_<Morse, Potential, PyMorse<> >(m, "Morse")
        /*
        ** Constructors
        */
        .def(py::init<>())
        .def(py::init<double /*re*/, double /*De*/, double /*a*/, double /*cutoff*/>())

        /*
        ** Operators
        */

        /*
        ** Methods
        */
        .def("cleanMemory", &Morse::cleanMemory)
        .def("energy_and_forces", &Morse::energy_and_forces,
             py::arg("nAtoms"), "positions"_a, "box"_a)
        .def("force",
             py::overload_cast<long /*N*/,
                               const double* /*R*/,
                               const int* /*atomicNrs*/,
                               double* /*F*/,
                               double* /*U*/,
                               const double* /*box*/,
                               int /*nImages*/
                               >(&Morse::force),
             py::arg("N"), "R"_a, "atomicNrs"_a, "F"_a, "U"_a, "box"_a, "nImages"_a)
        .def("setParameters", &Morse::setParameters)
        .def("initialze", &Morse::initialize)
        /*
        ** Parameters
        */

        /*
        ** Python helpers
        */

        .def("__repr__", [](const Morse &a) { return "<Morse object>"; });
}

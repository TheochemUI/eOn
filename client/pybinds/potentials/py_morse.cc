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
        /** TODO: This should also be part of Potential, no need to make separate bindings **/
        .def("ef_matter", [](Morse &mpot, Matter mat){
            Parameters params{mat.getParameters()};
            params.potential = "morse_pt";
            mpot.setParams(&params);
            mat.setPotential(&mpot);
            return std::make_pair(mat.getPotentialEnergy(), mat.getForces());
        }, py::arg("matter"))
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

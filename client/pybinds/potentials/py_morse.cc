// clang-format off
#include "../py_wrapper.hpp"
// Additional
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
// clang-format on

void py_morse(py::module_ &m) {
    py::class_<Morse, Potential>(m, "Morse")
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

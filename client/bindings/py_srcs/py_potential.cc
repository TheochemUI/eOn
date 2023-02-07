#include "py_potential.hpp"

void py_potential(py::module &m) {
    py::class_<Potential, PyPotential<>>(m, "Potential")
        .def(py::init())
        /* Functions */
        .def_static("getPotential", &Potential::getPotential)
        .def("initialize", &Potential::initialize)
        /* Book-keeping */
        /* TODO: This is static, so wrong when multiple potentials are used? */
        .def_readwrite_static("fcalls", &Potential::fcalls, "Number of force calls")
        .def_readwrite_static("fcallsTotal", &Potential::fcallsTotal, "Total force calls")
        .def_readwrite_static("totalUserTime", &Potential::fcalls, "Time taken")
        /* Objects */
        /* TODO: This is likely to be broken, params is a pointer */
        .def_readwrite("params", &Potential::params, "Parameters")
        .def_readwrite_static("pot", &Potential::pot, "Potential object")
        /* Helpers */
        .def("__repr__", [](const Potential &pot) { return pot.getName(); });
};

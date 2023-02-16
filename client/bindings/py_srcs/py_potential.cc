#include "py_potential.hpp"

void py_potential(py::module &m) {
    py::class_<Potential, PyPotential<>>(m, "Potential")
        .def(py::init<Parameters*>())
        .def("get_ef", &Potential::get_ef)
        .def("getType", &Potential::getType)
        /* Helpers */
        .def("__repr__", [](Potential &pot) { return helper_functions::getPotentialName( pot.getType() ); });
    m.def("makePotential", &helper_functions::makePotential);
};

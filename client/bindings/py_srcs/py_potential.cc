#include "py_potential.hpp"

void py_potential(py::module &m) {
    py::class_<Potential, PyPotential<>>(m, "Potential")
        .def(py::init<Parameters*>())
        .def("get_ef", &Potential::get_ef)
        .def("getType", &Potential::getType)
        .def(
            "callPot",
            [](Potential *pot, Matter &mat) {
                // Store and restore
                Potential* oldPot = mat.getPotential();
                mat.setPotential(pot);
                double e_pot{mat.getPotentialEnergy()};
                AtomMatrix f_pot {mat.getForces()};
                mat.setPotential(oldPot);
                return std::make_pair(e_pot, f_pot);
            },
            py::arg("matter"))
        /* Helpers */
        .def("__repr__", [](Potential &pot) { return helper_functions::getPotentialName( pot.getType() ); });
    m.def("makePotential", &helper_functions::makePotential);
};

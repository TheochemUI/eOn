#include "py_potential.hpp"

void py_potential(py::module_ &m) {
    py::class_<Potential, PyPotential<>>(m, "Potential")
        .def(py::init())
        /*
        ** Functions
        */
        .def_static("getPotential", &Potential::getPotential)
        .def("initialize", &Potential::initialize)
        .def(
            "force",
            [](Potential &pot,
               long natms,
               AtomMatrix pos,
               VectorXi atmnrs,
               Matrix3d box,
               int nimgs) {
                AtomMatrix forces{AtomMatrix::Zero(natms, 3)};
                double energy = 0;
                pot.force(
                    natms, pos.data(), atmnrs.data(), forces.data(), &energy, box.data(), nimgs);
                return std::make_pair(energy, forces);
            },
            "Get forces",
            py::arg("nAtoms"),
            "positions"_a,
            "atomicNrs"_a,
            "box"_a,
            "nImages"_a)

        /*
        ** Parameters
        */
        // Potential Names

        .def_property_readonly_static(
            "POT_EMT",
            [](py::object /*self*/) { return PotentialStrings::POT_EMT; },
            "EMT potential string")
        .def_property_readonly_static(
            "POT_EXT",
            [](py::object /*self*/) { return PotentialStrings::POT_EXT; },
            "EXT potential string")
        .def_property_readonly_static(
            "POT_LJ",
            [](py::object /*self*/) { return PotentialStrings::POT_LJ; },
            "LJ potential string")
        .def_property_readonly_static(
            "POT_LJCLUSTER",
            [](py::object /*self*/) { return PotentialStrings::POT_LJCLUSTER; },
            "LJ cluster potential string")
        .def_property_readonly_static(
            "POT_MORSE_PT",
            [](py::object /*self*/) { return PotentialStrings::POT_MORSE_PT; },
            "Morse potential string")
        .def_property_readonly_static(
            "POT_NEW",
            [](py::object /*self*/) { return PotentialStrings::POT_NEW; },
            "New potential string")
#ifdef IMD_POT
        .def_property_readonly_static(
            "POT_IMD",
            [](py::object /*self*/) { return PotentialStrings::POT_IMD; },
            "IMD potential string")
#endif
#ifdef GPR_POT
        .def_property_readonly_static(
            "POT_GPR",
            [](py::object /*self*/) { return PotentialStrings::POT_GPR; },
            "GPR potential string")
#endif
#ifdef WITH_WATER
        .def_property_readonly_static(
            "POT_TIP4P",
            [](py::object /*self*/) { return PotentialStrings::POT_TIP4P; },
            "TIP4P potential string")
        .def_property_readonly_static(
            "POT_TIP4P_PT",
            [](py::object /*self*/) { return PotentialStrings::POT_TIP4P_PT; },
            "TIP4P_PT potential string")
#    ifdef WITH_FORTRAN
        .def_property_readonly_static(
            "POT_TIP4P_H",
            [](py::object /*self*/) { return PotentialStrings::POT_TIP4P_H; },
            "TIP4P_H potential string")
#    endif
        .def_property_readonly_static(
            "POT_SPCE",
            [](py::object /*self*/) { return PotentialStrings::POT_SPCE; },
            "SPCE potential string")
#endif
#ifdef WITH_FORTRAN
        .def_property_readonly_static(
            "POT_EAM_AL",
            [](py::object /*self*/) { return PotentialStrings::POT_EAM_AL; },
            "EAM Aluminum potential string")
        .def_property_readonly_static(
            "POT_EDIP",
            [](py::object /*self*/) { return PotentialStrings::POT_EDIP; },
            "EDIP potential string")
        .def_property_readonly_static(
            "POT_FEHE",
            [](py::object /*self*/) { return PotentialStrings::POT_FEHE; },
            "FEHE potential string")
        .def_property_readonly_static(
            "POT_LENOSKY_SI",
            [](py::object /*self*/) { return PotentialStrings::POT_LENOSKY_SI; },
            "LENOSKY SI potential string")
        .def_property_readonly_static(
            "POT_SW_SI",
            [](py::object /*self*/) { return PotentialStrings::POT_SW_SI; },
            "SW SI potential string")
        .def_property_readonly_static(
            "POT_TERSOFF_SI",
            [](py::object /*self*/) { return PotentialStrings::POT_TERSOFF_SI; },
            "TERSOFF SI potential string")
#endif
#ifndef WIN32
#    ifdef WITH_VASP
        .def_property_readonly_static(
            "POT_VASP",
            [](py::object /*self*/) { return PotentialStrings::POT_VASP; },
            "VASP potential string")
#    endif
#endif
#ifdef LAMMPS_POT
        .def_property_readonly_static(
            "POT_LAMMPS",
            [](py::object /*self*/) { return PotentialStrings::POT_LAMMPS; },
            "LAMMPS potential string")
#endif
#ifdef EONMPI
        .def_property_readonly_static(
            "POT_MPI",
            [](py::object /*self*/) { return PotentialStrings::POT_MPI; },
            "MPI potential string")
#endif
#ifdef WITH_PYTHON
#    ifdef PYAMFF_POT
        .def_property_readonly_static(
            "POT_PYAMFF",
            [](py::object /*self*/) { return PotentialStrings::POT_PYAMFF; },
            "PYAMFF potential string")
#    endif
        .def_property_readonly_static(
            "POT_QSC",
            [](py::object /*self*/) { return PotentialStrings::POT_QSC; },
            "QSC potential string")
#endif

#ifdef WITH_AMS
        .def_property_readonly_static(
            "POT_AMS",
            [](py::object /*self*/) { return PotentialStrings::POT_AMS; },
            "AMS potential string")
        .def_property_readonly_static(
            "POT_AMS_IO",
            [](py::object /*self*/) { return PotentialStrings::POT_AMS_IO; },
            "AMS_IO potential string")
#endif
        // Book-keeping
        .def_readwrite_static("fcalls", &Potential::fcalls, "Number of force calls")
        .def_readwrite_static("fcallsTotal", &Potential::fcallsTotal, "Total force calls")
        .def_readwrite_static("totalUserTime", &Potential::totalUserTime, "Time taken")

        // Objects
        .def_readwrite("params", &Potential::params, "Parameters")
        .def_readwrite_static("pot", &Potential::pot, "Potential object")

        /*
        ** Python helpers
        */

        .def("__repr__", [](const Potential &a) { return "<Potential object>"; });
}

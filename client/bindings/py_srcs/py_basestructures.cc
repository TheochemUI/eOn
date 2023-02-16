// clang-format off
#include "py_wrapper.hpp"
// clang-format on

void py_basestructures(py::module_ &m) {
    py::enum_<PotType>(m, "PotType")
        .value("EMT", PotType::EMT)
        .value("EXT", PotType::EXT)
        .value("LJ", PotType::LJ)
        .value("LJCLUSTER", PotType::LJCLUSTER)
        .value("MORSE_PT", PotType::MORSE_PT)
        .value("NEW", PotType::NEW)
        .value("CUH2", PotType::CUH2)
        .value("IMD", PotType::IMD)
        .value("TIP4P", PotType::TIP4P)
        .value("TIP4P_PT", PotType::TIP4P_PT)
        .value("TIP4P_H", PotType::TIP4P_H)
        .value("SPCE", PotType::SPCE)
        .value("EAM_AL", PotType::EAM_AL)
        .value("EDIP", PotType::EDIP)
        .value("FEHE", PotType::FEHE)
        .value("LENOSKY_SI", PotType::LENOSKY_SI)
        .value("SW_SI", PotType::SW_SI)
        .value("TERSOFF_SI", PotType::TERSOFF_SI)
        .value("VASP", PotType::VASP)
        .value("LAMMPS", PotType::LAMMPS)
        .value("MPI", PotType::MPI)
        .value("PYAMFF", PotType::PYAMFF)
        .value("QSC", PotType::QSC)
        .value("UNKNOWN", PotType::UNKNOWN)
        .value("AMS", PotType::AMS)
        .value("AMS_IO", PotType::AMS)
        .value("GPR", PotType::GPR)
        .value("PYTHON", PotType::PYTHON);

    m.def("getPotentialType", &helper_functions::getPotentialType);
    m.def("getPotentialName", &helper_functions::getPotentialName);
}

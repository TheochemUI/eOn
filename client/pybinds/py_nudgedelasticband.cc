// clang-format off
#include "py_wrapper.hpp"
#include "../NudgedElasticBand.h"
// clang-format on

void py_nudgedelasticband(py::module_ &m) {
    py::class_<NudgedElasticBand>(m, "NudgedElasticBand")
        /*
        ** Constructors
        */
        .def(py::init<Matter * /*initialPassed*/,
                      Matter * /*finalPassed*/,
                      Parameters * /*parametersPassed*/>())

        /*
        ** Methods
        */
        .def("compute", &NudgedElasticBand::compute)
        .def("updateForces", &NudgedElasticBand::updateForces)
        .def("convergenceForce", &NudgedElasticBand::convergenceForce)
        .def("findExtrema", &NudgedElasticBand::findExtrema)
        .def("printImageData", &NudgedElasticBand::printImageData)

        /*
        ** Parameters
        */
        .def_readwrite("atoms", &NudgedElasticBand::atoms)
        .def_readwrite("nimages", &NudgedElasticBand::nimages)
        .def_readwrite("climbingImage", &NudgedElasticBand::climbingImage)
        .def_readwrite("numExtrema", &NudgedElasticBand::numExtrema)
        .def_readwrite("neb_images", &NudgedElasticBand::neb_images)
        .def_readwrite("neb_tangents", &NudgedElasticBand::neb_tangents)
        .def_readwrite("projectedForce", &NudgedElasticBand::projectedForce)
        .def_readwrite("movedAfterForceCall", &NudgedElasticBand::movedAfterForceCall)
        .def_readwrite("extremumEnergy", &NudgedElasticBand::extremumEnergy)
        .def_readwrite("extremumPosition", &NudgedElasticBand::extremumPosition)
        .def_readwrite("extremumCurvature", &NudgedElasticBand::extremumCurvature)
        .def_readwrite("maxEnergyImage", &NudgedElasticBand::maxEnergyImage)

        /*
        ** Python helpers
        */

        .def("__repr__", [](const NudgedElasticBand &a) { return "<NudgedElasticBand object>"; });

    /* Internal Types and Enums */
    py::enum_<NudgedElasticBand::nebStatus>(m, "nebStatus")
        .value("STATUS_GOOD", NudgedElasticBand::nebStatus::STATUS_GOOD)
        .value("STATUS_INIT", NudgedElasticBand::nebStatus::STATUS_INIT)
        .value("STATUS_BAD_MAX_ITERATIONS", NudgedElasticBand::nebStatus::STATUS_BAD_MAX_ITERATIONS);
        // .export_values(); // exports entries into the parent which makes no sense
}

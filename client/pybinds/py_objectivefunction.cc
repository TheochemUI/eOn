#include "py_objectivefunction.hpp"

void py_objectivefunction(py::module_ &m) {
    py::class_<ObjectiveFunction, PyObjectiveFunction<> >(m, "ObjectiveFunction")
        .def(py::init())
        /*
        ** Functions
        */
        .def("getEnergy", &ObjectiveFunction::getEnergy)
        .def("getConvergence", &ObjectiveFunction::getConvergence)
        .def("isConverged", &ObjectiveFunction::isConverged)
        .def("degreesOfFreedom", &ObjectiveFunction::degreesOfFreedom)
        .def("getPositions", &ObjectiveFunction::getPositions)
        .def("getGradient", &ObjectiveFunction::getGradient, py::arg("fdstep")=false)
        .def("difference", &ObjectiveFunction::difference, py::arg("a"),"b"_a)
        .def("setPositions", &ObjectiveFunction::setPositions)

        /*
        ** Python helpers
        */

        .def("__repr__", [](const ObjectiveFunction &a) { return "<ObjectiveFunction object>"; });
}

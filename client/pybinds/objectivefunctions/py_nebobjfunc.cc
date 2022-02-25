// clang-format off
#include "../py_wrapper.hpp"
#include "../py_objectivefunction.hpp"
#include "../../NudgedElasticBand.h"
// clang-format on

template <class NEBObjectiveFunctionBase = NEBObjectiveFunction>
class PyNEBObjectiveFunction : public PyObjectiveFunction<NEBObjectiveFunctionBase> {
public:
    using PyObjectiveFunction<NEBObjectiveFunctionBase>::PyObjectiveFunction; // Inherit
                                                                                 // constructor
    // Override pure virtual with non-pure
    double getEnergy() override {
        PYBIND11_OVERRIDE(double, NEBObjectiveFunctionBase, getEnergy, );
    };
    double getConvergence() override {
        PYBIND11_OVERRIDE(double, NEBObjectiveFunctionBase, getConvergence, );
    };
    bool isConverged() override {
        PYBIND11_OVERRIDE(bool, NEBObjectiveFunctionBase, isConverged, );
    };
    int degreesOfFreedom() override {
        PYBIND11_OVERRIDE(int, NEBObjectiveFunctionBase, degreesOfFreedom, );
    };
    VectorXd getPositions() override {
        PYBIND11_OVERRIDE(VectorXd, NEBObjectiveFunctionBase, getPositions, );
    };
    VectorXd getGradient(bool fdstep = false) override {
        PYBIND11_OVERRIDE(VectorXd, NEBObjectiveFunctionBase, getGradient, fdstep = false);
    };
    VectorXd difference(VectorXd a, VectorXd b) override {
        PYBIND11_OVERRIDE(VectorXd, NEBObjectiveFunctionBase, difference, a, b);
    };
    void setPositions(VectorXd x) override {
        PYBIND11_OVERRIDE(void, NEBObjectiveFunctionBase, setPositions, x);
    };
};

void py_nebobjfunc(py::module_ &m) {
    py::class_<NEBObjectiveFunction, ObjectiveFunction, PyNEBObjectiveFunction<>>(
        m, "NEBObjectiveFunction")
        /*
        ** Constructors
        */
        .def(py::init<NudgedElasticBand* /*nebPassed*/, Parameters * /*parametersPassed*/>())

        /*
        ** Methods
        */
        .def("getEnergy", &NEBObjectiveFunction::getEnergy)
        .def("getConvergence", &NEBObjectiveFunction::getConvergence)
        .def("isConverged", &NEBObjectiveFunction::isConverged)
        .def("degreesOfFreedom", &NEBObjectiveFunction::degreesOfFreedom)
        .def("getPositions", &NEBObjectiveFunction::getPositions)
        .def("getGradient", &NEBObjectiveFunction::getGradient, py::arg("fdstep") = false)
        .def("difference", &NEBObjectiveFunction::difference, py::arg("a"), "b"_a)
        .def("setPositions", &NEBObjectiveFunction::setPositions)

        /*
        ** Python helpers
        */

        .def("__repr__",
             [](const NEBObjectiveFunction &a) { return "<NEBObjectiveFunction object>"; });
}

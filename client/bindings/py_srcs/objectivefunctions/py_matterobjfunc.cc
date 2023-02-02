// clang-format off
#include "../py_wrapper.hpp"
#include "../py_objectivefunction.hpp"
// clang-format on

template <class MatterObjectiveFunctionBase = MatterObjectiveFunction>
class PyMatterObjectiveFunction : public PyObjectiveFunction<MatterObjectiveFunctionBase> {
public:
    using PyObjectiveFunction<MatterObjectiveFunctionBase>::PyObjectiveFunction; // Inherit
                                                                                 // constructor
    // Override pure virtual with non-pure
    double getEnergy() override {
        PYBIND11_OVERRIDE(double, MatterObjectiveFunctionBase, getEnergy, );
    };
    double getConvergence() override {
        PYBIND11_OVERRIDE(double, MatterObjectiveFunctionBase, getConvergence, );
    };
    bool isConverged() override {
        PYBIND11_OVERRIDE(bool, MatterObjectiveFunctionBase, isConverged, );
    };
    int degreesOfFreedom() override {
        PYBIND11_OVERRIDE(int, MatterObjectiveFunctionBase, degreesOfFreedom, );
    };
    VectorXd getPositions() override {
        PYBIND11_OVERRIDE(VectorXd, MatterObjectiveFunctionBase, getPositions, );
    };
    VectorXd getGradient(bool fdstep = false) override {
        PYBIND11_OVERRIDE(VectorXd, MatterObjectiveFunctionBase, getGradient, fdstep = false);
    };
    VectorXd difference(VectorXd a, VectorXd b) override {
        PYBIND11_OVERRIDE(VectorXd, MatterObjectiveFunctionBase, difference, a, b);
    };
    void setPositions(VectorXd x) override {
        PYBIND11_OVERRIDE(void, MatterObjectiveFunctionBase, setPositions, x);
    };
};

void py_matterobjfunc(py::module_ &m) {
    py::class_<MatterObjectiveFunction, ObjectiveFunction, PyMatterObjectiveFunction<>>(
        m, "MatterObjectiveFunction")
        /*
        ** Constructors
        */
        .def(py::init<Matter * /*matterPassed*/, Parameters * /*parametersPassed*/>())

        /*
        ** Methods
        */
        .def("getEnergy", &MatterObjectiveFunction::getEnergy)
        .def("getConvergence", &MatterObjectiveFunction::getConvergence)
        .def("isConverged", &MatterObjectiveFunction::isConverged)
        .def("degreesOfFreedom", &MatterObjectiveFunction::degreesOfFreedom)
        .def("getGradient", &MatterObjectiveFunction::getGradient, py::arg("fdstep") = false)
        .def("difference", &MatterObjectiveFunction::difference, py::arg("a"), "b"_a)
        .def_property("positions",
                      &MatterObjectiveFunction::getPositions,
                      &MatterObjectiveFunction::setPositions)

        /*
        ** Python helpers
        */

        .def("__repr__",
             [](const MatterObjectiveFunction &a) { return "<MatterObjectiveFunction object>"; });
}

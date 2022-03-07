// clang-format off
#include "../py_optimizer.hpp"
// Binding
#include "../../LBFGS.h"
// clang-format on

template <class LBFGSBase = LBFGS>
class PyLBFGS : public PyOptimizer<LBFGSBase> {
public:
    using PyOptimizer<LBFGSBase>::PyOptimizer; // Inherit
    // Override pure virtual with non-pure
    int step(double maxMove) override {
        PYBIND11_OVERRIDE(int, LBFGSBase, step, maxMove);
    };
    int run(int maxIterations, double maxMove) override {
        PYBIND11_OVERRIDE(int, LBFGSBase, run, maxIterations, maxMove);
    };
};

// TODO: Generalize, see Prior in gprd
void py_lbfgs(py::module_ &m) {
    py::class_<LBFGS, Optimizer, PyLBFGS<>>(
        m, "LBFGS")
        /*
        ** Constructors
        */
        .def(py::init<ObjectiveFunction* /*objf*/, Parameters * /*parameters*/>())

        /*
        ** Methods
        */
        .def("step", &LBFGS::step, py::arg("maxMove"))
        .def("run", &LBFGS::run, py::arg("maxIterations"), "maxMove"_a)
        .def("update", &LBFGS::update, py::arg("r1"), "r0"_a, "f1"_a, "f0"_a)
        .def("reset", &LBFGS::reset)

        /*
        ** Python helpers
        */

        .def("__repr__",
             [](const LBFGS &a) { return "<LBFGS optimizer>"; });
}

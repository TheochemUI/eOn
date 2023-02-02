#include "py_optimizer.hpp"

void py_optimizer(py::module_ &m) {
    py::class_<Optimizer, PyOptimizer<> >(m, "Optimizer")
        .def(py::init())
        /*
        ** Functions
        */
        .def_static("getOptimizer", &Optimizer::getOptimizer)
        .def("step", &Optimizer::step, py::arg("maxMove"))
        .def("run", &Optimizer::run, py::arg("maxIterations"), "maxMove"_a)

        /*
        ** Python helpers
        */

        .def("__repr__", [](const Optimizer &a) { return "<Optimizer object>"; });
}

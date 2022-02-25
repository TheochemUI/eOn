// clang-format off
#include "py_wrapper.hpp"
#include "../Log.h"
// Additional
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
// clang-format on

void py_log(py::module_ &m) {
    // TODO: Make this a class so this can be initialized multiple times
    m.def("log_init", &log_init, py::arg("p"), "filename"_a);
    m.def("log_init", &log);
    m.def("log_close", &log_close);
}

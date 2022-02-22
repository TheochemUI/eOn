#ifndef PY_WRAPPER_H
#define PY_WRAPPER_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <fstream>
#include <iostream>

using namespace std::string_literals; // For ""s
using namespace pybind11::literals;   // For ""_a
namespace py = pybind11;              // Convention

void py_parameters(py::module_ &m);
void py_potential(py::module_ &m);

#endif /* PY_WRAPPER_H */

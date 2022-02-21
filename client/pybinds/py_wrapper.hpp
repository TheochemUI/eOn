#ifndef PY_WRAPPER_H
#define PY_WRAPPER_H

#include <pybind11/pybind11.h>
#include <fstream>
#include <iostream>

using namespace std::string_literals; // For ""s
namespace py = pybind11; // Convention

void py_parameters(py::module_ &);

#endif /* PY_WRAPPER_H */

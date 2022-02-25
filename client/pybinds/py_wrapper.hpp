#ifndef PY_WRAPPER_H
#define PY_WRAPPER_H

// Standard libraries
#include <fstream>
#include <iostream>
#include <string>

// Bindings
#include <pybind11/pybind11.h>
// TODO: Setup bindings for this
PYBIND11_MAKE_OPAQUE(std::vector<std::string>)
// Basics
#include "../Matter.h"
#include "../Parameters.h"

// Namespaces
using namespace std::string_literals; // For ""s
using namespace pybind11::literals;   // For ""_a
namespace py = pybind11;              // Convention

// Forward declarations
void py_parameters(py::module_ &m);
void py_matter(py::module_ &m);
void py_log(py::module_ &m);
// Jobs
void py_job(py::module_ &m);
void py_saddlesearchjob(py::module_ &m);
// Optimizers
void py_optimizer(py::module_ &m);
// Objective Functions
void py_objectivefunction(py::module_ &m);
// Potentials
void py_potential(py::module_ &m);
void py_morse(py::module_ &m);

#endif /* PY_WRAPPER_H */

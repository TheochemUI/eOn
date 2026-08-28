/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/

#include "ASE.h"
#include "../../PyGuard.h"
#include "../../fpe_handler.h"
#include <cstdlib> // for exit()
#include <pybind11/embed.h>
#include <pybind11/numpy.h> // for py::array_t
#include <pybind11/pybind11.h>
#include <string>
#include <tuple>
#include <vector>

namespace py = pybind11;

ASE::ASE(const Parameters &a_params)
    : Potential(PotType::ASE_POT, a_params) {
  eonc::ensure_interpreter();
  counter = 1;
  std::string py_file = a_params.potential_options.extPotPath;

  // import
  try {
    // must briefly disable FPE because Python packages like Numpy causes it
    // during import
    eonc::FPEHandler fpeh;
    fpeh.eat_fpe();

    // Create a Python script to use importlib.util to load the module
    py::exec(R"(
            import sys
            import importlib.util

            def load_module_from_path(module_name, file_path):
                spec = importlib.util.spec_from_file_location(module_name, file_path)
                module = importlib.util.module_from_spec(spec)
                sys.modules[module_name] = module
                spec.loader.exec_module(module)
                return module
        )");

    // Prepare the module name and file path
    std::string module_name = "ase_eon";
    py::object load_module = py::globals()["load_module_from_path"];
    py_module = load_module(module_name, py_file);

    fpeh.restore_fpe();

    calculator = py_module.attr("ase_calc")();
    _calculate = py_module.attr("_calculate");
    if (py::hasattr(py_module, "batch_calculate")) {
      batch_calculate = py_module.attr("batch_calculate");
    }

  } catch (const std::exception &e) {
    fprintf(stderr,
            "ASE Calculator: Exception during Python module import: %s\n",
            e.what());
    fprintf(stderr, "%s should exist and have no errors on the Python side.\n",
            py_file.c_str());
    exit(1);
  }
  return;
}

void ASE::force(long nAtoms, const double *R, const int *atomicNrs, double *F,
                double *U, double *variance, const double *box) {
  variance = nullptr;
  try {
    // TODO(rg): This is easier on the type system if Eigen::Map is used like in
    // ASE_ORCA convert arrays to Numpy arrays
    std::vector<size_t> R_shape = {static_cast<size_t>(nAtoms), 3};
    py::array_t<double> R_np(R_shape, R);
    py::array_t<int> atomicNrs_np(nAtoms, atomicNrs);
    py::array_t<double> box_np({3, 3}, box);

    // get energy and forces (in this order) from Python
    std::tuple<double, py::array_t<double>> py_result =
        _calculate(R_np, atomicNrs_np, box_np, calculator)
            .cast<std::tuple<double, py::array_t<double>>>();

    // copy the results to the output arrays
    *U = std::get<0>(py_result);
    py::array_t<double> forces = std::get<1>(py_result);
    auto buffer = forces.request();
    double *ptr = static_cast<double *>(buffer.ptr);
    std::copy(ptr, ptr + buffer.size, F);

  } catch (py::error_already_set &e) {
    fprintf(stderr, "ASE calculator: Python error: %s\n", e.what());
    exit(1);
  } catch (const std::exception &e) {
    fprintf(stderr, "ASE calculator: C++ exception: %s\n", e.what());
    exit(1);
  }

  counter++;
  return;
}

void ASE::forceBatch(long nSystems, long nAtoms, const double *const *R,
                     const int *const *atomicNrs, double *const *forces,
                     double *energies, double *variances,
                     const double *const *boxes) {
  if (!batch_calculate) {
    for (int i = 0; i < nSystems; ++i) {
      force(nAtoms, R[i], atomicNrs[i], forces[i], &energies[i], &variances[i],
            boxes[i]);
    }
    return;
  }

  try {
    variances = nullptr;

    std::vector<double> R_data(static_cast<size_t>(nSystems) *
                               static_cast<size_t>(nAtoms) * 3);
    std::vector<int> atomicNrs_data(static_cast<size_t>(nSystems) *
                                    static_cast<size_t>(nAtoms));
    std::vector<double> boxes_data(static_cast<size_t>(nSystems) * 9);

    for (long i = 0; i < nSystems; ++i) {
      std::copy(R[i], R[i] + nAtoms * 3, R_data.begin() + i * nAtoms * 3);
      std::copy(atomicNrs[i], atomicNrs[i] + nAtoms,
                atomicNrs_data.begin() + i * nAtoms);
      std::copy(boxes[i], boxes[i] + 9, boxes_data.begin() + i * 9);
    }

    std::vector<size_t> R_shape = {static_cast<size_t>(nSystems),
                                   static_cast<size_t>(nAtoms), 3};
    py::array_t<double> R_np(R_shape, R_data.data());

    std::vector<size_t> atomicNrs_shape = {static_cast<size_t>(nSystems),
                                           static_cast<size_t>(nAtoms)};
    py::array_t<int> atomicNrs_np(atomicNrs_shape, atomicNrs_data.data());

    std::vector<size_t> boxes_shape = {static_cast<size_t>(nSystems), 3, 3};
    py::array_t<double> boxes_np(boxes_shape, boxes_data.data());

    std::tuple<py::array_t<double>, py::array_t<double>> py_result =
        (*batch_calculate)(R_np, atomicNrs_np, boxes_np, calculator)
            .cast<std::tuple<py::array_t<double>, py::array_t<double>>>();

    // copy the results to the output arrays
    py::array_t<double> E = std::get<0>(py_result);
    auto buffer_E = E.request();
    double *ptr_E = static_cast<double *>(buffer_E.ptr);
    std::copy(ptr_E, ptr_E + buffer_E.size, energies);

    py::array_t<double> F = std::get<1>(py_result);
    auto buffer_F = F.request();
    double *ptr_F = static_cast<double *>(buffer_F.ptr);
    for (long i = 0; i < nSystems; ++i) {
      std::copy(ptr_F + i * nAtoms * 3, ptr_F + (i + 1) * nAtoms * 3,
                forces[i]);
    }
  } catch (py::error_already_set &e) {
    fprintf(stderr, "ASE calculator: Python error: %s\n", e.what());
    exit(1);
  } catch (const std::exception &e) {
    fprintf(stderr, "ASE calculator: C++ exception: %s\n", e.what());
    exit(1);
  }

  counter += nSystems;
  return;
}

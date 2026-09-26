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

#include "eon/potentials/ASE/ASE.h"
#include "eon/Eigen.h"
#include "eon/EonLogger.h"
#include "eon/Parameters.h"
#include "eon/PyGuard.h"
#include "eon/fpe_handler.h"
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace py = pybind11;

ASE::ASE(const eonc::Parameters &a_params)
    : eonc::Potential(eonc::PotType::ASE_POT, a_params) {
  eonc::ensure_interpreter();
  counter = 1;
  std::string py_file = a_params.potential_options().extPotPath;

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
    has_batch_ = py::hasattr(py_module, "batch_calculate");
    if (has_batch_) {
      batch_calculate_ = py_module.attr("batch_calculate");
    }

  } catch (const std::exception &e) {
    EONC_LOG_ERROR("ASE calculator import failed for {}: {}", py_file,
                   e.what());
    throw std::runtime_error(std::string("ASE calculator import failed: ") +
                             e.what());
  }
  return;
}

void ASE::force(long nAtoms, const double *R, const int *atomicNrs, double *F,
                double *U, double *variance, const double *box) {
  if (variance != nullptr) {
    *variance = 0.0;
  }
  py::gil_scoped_acquire gil;
  try {
    const Eigen::Map<const AtomMatrix> positions(R, nAtoms, 3);
    const Eigen::Map<const RotationMatrix> boxx(box);
    const Eigen::Map<const Eigen::VectorXi> atmnmrs(atomicNrs, nAtoms);

    std::tuple<double, py::array_t<double>> py_result =
        _calculate(positions, atmnmrs, boxx, calculator)
            .cast<std::tuple<double, py::array_t<double>>>();

    *U = std::get<0>(py_result);
    py::array_t<double> forces = std::get<1>(py_result);
    auto buffer = forces.request();
    if (buffer.size < nAtoms * 3) {
      throw std::runtime_error(
          "ASE _calculate returned forces of the wrong size");
    }
    Eigen::Map<AtomMatrix>(F, nAtoms, 3) = Eigen::Map<const AtomMatrix>(
        static_cast<const double *>(buffer.ptr), nAtoms, 3);

  } catch (py::error_already_set &e) {
    EONC_LOG_ERROR("ASE calculator Python error: {}", e.what());
    throw std::runtime_error(std::string("ASE calculator Python error: ") +
                             e.what());
  } catch (const std::exception &e) {
    EONC_LOG_ERROR("ASE calculator C++ exception: {}", e.what());
    throw std::runtime_error(std::string("ASE calculator C++ exception: ") +
                             e.what());
  }

  counter++;
  return;
}

void ASE::forceBatch(long nSystems, long nAtoms, const double *const *positions,
                     const int *const *atomicNrs, double *const *forces,
                     double *energies, double *variances,
                     const double *const *boxes) {
  if (!has_batch_) {
    eonc::Potential::forceBatch(nSystems, nAtoms, positions, atomicNrs, forces,
                                energies, variances, boxes);
    return;
  }
  py::gil_scoped_acquire gil;
  try {
    const auto ns = static_cast<size_t>(nSystems);
    const auto na = static_cast<size_t>(nAtoms);
    std::vector<double> R_data(ns * na * 3);
    std::vector<int> Z_data(ns * na);
    std::vector<double> box_data(ns * 9);
    for (long s = 0; s < nSystems; ++s) {
      std::copy(positions[s], positions[s] + static_cast<long>(na) * 3,
                R_data.data() + static_cast<size_t>(s) * na * 3);
      std::copy(atomicNrs[s], atomicNrs[s] + nAtoms,
                Z_data.data() + static_cast<size_t>(s) * na);
      std::copy(boxes[s], boxes[s] + 9,
                box_data.data() + static_cast<size_t>(s) * 9);
    }
    py::array_t<double> R_np({ns, na, size_t{3}}, R_data.data());
    py::array_t<int> Z_np({ns, na}, Z_data.data());
    py::array_t<double> box_np({ns, size_t{3}, size_t{3}}, box_data.data());
    auto py_result =
        batch_calculate_(R_np, Z_np, box_np, calculator)
            .cast<std::tuple<py::array_t<double>, py::array_t<double>>>();
    py::array_t<double> E = std::get<0>(py_result);
    py::array_t<double> F = std::get<1>(py_result);
    auto bufE = E.request();
    auto bufF = F.request();
    if (bufE.size < nSystems || bufF.size < nSystems * nAtoms * 3) {
      throw std::runtime_error(
          "ASE batch_calculate returned energies/forces of the wrong size");
    }
    auto *ePtr = static_cast<double *>(bufE.ptr);
    auto *fPtr = static_cast<double *>(bufF.ptr);
    std::copy(ePtr, ePtr + nSystems, energies);
    for (long s = 0; s < nSystems; ++s) {
      std::copy(fPtr + s * nAtoms * 3, fPtr + (s + 1) * nAtoms * 3, forces[s]);
      if (variances) {
        variances[s] = 0.0;
      }
    }
  } catch (py::error_already_set &e) {
    EONC_LOG_ERROR("ASE calculator Python error: {}", e.what());
    throw std::runtime_error(std::string("ASE calculator Python error: ") +
                             e.what());
  } catch (const std::exception &e) {
    EONC_LOG_ERROR("ASE calculator C++ exception: {}", e.what());
    throw std::runtime_error(std::string("ASE calculator C++ exception: ") +
                             e.what());
  }
  counter += static_cast<size_t>(nSystems);
}

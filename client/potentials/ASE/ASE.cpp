//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "ASE.h"
#include "../../ExceptionsEON.h"
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>  // for py::array_t
#include <vector>
#include <string>
#include <tuple>
#include <cstdlib>  // for exit()


namespace py = pybind11;

ASE::ASE(Parameters *p) : guard{} {
    parameters = p;
    
    // strip the .py extension
    std::string py_file = parameters->ASEScript;
    const std::string extension = ".py";
    size_t pos = py_file.rfind(extension);
    if (pos != std::string::npos && pos == py_file.length() - extension.length()) {
        // Erase the ".py" extension
        py_file.erase(pos, extension.length());
    }
    
    // import
    try {
        // must briefly disable FPE because Python packages like Numpy causes it during import
        bool FPE_WAS_ENABLED = isFPEEnabled();
        if (FPE_WAS_ENABLED) {
            disableFPE();
        }
        py_module = py::module_::import(py_file.c_str());
        if (FPE_WAS_ENABLED) {
            enableFPE();
        }

        calculator = py_module.attr("ase_calc")();
        _calculate = py_module.attr("_calculate");

    } catch (const std::exception &e) {
        fprintf(stderr, "ASE Calculator: Exception during Python module import: %s\n", e.what());
        fprintf(stderr, "%s.py should be a valid script in the same directory"
                         " as where eonclient is run.\n", py_file.c_str());
        exit(1);
    }
    return;
}

void ASE::cleanMemory(void){
    return;
}

ASE::~ASE()
{
    cleanMemory();
}


// number of atoms, pointer to array of positions, pointer to array of atomic numbers,
// pointer to array of forces, pointer to internal energy,
// pointer to cell array
void ASE::force(long N, const double *R, const int *atomicNrs,
                double *F, double *U, const double *box) {
    try {
        // convert arrays to Numpy arrays
        std::vector<size_t> R_shape = {static_cast<size_t>(N), 3};
        py::array_t<double> R_np(R_shape, R);
        py::array_t<int> atomicNrs_np(N, atomicNrs);
        py::array_t<double> box_np({3, 3}, box);

        // get energy and forces (in this order) from Python
        std::tuple<double, py::array_t<double>> py_result =
            _calculate(R_np, atomicNrs_np, box_np, calculator)
            .cast<std::tuple<double, py::array_t<double>>>();

        // copy the results to the output arrays
        *U = std::get<0>(py_result);
        py::array_t<double> forces = std::get<1>(py_result);
        auto buffer = forces.request();
        double* ptr = static_cast<double*>(buffer.ptr);
        std::copy(ptr, ptr + buffer.size, F);

    } catch (py::error_already_set& e) {
        fprintf(stderr, "ASE calculator: Python error: %s\n", e.what());
        exit(1);
    } catch (const std::exception& e) {
        fprintf(stderr, "ASE calculator: C++ exception: %s\n", e.what());
        exit(1);
    }

    return;
}


//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef PYSURROGATE_INTERFACE
#define PYSURROGATE_INTERFACE

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include "../../Potential.h"

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

namespace py = pybind11;
using namespace pybind11::literals; // to bring in the `_a` literal

class PySurrogate : public Potential {

private:
  // Variables
  double fake1;
  double fake2;
  py::object hpfit;
  py::object gpmod;

public:
  // Functions
  // constructor and destructor
  PySurrogate(Parameters *p) : Potential(p), fake1{0}, fake2{0} {
    py::module_ sys = py::module_::import("sys");
    py::exec(
        R"(sys.path.insert(0, "/home/rgoswami/Git/Github/Python/DTU_CatLearn"))");

    // Import the required modules
    py::module_ hpfitter_module = py::module_::import("catlearn.regression.gaussianprocess.hpfitter");
    py::module_ objectfunctions_module = py::module_::import("catlearn.regression.gaussianprocess.objectfunctions.factorized_likelihood");
    py::module_ optimizers_module = py::module_::import("catlearn.regression.gaussianprocess.optimizers");
    py::module_ gp_module = py::module_::import("catlearn.regression.gaussianprocess.gp.gp");
    py::module_ prior_median_module = py::module_::import("catlearn.regression.gaussianprocess.means.median");

    // Get the classes from the imported modules
    py::object hpfitter_class = hpfitter_module.attr("HyperparameterFitter");
    py::object objectfunctions_class = objectfunctions_module.attr("FactorizedLogLikelihood");
    py::object optimizers_class = optimizers_module.attr("run_golden");
    py::object line_search_scale_class = optimizers_module.attr("line_search_scale");
    py::object gp_class = gp_module.attr("GaussianProcess");
    py::object prior_median_class = prior_median_module.attr("Prior_median");

    // Create the objects and set the arguments
    py::dict local_kwargs;
    local_kwargs["tol"] = 1e-5;
    local_kwargs["optimize"] = true;
    local_kwargs["multiple_max"] = true;

    py::dict kwargs_optimize;
    kwargs_optimize["local_run"] = optimizers_class;
    kwargs_optimize["maxiter"] = 5000;
    kwargs_optimize["jac"] = false;
    kwargs_optimize["bounds"] = py::none(); // None
    kwargs_optimize["ngrid"] = 80;
    kwargs_optimize["use_bounds"] = true;
    kwargs_optimize["local_kwargs"] = local_kwargs;

    this->hpfit = hpfitter_class(objectfunctions_class(),
                                      line_search_scale_class,
                                      kwargs_optimize,
                                      true); // distance_matrix = true

    this->gpmod = gp_class(py::arg("prior") = prior_median_class(),
                                 py::arg("use_derivatives") = true,
                                 py::arg("hpfitter") = hpfit);
};

  // To satisfy interface
  void cleanMemory(void);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, const double *box) override;
};
#endif

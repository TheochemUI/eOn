//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "CatLearnPot.h"
#include "Eigen/src/Core/Matrix.h"

CatLearnPot::CatLearnPot(shared_ptr<Parameters> a_params)
    : Potential(PotType::CatLearn, a_params) {
  py::module_ sys = py::module_::import("sys");
  py::exec(fmt::format("sys.path.insert(0, {})", a_params->catl_path));

  // Import the required modules
  py::module np = py::module::import("numpy");
  py::module_ hpfitter_module =
      py::module_::import("catlearn.regression.gaussianprocess.hpfitter");
  py::module_ objectfunctions_module =
      py::module_::import("catlearn.regression.gaussianprocess."
                          "objectfunctions.factorized_likelihood");
  py::module_ optimizers_module =
      py::module_::import("catlearn.regression.gaussianprocess.optimizers");
  py::module_ gp_module =
      py::module_::import("catlearn.regression.gaussianprocess.gp.gp");
  py::module_ prior_max_module =
      py::module_::import("catlearn.regression.gaussianprocess.means.max");
  py::module_ _kernel =
      py::module_::import("catlearn.regression.gaussianprocess.kernel.se");
  py::module_ normal_module = py::module_::import(
      "catlearn.regression.gaussianprocess.pdistributions.normal");

  // Get the classes from the imported modules
  py::object hpfitter_class = hpfitter_module.attr("HyperparameterFitter");
  py::object objectfunctions_class =
      objectfunctions_module.attr("FactorizedLogLikelihood");
  py::object optimizers_class = optimizers_module.attr("run_golden");
  py::object line_search_scale_class =
      optimizers_module.attr("line_search_scale");
  py::object gp_class = gp_module.attr("GaussianProcess");
  py::object prior_max_class = prior_max_module.attr("Prior_max");
  py::object kernel_class = _kernel.attr("SE_Derivative");

  // Create the objects and set the arguments
  py::dict local_kwargs;
  local_kwargs["tol"] = 1e-5;
  local_kwargs["optimize"] = true;
  local_kwargs["multiple_max"] = true;

  py::dict kwargs_optimize;
  kwargs_optimize["local_run"] = optimizers_class;
  kwargs_optimize["maxiter"] = 1000;
  kwargs_optimize["jac"] = false;
  kwargs_optimize["bounds"] = py::none(); // None
  kwargs_optimize["ngrid"] = 80;
  kwargs_optimize["use_bounds"] = true;
  kwargs_optimize["local_kwargs"] = local_kwargs;

  // Hyperparameter Prior
  py::object normal_class = normal_module.attr("Normal_prior");
  py::list lenlist, noiselist;
  lenlist.append(normal_class(0.0, 2.0));
  noiselist.append(normal_class(-9.0, 2.0));
  _prior["length"] = np.attr("array")(lenlist);
  _prior["noise"] = np.attr("array")(noiselist);

  // Hyperparameter Fitter
  this->hpfit = hpfitter_class(objectfunctions_class(), line_search_scale_class,
                               py::arg("opt_kwargs") = kwargs_optimize,
                               py::arg("distance_matrix") = true);
  // Kernel
  this->kernel = kernel_class(py::arg("use_fingerprint") = false);
  // GP Model
  this->gpmod =
      gp_class(py::arg("prior") = prior_max_class(), py::arg("kernel") = kernel,
               py::arg("use_derivatives") = true, py::arg("hpfitter") = hpfit);
};

void CatLearnPot::train_optimize(Eigen::MatrixXd features,
                                 Eigen::MatrixXd targets) {

  gpmod.attr("optimize")(features, targets, py::arg("retrain") = true,
                         py::arg("prior") = _prior);
  return;
}
void CatLearnPot::cleanMemory(void) { return; }

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
void CatLearnPot::force(long N, const double *R, const int *atomicNrs,
                        double *F, double *U, double *variance,
                        const double *box) {
  Eigen::MatrixXd features =
      Eigen::Map<Eigen::MatrixXd>(const_cast<double *>(R), 1, N * 3);
  py::tuple ef_and_unc =
      (this->gpmod.attr("predict")(features, "get_variance"_a = true));
  auto ef_dat = ef_and_unc[0].cast<Eigen::MatrixXd>();
  auto vari = ef_and_unc[1].cast<Eigen::MatrixXd>();
  auto gradients = ef_dat.block(0, 1, 1, N * 3);
  for (int idx = 0; idx < N; idx++) {
    F[3 * idx] = gradients(0, 3 * idx) * -1;
    F[3 * idx + 1] = gradients(0, 3 * idx + 1) * -1;
    F[3 * idx + 2] = gradients(0, 3 * idx + 2) * -1;
  }
  for (int idx = 0; idx < 1 + (N * 3); idx++) {
    variance[idx] = vari(0, idx);
  }
  *U = ef_dat(0, 0);
  // SPDLOG_TRACE("Energy and Forces: {}\nVariance: {}", fmt::streamed(ef_dat),
  // fmt::streamed(vari));
  return;
}

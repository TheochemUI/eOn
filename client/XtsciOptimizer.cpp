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
#include "eon/XtsciOptimizer.h"
#include "eon/PairHessian.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <xts.h>

#if XTS_ABI_VERSION_MINOR < 5
#error                                                                         \
    "xtsci-optimize ABI minor 5 or newer is required for atom maxmove and rematch LBFGS policy"
#endif

namespace {

struct XtsObjectiveContext {
  eonc::ObjectiveFunction *objective;
  eonc::PotType pot;
  std::string precon;
  double precon_A;
  double precon_mu;
  double precon_rcut;
};

Eigen::Map<const Eigen::VectorXd>
map_input(const DLManagedTensorVersioned *tensor) {
  const auto &dl = tensor->dl_tensor;
  if (dl.ndim != 1 || dl.dtype.code != kDLFloat || dl.dtype.bits != 64 ||
      dl.dtype.lanes != 1 || dl.device.device_type != kDLCPU ||
      dl.shape == nullptr || dl.data == nullptr) {
    throw std::runtime_error("xtsci objective requires a CPU f64 vector");
  }
  return {static_cast<const double *>(dl.data) +
              dl.byte_offset / sizeof(double),
          static_cast<Eigen::Index>(dl.shape[0])};
}

Eigen::Map<Eigen::VectorXd> map_output(DLManagedTensorVersioned *tensor) {
  const auto &dl = tensor->dl_tensor;
  if (dl.ndim != 1 || dl.dtype.code != kDLFloat || dl.dtype.bits != 64 ||
      dl.dtype.lanes != 1 || dl.device.device_type != kDLCPU ||
      dl.shape == nullptr || dl.data == nullptr) {
    throw std::runtime_error("xtsci tensor requires a CPU f64 vector");
  }
  return {static_cast<double *>(dl.data) + dl.byte_offset / sizeof(double),
          static_cast<Eigen::Index>(dl.shape[0])};
}

// Matter::setPositions always dirties the PES cache. Energy and
// forces come from one potential->force() call, so eval then grad at
// the same x must not setPositions twice.
void set_positions_if_changed(eonc::ObjectiveFunction *obj,
                              const Eigen::VectorXd &x) {
  const auto cur = obj->getPositions();
  if (cur.size() == x.size() && (cur - x).isZero(0.0)) {
    return;
  }
  obj->setPositions(x);
}

// One Matter::computePotential. Energy and forces share that call.
xts_status_t evaluate_gradient(void *user, const DLManagedTensorVersioned *x,
                               double *value_out,
                               DLManagedTensorVersioned *gradient_out) {
  try {
    auto *context = static_cast<XtsObjectiveContext *>(user);
    const auto positions = map_input(x);
    set_positions_if_changed(context->objective, positions);
    *value_out = context->objective->getEnergy();
    const auto gradient_value = context->objective->getGradient();
    auto output = map_output(gradient_out);
    if (output.size() != gradient_value.size()) {
      return XTS_INVALID_PARAMETER;
    }
    output = gradient_value;
    return XTS_SUCCESS;
  } catch (...) {
    return XTS_INTERNAL_ERROR;
  }
}

xts_status_t hessian(void *user, const DLManagedTensorVersioned *x,
                     DLManagedTensorVersioned *hess_out) {
  try {
    auto *context = static_cast<XtsObjectiveContext *>(user);
    const auto positions = map_input(x);
    const Eigen::MatrixXd H = eonc::pairhess::build(
        positions, context->precon, context->pot, context->precon_A,
        context->precon_mu, context->precon_rcut, *context->objective);
    auto output = map_output(hess_out);
    const Eigen::Index n = H.rows();
    if (output.size() != n * n) {
      return XTS_INVALID_PARAMETER;
    }
    for (Eigen::Index i = 0; i < n; ++i) {
      for (Eigen::Index j = 0; j < n; ++j) {
        output[i * n + j] = H(i, j);
      }
    }
    return XTS_SUCCESS;
  } catch (...) {
    return XTS_INTERNAL_ERROR;
  }
}

xts_method_t method_from_name(const std::string &name) {
  if (name == "lbfgs") {
    return XTS_LBFGS;
  }
  if (name == "bfgs") {
    return XTS_BFGS;
  }
  if (name == "sr1") {
    return XTS_SR1;
  }
  if (name == "sr2") {
    return XTS_SR2;
  }
  if (name == "newton") {
    return XTS_NEWTON;
  }
  if (name == "rfo") {
    return XTS_RFO;
  }
  if (name == "steepest") {
    return XTS_STEEPEST;
  }
  if (name == "adam") {
    return XTS_ADAM;
  }
  if (name == "pso") {
    return XTS_PSO;
  }
  if (name == "polak_ribiere" || name == "nlcg" || name == "pr") {
    return XTS_POLAK_RIBIERE;
  }
  if (name == "fletcher_reeves") {
    return XTS_FLETCHER_REEVES;
  }
  if (name == "hestenes_stiefel") {
    return XTS_HESTENES_STIEFEL;
  }
  if (name == "dai_yuan") {
    return XTS_DAI_YUAN;
  }
  if (name == "conjugate_descent") {
    return XTS_CONJUGATE_DESCENT;
  }
  if (name == "hager_zhang") {
    return XTS_HAGER_ZHANG;
  }
  if (name == "liu_storey") {
    return XTS_LIU_STOREY;
  }
  if (name == "fr_pr") {
    return XTS_FR_PR;
  }
  throw std::invalid_argument(
      "unknown Xtsci.method '" + name +
      "' (lbfgs, bfgs, sr1, sr2, newton, rfo, steepest, adam, pso, "
      "polak_ribiere, fletcher_reeves, hestenes_stiefel, dai_yuan, "
      "conjugate_descent, hager_zhang, liu_storey, fr_pr)");
}

bool is_host_precon(const std::string &name) {
  return name == "pair" || name == "pair_abs" || name == "pair_full" ||
         name == "exp" || name == "c1" || name == "lindh";
}

// [Xtsci] qn_step / precon are the natural knobs. [LBFGS] lbfgs_step
// and lbfgs_precon stay on the native optimizer and fill in when the
// Xtsci fields are still at their defaults.
std::string
resolved_qn_step(const eonc::Parameters::optimizer_options_t &opts) {
  if (opts.xtsci.qn_step != "lbfgs") {
    return opts.xtsci.qn_step;
  }
  if (opts.lbfgs.step == "newton" || opts.lbfgs.step == "rfo") {
    return opts.lbfgs.step;
  }
  return opts.xtsci.qn_step;
}

std::string resolved_precon(const eonc::Parameters::optimizer_options_t &opts) {
  if (opts.xtsci.precon != "none") {
    return opts.xtsci.precon;
  }
  return opts.lbfgs.precon;
}

std::string resolved_accept(const eonc::Parameters::optimizer_options_t &opts) {
  if (opts.xtsci.accept != "none") {
    return opts.xtsci.accept;
  }
  if (opts.lbfgs.accept == "energy" || opts.lbfgs.accept == "nonmonotone") {
    return opts.lbfgs.accept;
  }
  return opts.xtsci.accept;
}

} // namespace

namespace eonc {

XtsciOptimizer::~XtsciOptimizer() { xts_solver_free(m_solver); }

void XtsciOptimizer::ensureSolver(double a_maxMove) {
  if (m_solver != nullptr) {
    return;
  }
  const auto dim = static_cast<size_t>(m_objf->degreesOfFreedom());
  if (dim == 0) {
    return;
  }
  const auto stamp = xts_abi_stamp();
  if (xts_abi_compatible(&stamp) == 0) {
    throw std::runtime_error("incompatible xtsci-optimize ABI");
  }
  const auto &lbfgs = m_optConfig.opts.lbfgs;
  xts_control_t control{
      static_cast<size_t>(std::max<long>(1, m_optConfig.opts.max_iterations)),
      m_optConfig.opts.converged_force,
      std::max(a_maxMove, 1.0e-12),
      static_cast<size_t>(std::max<long>(0, lbfgs.memory)),
      0.0,
  };
  m_solver = xts_solver_create(method_from_name(m_optConfig.opts.xtsci.method),
                               &control, dim);
  if (m_solver == nullptr) {
    throw std::runtime_error(xts_last_error());
  }
  const auto step = resolved_qn_step(m_optConfig.opts);
  if (step == "newton") {
    xts_solver_set_qn_step(m_solver, XTS_QN_NEWTON);
  } else if (step == "rfo") {
    xts_solver_set_qn_step(m_solver, XTS_QN_RFO);
  } else {
    xts_solver_set_qn_step(m_solver, XTS_QN_LBFGS);
  }
  const auto accept = resolved_accept(m_optConfig.opts);
  if (accept == "energy") {
    xts_solver_set_accept(m_solver, XTS_ACCEPT_ENERGY);
  } else if (accept == "nonmonotone") {
    xts_solver_set_accept(m_solver, XTS_ACCEPT_NONMONOTONE);
  } else {
    xts_solver_set_accept(m_solver, XTS_ACCEPT_NONE);
  }
  xts_solver_set_project_rigid(m_solver, lbfgs.project_rigid ? 1 : 0);
  xts_solver_set_extra_updates(
      m_solver, static_cast<size_t>(std::max<long>(0, lbfgs.extra_updates)));
  if (lbfgs.curvature == "cautious") {
    xts_solver_set_cautious(m_solver, lbfgs.cautious_eps, lbfgs.cautious_alpha);
  } else {
    xts_solver_set_cautious(m_solver, 0.0, lbfgs.cautious_alpha);
  }
}

int XtsciOptimizer::step(double a_maxMove) {
  if (m_objf->degreesOfFreedom() <= 0) {
    return m_objf->isConverged() ? 1 : 0;
  }
  ensureSolver(a_maxMove);
  if (m_solver == nullptr) {
    return m_objf->isConverged() ? 1 : 0;
  }
  // Native LBFGS clips by the largest per-atom move, not ||d||_2.
  xts_solver_set_maxmove(m_solver, 0.0);
  xts_solver_set_atom_maxmove(m_solver, std::max(a_maxMove, 0.0));

  auto positions = m_objf->getPositions();
  if (positions.size() != m_objf->degreesOfFreedom()) {
    throw std::runtime_error("xtsci objective position dimension mismatch");
  }
  auto *tensor = xts_tensor_borrow_cpu_f64(
      positions.data(), static_cast<size_t>(positions.size()));
  if (tensor == nullptr) {
    throw std::runtime_error("could not allocate xtsci objective tensor");
  }
  const auto &lbfgs = m_optConfig.opts.lbfgs;
  const auto precon = resolved_precon(m_optConfig.opts);
  XtsObjectiveContext context{
      m_objf.get(),   m_optConfig.potential, precon,
      lbfgs.precon_A, lbfgs.precon_mu,       lbfgs.precon_rcut,
  };
  xts_report_t report{};
  const xts_method_t method = method_from_name(m_optConfig.opts.xtsci.method);
  const bool want_hess =
      method == XTS_NEWTON || method == XTS_RFO || is_host_precon(precon);
  xts_status_t status;
  if (want_hess) {
    status = xts_solver_step_hess_fg(m_solver, evaluate_gradient, hessian,
                                     &context, tensor, &report);
  } else {
    status = xts_solver_step_fg(m_solver, evaluate_gradient, &context, tensor,
                                &report);
  }
  xts_tensor_free(tensor);
  if (status != XTS_SUCCESS) {
    throw std::runtime_error(xts_last_error());
  }
  set_positions_if_changed(m_objf.get(), positions);
  return m_objf->isConverged() ? 1 : 0;
}

int XtsciOptimizer::run(size_t a_maxIterations, double a_maxMove) {
  if (m_objf->degreesOfFreedom() <= 0 || a_maxIterations == 0) {
    return m_objf->isConverged() ? 1 : 0;
  }
  for (size_t i = 0; i < a_maxIterations; ++i) {
    if (step(a_maxMove) != 0) {
      return 1;
    }
  }
  return m_objf->isConverged() ? 1 : 0;
}

} // namespace eonc

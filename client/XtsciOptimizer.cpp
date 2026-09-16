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

#if XTS_ABI_VERSION_MINOR < 10
#error "xtsci-optimize ABI minor 10 or newer is required for set_periodic"
#endif

namespace {

struct XtsObjectiveContext {
  eonc::ObjectiveFunction *objective;
  eonc::PotType pot;
  std::string precon;
  double precon_A;
  double precon_mu;
  double precon_rcut;
  Eigen::VectorXd *cached_x;
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
                              const Eigen::VectorXd &x,
                              Eigen::VectorXd *cached) {
  if (cached != nullptr && cached->size() == x.size() &&
      (*cached - x).isZero(0.0)) {
    return;
  }
  // First host iteration: relaxMatter already called isConverged()
  // at this geometry. Do not dirty that cache.
  const auto cur = obj->getPositions();
  if (cur.size() == x.size() && (cur - x).isZero(0.0)) {
    if (cached != nullptr) {
      *cached = x;
    }
    return;
  }
  obj->setPositions(x);
  if (cached != nullptr) {
    *cached = x;
  }
}

// One Matter::computePotential. Energy and forces share that call.
xts_status_t evaluate_gradient(void *user, const DLManagedTensorVersioned *x,
                               double *value_out,
                               DLManagedTensorVersioned *gradient_out) {
  try {
    auto *context = static_cast<XtsObjectiveContext *>(user);
    const auto positions = map_input(x);
    set_positions_if_changed(context->objective, positions, context->cached_x);
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
  if (name == "fire") {
    return XTS_FIRE;
  }
  if (name == "bb" || name == "barzilai_borwein") {
    return XTS_BB;
  }
  if (name == "dogleg") {
    return XTS_DOGLEG;
  }
  if (name == "fire2") {
    return XTS_FIRE2;
  }
  throw std::invalid_argument(
      "unknown Xtsci.method '" + name +
      "' (lbfgs, bfgs, sr1, sr2, newton, rfo, steepest, adam, pso, "
      "polak_ribiere, fletcher_reeves, hestenes_stiefel, dai_yuan, "
      "conjugate_descent, hager_zhang, liu_storey, fr_pr, fire, bb, "
      "dogleg, fire2)");
}

bool is_host_precon(const std::string &name) {
  return name == "pair" || name == "pair_abs" || name == "pair_full" ||
         name == "exp" || name == "c1" || name == "lindh" ||
         name == "lindh_full" || name == "fischer" || name == "schlegel" ||
         name == "swart";
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
  const xts_method_t method = method_from_name(m_optConfig.opts.xtsci.method);
  double istep = std::max(a_maxMove, 1.0e-12);
  if (method == XTS_FIRE || method == XTS_FIRE2) {
    istep = m_optConfig.opts.time_step;
    if (istep <= 0.0) {
      istep = m_optConfig.opts.time_step_input;
    }
    if (istep <= 0.0) {
      istep = 0.1;
    }
  }
  xts_control_t control{
      static_cast<size_t>(std::max<long>(1, m_optConfig.opts.max_iterations)),
      m_optConfig.opts.converged_force,
      istep,
      static_cast<size_t>(std::max<long>(0, lbfgs.memory)),
      0.0,
  };
  m_solver = xts_solver_create(method, &control, dim);
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
  if (m_optConfig.opts.xtsci.highs) {
    if (xts_solver_set_highs(m_solver, 1) != 0) {
      throw std::runtime_error(
          "Xtsci.highs needs xtsci-optimize built with --features highs");
    }
  }
  const auto &mani = m_optConfig.opts.xtsci.manifold;
  if (mani == "so3" && dim != 9) {
    throw std::runtime_error(
        "Xtsci.manifold=so3 needs length 9; a 3N cluster is "
        "rigid_quotient (Sella R^{3N}/SE(3)), not a packed rotation");
  }
  if (mani == "se3" && dim != 12) {
    throw std::runtime_error(
        "Xtsci.manifold=se3 needs length 12; a 3N cluster is "
        "rigid_quotient (Sella R^{3N}/SE(3)), not an SE(3) prefix");
  }
  if ((mani == "rigid_quotient" || mani == "mw_rigid" || mani == "sella" ||
       mani == "eckart" || mani == "irc") &&
      (dim < 6 || dim % 3 != 0)) {
    throw std::runtime_error("Xtsci.manifold=" + mani +
                             " needs 3N Cartesians with N >= 2");
  }
  if (mani == "sphere") {
    xts_solver_set_manifold(m_solver, XTS_MANIFOLD_SPHERE);
  } else if (mani == "so3") {
    xts_solver_set_manifold(m_solver, XTS_MANIFOLD_SO3);
  } else if (mani == "stiefel") {
    xts_solver_set_manifold(m_solver, XTS_MANIFOLD_STIEFEL);
  } else if (mani == "se3") {
    xts_solver_set_manifold(m_solver, XTS_MANIFOLD_SE3);
  } else if (mani == "rigid_quotient" || mani == "sella") {
    xts_solver_set_manifold(m_solver, XTS_MANIFOLD_RIGID_QUOTIENT);
  } else if (mani == "mw_rigid" || mani == "eckart" || mani == "irc") {
    xts_solver_set_manifold(m_solver, XTS_MANIFOLD_MW_RIGID);
    const auto masses = m_objf->getMasses();
    if (masses.size() > 0) {
      xts_solver_set_masses(m_solver, masses.data(),
                            static_cast<size_t>(masses.size()));
    }
  } else {
    xts_solver_set_manifold(m_solver, XTS_MANIFOLD_EUCLIDEAN);
  }
  if (mani == "rigid_quotient" || mani == "mw_rigid" || mani == "sella" ||
      mani == "eckart" || mani == "irc") {
    xts_solver_set_periodic(m_solver, m_objf->getPeriodic() ? 1 : 0);
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

  if (m_x.size() == 0) {
    m_x = m_objf->getPositions();
  }
  if (m_x.size() != m_objf->degreesOfFreedom()) {
    throw std::runtime_error("xtsci objective position dimension mismatch");
  }
  auto *tensor =
      xts_tensor_borrow_cpu_f64(m_x.data(), static_cast<size_t>(m_x.size()));
  if (tensor == nullptr) {
    throw std::runtime_error("could not allocate xtsci objective tensor");
  }
  const auto &lbfgs = m_optConfig.opts.lbfgs;
  std::string precon = resolved_precon(m_optConfig.opts);
  const xts_method_t method = method_from_name(m_optConfig.opts.xtsci.method);
  if (method == XTS_DOGLEG && !is_host_precon(precon)) {
    precon = "pair";
  }
  XtsObjectiveContext context{
      m_objf.get(),    m_optConfig.potential, precon,      lbfgs.precon_A,
      lbfgs.precon_mu, lbfgs.precon_rcut,     &m_cached_x,
  };
  xts_report_t report{};
  const bool want_hess = method == XTS_NEWTON || method == XTS_RFO ||
                         method == XTS_DOGLEG || is_host_precon(precon);
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

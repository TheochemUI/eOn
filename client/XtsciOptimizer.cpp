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

#if XTS_ABI_VERSION_MINOR < 1
#error "xtsci-optimize ABI minor 1 or newer is required for Newton/RFO"
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
  return {static_cast<const double *>(dl.data) + dl.byte_offset / sizeof(double),
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

xts_status_t evaluate(void *user, const DLManagedTensorVersioned *x,
                      double *value_out) {
  try {
    auto *context = static_cast<XtsObjectiveContext *>(user);
    const auto positions = map_input(x);
    context->objective->setPositions(positions);
    *value_out = context->objective->getEnergy();
    return XTS_SUCCESS;
  } catch (...) {
    return XTS_INTERNAL_ERROR;
  }
}

xts_status_t gradient(void *user, const DLManagedTensorVersioned *x,
                      DLManagedTensorVersioned *gradient_out) {
  try {
    auto *context = static_cast<XtsObjectiveContext *>(user);
    const auto positions = map_input(x);
    context->objective->setPositions(positions);
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
  if (name == "newton") {
    return XTS_NEWTON;
  }
  if (name == "rfo") {
    return XTS_RFO;
  }
  return XTS_LBFGS;
}

} // namespace

namespace eonc {

int XtsciOptimizer::step(double a_maxMove) { return run(1, a_maxMove); }

int XtsciOptimizer::run(size_t a_maxIterations, double a_maxMove) {
  if (m_objf->degreesOfFreedom() <= 0 || a_maxIterations == 0) {
    return m_objf->isConverged() ? 1 : 0;
  }

  auto positions = m_objf->getPositions();
  if (positions.size() != m_objf->degreesOfFreedom()) {
    throw std::runtime_error("xtsci objective position dimension mismatch");
  }
  const auto stamp = xts_abi_stamp();
  if (xts_abi_compatible(&stamp) == 0) {
    throw std::runtime_error("incompatible xtsci-optimize ABI");
  }
  auto *tensor = xts_tensor_borrow_cpu_f64(positions.data(),
                                           static_cast<size_t>(positions.size()));
  if (tensor == nullptr) {
    throw std::runtime_error("could not allocate xtsci objective tensor");
  }
  const auto &lbfgs = m_optConfig.opts.lbfgs;
  XtsObjectiveContext context{
      m_objf.get(),
      m_optConfig.potential,
      lbfgs.precon,
      lbfgs.precon_A,
      lbfgs.precon_mu,
      lbfgs.precon_rcut,
  };
  xts_control_t control{
      a_maxIterations,
      m_optConfig.opts.converged_force,
      std::max(a_maxMove, 1.0e-12),
      static_cast<size_t>(std::max<long>(0, lbfgs.memory)),
      std::max(a_maxMove, 0.0),
  };
  xts_report_t report{};
  const xts_method_t method = method_from_name(m_optConfig.opts.xtsci_method);
  xts_status_t status;
  if (method == XTS_NEWTON || method == XTS_RFO) {
    status = xts_minimize_hess(evaluate, gradient, hessian, &context, tensor,
                               &control, method, &report);
  } else {
    status = xts_minimize(evaluate, gradient, &context, tensor, &control,
                          method, &report);
  }
  xts_tensor_free(tensor);
  if (status != XTS_SUCCESS) {
    throw std::runtime_error(xts_last_error());
  }
  m_objf->setPositions(positions);
  return m_objf->isConverged() ? 1 : 0;
}

} // namespace eonc

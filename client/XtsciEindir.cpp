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
#include "eon/XtsciEindir.h"

#include <format>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#ifdef WITH_EINDIR
#include <eindir-core.h>
#include <xts.h>
#endif

namespace eonc::xtsci_eindir {

namespace {

std::string &provenance_slot() {
  static std::string text;
  return text;
}

bool scalar_ok(long actual, long required) {
  return required == 0 || actual == required;
}

bool text_ok(const std::string &actual, const std::string &required) {
  return required.empty() || actual == required;
}

} // namespace

#ifdef WITH_EINDIR

struct State {
  ObjectiveFunction *objective{nullptr};
  Eigen::VectorXd *cached{nullptr};
  std::string schema;
  std::string producer;
  std::string length_unit;
  std::string energy_unit;
  std::vector<double> low;
  std::vector<double> high;
  eindir_objective_descriptor_t desc{};
  eindir_objective_t obj{};
};

void set_positions_if_changed(State *state, const Eigen::VectorXd &x) {
  auto *cached = state->cached;
  if (cached != nullptr && cached->size() == x.size() &&
      (*cached - x).isZero(0.0)) {
    return;
  }
  const auto cur = state->objective->getPositions();
  if (cur.size() == x.size() && (cur - x).isZero(0.0)) {
    if (cached != nullptr) {
      *cached = x;
    }
    return;
  }
  state->objective->setPositions(x);
  if (cached != nullptr) {
    *cached = x;
  }
}

Eigen::Map<const Eigen::VectorXd>
map_input(const DLManagedTensorVersioned *tensor) {
  const auto &dl = tensor->dl_tensor;
  if (dl.ndim != 1 || dl.dtype.code != kDLFloat || dl.dtype.bits != 64 ||
      dl.dtype.lanes != 1 || dl.device.device_type != kDLCPU ||
      dl.shape == nullptr || dl.data == nullptr) {
    throw std::runtime_error("eindir objective requires a CPU f64 vector");
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
    throw std::runtime_error("eindir gradient requires a CPU f64 vector");
  }
  return {static_cast<double *>(dl.data) + dl.byte_offset / sizeof(double),
          static_cast<Eigen::Index>(dl.shape[0])};
}

eindir_status_t eval_cb(void *user, const DLManagedTensorVersioned *x,
                        double *value_out) {
  try {
    auto *state = static_cast<State *>(user);
    set_positions_if_changed(state, map_input(x));
    *value_out = state->objective->getEnergy();
    return EINDIR_SUCCESS;
  } catch (...) {
    return EINDIR_INTERNAL_ERROR;
  }
}

eindir_status_t grad_cb(void *user, const DLManagedTensorVersioned *x,
                        DLManagedTensorVersioned *grad_out) {
  try {
    auto *state = static_cast<State *>(user);
    set_positions_if_changed(state, map_input(x));
    const auto gradient = state->objective->getGradient();
    auto output = map_output(grad_out);
    if (output.size() != gradient.size()) {
      return EINDIR_INVALID_PARAMETER;
    }
    output = gradient;
    return EINDIR_SUCCESS;
  } catch (...) {
    return EINDIR_INTERNAL_ERROR;
  }
}

#endif

View eon_requirement() {
  View view;
  view.schema_id = kSchemaId;
  view.length_unit = "angstrom";
  view.energy_unit = "eV";
  view.energy_sign = 1;
  view.gradient_sign = 1;
  view.operations = kOpEnergy | kOpForces;
  view.tensor_device_type = kDeviceCpu;
  view.tensor_dtype_code = kDtypeFloat;
  view.tensor_dtype_bits = 64;
  view.tensor_dtype_lanes = 1;
  view.tensor_layout = kLayoutContiguous;
  view.callback_lifetime = kLifetimeBorrowed;
  return view;
}

View rgpot_ev_angstrom() {
  View view = eon_requirement();
  view.producer_id = "rgpot";
  return view;
}

bool compatible(const View &actual, const View &required) {
  if (actual.schema_id.empty() || actual.schema_id != required.schema_id) {
    return false;
  }
  if (!text_ok(actual.producer_id, required.producer_id) ||
      !text_ok(actual.length_unit, required.length_unit) ||
      !text_ok(actual.energy_unit, required.energy_unit)) {
    return false;
  }
  if (!scalar_ok(actual.energy_sign, required.energy_sign) ||
      !scalar_ok(actual.gradient_sign, required.gradient_sign)) {
    return false;
  }
  if ((actual.operations & required.operations) != required.operations) {
    return false;
  }
  return scalar_ok(actual.tensor_device_type, required.tensor_device_type) &&
         scalar_ok(actual.tensor_dtype_code, required.tensor_dtype_code) &&
         scalar_ok(actual.tensor_dtype_bits, required.tensor_dtype_bits) &&
         scalar_ok(actual.tensor_dtype_lanes, required.tensor_dtype_lanes) &&
         scalar_ok(static_cast<long>(actual.tensor_layout),
                   static_cast<long>(required.tensor_layout)) &&
         scalar_ok(static_cast<long>(actual.callback_lifetime),
                   static_cast<long>(required.callback_lifetime));
}

bool abi_accepts_gradient(std::uint32_t major, std::uint32_t expected_major,
                          std::uint64_t features) {
  return major == expected_major && (features & kFeatureGradient) != 0;
}

const std::string &provenance() { return provenance_slot(); }

State *bind(ObjectiveFunction *objective, Eigen::VectorXd *cached) {
#ifndef WITH_EINDIR
  (void)objective;
  (void)cached;
  return nullptr;
#else
  if (objective == nullptr || objective->degreesOfFreedom() <= 0) {
    return nullptr;
  }
  auto state = std::make_unique<State>();
  state->objective = objective;
  state->cached = cached;
  const auto dim = static_cast<std::size_t>(objective->degreesOfFreedom());
  const double inf = std::numeric_limits<double>::infinity();
  state->low.assign(dim, -inf);
  state->high.assign(dim, inf);
  View actual = eon_requirement();
  actual.producer_id = "eon.objective";
  state->schema = actual.schema_id;
  state->producer = actual.producer_id;
  state->length_unit = actual.length_unit;
  state->energy_unit = actual.energy_unit;
  state->desc.schema_id = state->schema.c_str();
  state->desc.producer_id = state->producer.c_str();
  state->desc.length_unit = state->length_unit.c_str();
  state->desc.energy_unit = state->energy_unit.c_str();
  state->desc.energy_sign = actual.energy_sign;
  state->desc.gradient_sign = actual.gradient_sign;
  state->desc.operations = actual.operations;
  state->desc.tensor_device_type = actual.tensor_device_type;
  state->desc.tensor_dtype_code = actual.tensor_dtype_code;
  state->desc.tensor_dtype_bits =
      static_cast<std::uint8_t>(actual.tensor_dtype_bits);
  state->desc.tensor_dtype_lanes =
      static_cast<std::uint8_t>(actual.tensor_dtype_lanes);
  state->desc.tensor_layout = actual.tensor_layout;
  state->desc.callback_lifetime = actual.callback_lifetime;
  state->obj.dim = dim;
  state->obj.low = state->low.data();
  state->obj.high = state->high.data();
  state->obj.eval_fn = eval_cb;
  state->obj.grad_fn = grad_cb;
  state->obj.user_data = state.get();
  state->obj.free_fn = nullptr;
  state->obj.descriptor = &state->desc;
  if (!compatible(actual, eon_requirement())) {
    throw std::runtime_error("eOn eindir objective descriptor was rejected");
  }
  const auto stamp = eindir_core_abi_stamp();
  if (eindir_core_abi_compatible(&stamp) == 0 ||
      !abi_accepts_gradient(stamp.abi_major, 1, stamp.features)) {
    throw std::runtime_error("incompatible eindir ABI stamp");
  }
  if (eindir_objective_has_grad(&state->obj) == 0) {
    throw std::runtime_error("eindir objective has no gradient");
  }
  eindir_objective_descriptor_t required_desc = state->desc;
  required_desc.producer_id = "";
  if (eindir_objective_descriptor_compatible(&state->desc, &required_desc) ==
      0) {
    const char *err = eindir_last_error();
    throw std::runtime_error(err != nullptr ? err
                                            : "eindir descriptor mismatch");
  }
  const auto xts = xts_abi_stamp();
  provenance_slot() =
      std::format("xts-{}.{}.{}+eindir-{}.{}-layout-{}", xts.abi_major,
                  xts.abi_minor, xts.layout_revision, stamp.abi_major,
                  stamp.abi_minor, stamp.objective_layout);
  return state.release();
#endif
}

void release(State *state) {
#ifdef WITH_EINDIR
  delete state;
#else
  (void)state;
#endif
}

int eval_grad(State *state, const DLManagedTensorVersioned *x, double *value,
              DLManagedTensorVersioned *gradient) {
#ifndef WITH_EINDIR
  (void)state;
  (void)x;
  (void)value;
  (void)gradient;
  return 2;
#else
  if (state == nullptr) {
    return 1;
  }
  if (eindir_objective_eval(&state->obj, x, value) != EINDIR_SUCCESS ||
      eindir_objective_grad(&state->obj, x, gradient) != EINDIR_SUCCESS) {
    return 2;
  }
  return 0;
#endif
}

int minimize(State *state, double *x, std::size_t n, std::size_t maxiter,
             double gtol, double istep, std::size_t memory, int method,
             double *value_out) {
#ifndef WITH_EINDIR
  (void)state;
  (void)x;
  (void)n;
  (void)maxiter;
  (void)gtol;
  (void)istep;
  (void)memory;
  (void)method;
  (void)value_out;
  return 1;
#else
  if (state == nullptr || x == nullptr) {
    return 1;
  }
  auto *tensor = xts_tensor_borrow_cpu_f64(x, n);
  if (tensor == nullptr) {
    return 2;
  }
  const auto stamp = eindir_core_abi_stamp();
  xts_control_t control{maxiter, gtol, istep, memory, 0.0};
  xts_report_t report{};
  const auto status =
      xts_minimize_eindir(&state->obj, &stamp, tensor, &control,
                          static_cast<xts_method_t>(method), &report);
  xts_tensor_free(tensor);
  if (value_out != nullptr) {
    *value_out = report.value;
  }
  return static_cast<int>(status);
#endif
}

} // namespace eonc::xtsci_eindir

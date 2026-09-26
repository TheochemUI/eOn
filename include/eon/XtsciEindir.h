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
#pragma once

#include "eon/ObjectiveFunction.h"

#include <cstddef>
#include <cstdint>
#include <string>

struct DLManagedTensorVersioned;

namespace eonc::xtsci_eindir {

/// Semantic objective contract shared with eindir_objective_descriptor_t.
/// Empty strings and zero scalars are wildcards, matching that predicate.
struct View {
  std::string schema_id;
  std::string producer_id;
  std::string length_unit;
  std::string energy_unit;
  int energy_sign{0};
  int gradient_sign{0};
  std::uint64_t operations{0};
  int tensor_device_type{0};
  int tensor_dtype_code{0};
  int tensor_dtype_bits{0};
  int tensor_dtype_lanes{0};
  std::uint32_t tensor_layout{0};
  std::uint32_t callback_lifetime{0};
};

inline constexpr const char *kSchemaId = "eindir.objective-descriptor.v1";
inline constexpr std::uint64_t kOpEnergy = 1ull << 0;
inline constexpr std::uint64_t kOpForces = 1ull << 1;
inline constexpr int kDeviceCpu = 1;
inline constexpr int kDtypeFloat = 2;
inline constexpr std::uint32_t kLayoutContiguous = 1;
inline constexpr std::uint32_t kLifetimeBorrowed = 1;
inline constexpr std::uint64_t kFeatureGradient = 1ull << 0;

/// What an eOn minimization may consume: eV, angstrom, analytic forces.
View eon_requirement();

/// rgpot producer advertising the same units and tensor contract.
View rgpot_ev_angstrom();

/// Nonzero when actual satisfies required. Schema id is never a wildcard.
bool compatible(const View &actual, const View &required);

/// Major must match and the gradient feature bit must be set.
bool abi_accepts_gradient(std::uint32_t major, std::uint32_t expected_major,
                          std::uint64_t features);

/// Stamp text from the last successful eindir bind. Empty when unused.
const std::string &provenance();

struct State;

/// Borrow an ObjectiveFunction as an eindir objective. Null when eindir
/// is not linked. Does not take ownership of the objective or the cache.
State *bind(ObjectiveFunction *objective, Eigen::VectorXd *cached);

void release(State *state);

/// One fused energy and gradient through the borrowed eindir handle.
int eval_grad(State *state, const DLManagedTensorVersioned *x, double *value,
              DLManagedTensorVersioned *gradient);

/// xts_minimize_eindir on a caller-owned buffer. The caller keeps State.
int minimize(State *state, double *x, std::size_t n, std::size_t maxiter,
             double gtol, double istep, std::size_t memory, int method,
             double *value_out);

} // namespace eonc::xtsci_eindir

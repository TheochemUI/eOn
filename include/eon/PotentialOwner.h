// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <limits>

namespace eonc::potential_owner {
// Ring beads use their nonnegative index. Auxiliary geometry and derivative
// stencils each reuse one scratch slot, bounding retained electronic states.
inline constexpr long auxiliary = std::numeric_limits<long>::max() - 1;
inline constexpr long derivative = std::numeric_limits<long>::max();
inline constexpr long reactant = std::numeric_limits<long>::max() - 2;
inline constexpr long saddle = std::numeric_limits<long>::max() - 3;
} // namespace eonc::potential_owner

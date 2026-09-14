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

#include "Potential.h"
#include <concepts>
#include <type_traits>

namespace eonc {

/// Duck-typed thread-safety query. Virtual `Potential::isThreadSafe` stays
/// the override point; callers should prefer these nonmembers so a later
/// trait/concept cut does not touch every NEB/dimer TU.
template <class P>
concept ThreadSafeQueryable = requires(const P &p) {
  { p.isThreadSafe() } -> std::convertible_to<bool>;
};

template <class P>
  requires ThreadSafeQueryable<P>
[[nodiscard]] bool potIsThreadSafe(const P &p) noexcept {
  return p.isThreadSafe();
}

template <class P>
  requires requires(const P &p) {
    { p.isSharedInstanceThreadSafe() } -> std::convertible_to<bool>;
  }
[[nodiscard]] bool potAllowsSharedInstance(const P &p) noexcept {
  return p.isSharedInstanceThreadSafe();
}

} // namespace eonc

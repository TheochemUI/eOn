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

#include <pybind11/embed.h>

namespace eonc {

// Lazy singleton for the Python interpreter. The interpreter is only started on
// the first call and lives until program exit. This avoids paying the cost of
// Py_Initialize for runs that never use a Python-based potential.
inline void ensure_interpreter() {
  static pybind11::scoped_interpreter guard{};
  (void)guard;
}

} // namespace eonc

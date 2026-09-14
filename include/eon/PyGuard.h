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

// Lazy, one-shot Python interpreter init.  We intentionally never call
// Py_Finalize so that pybind11 py::object members in potentials can release
// their refcounts during normal destruction without racing against interpreter
// teardown.  The OS reclaims everything on process exit.
inline void ensure_interpreter() {
  if (!Py_IsInitialized()) {
    pybind11::initialize_interpreter();
  }
}

} // namespace eonc

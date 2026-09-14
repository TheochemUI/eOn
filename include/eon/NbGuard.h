/*
** Interpreter init without pybind11.
**
** Polarity:
**   * Preferred: Python owns the process (eon-server + pyeonclient).
*Interpreter
**     already exists; ensure_interpreter() is a no-op.
**   * Fallback: standalone eonclient needs a pot that calls into Python (ASE).
**     Use Py_InitializeEx — do NOT depend on pybind11::embed.
**
** Object/array marshalling for pots should use nanobind (nb::object,
*nb::ndarray)
** once the interpreter is up, matching pyeonclient.
*/
#pragma once

#include <Python.h>

namespace eonc {

inline void ensure_interpreter() {
  if (!Py_IsInitialized()) {
    // 0 = skip signal handlers (same idea as pybind11 embed defaults)
    Py_InitializeEx(0);
  }
}

} // namespace eonc

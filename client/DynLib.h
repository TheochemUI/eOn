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

/// Cross-platform dynamic library loading abstraction.
/// Used by LammpsLoader, XtbLoader, and any other runtime-loaded potentials.

#include <string>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace eonc::dynlib {

#ifdef _WIN32
using Handle = HMODULE;

inline Handle open(const char *name) noexcept { return LoadLibraryA(name); }

inline void *sym(Handle h, const char *name) noexcept {
  return reinterpret_cast<void *>(GetProcAddress(h, name));
}

inline void close(Handle h) noexcept {
  if (h)
    FreeLibrary(h);
}

inline std::string error() {
  DWORD err = GetLastError();
  if (err == 0)
    return {};
  LPSTR buf = nullptr;
  FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                     FORMAT_MESSAGE_IGNORE_INSERTS,
                 nullptr, err, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                 reinterpret_cast<LPSTR>(&buf), 0, nullptr);
  std::string msg(buf ? buf : "unknown error");
  if (buf)
    LocalFree(buf);
  return msg;
}

#else // POSIX (Linux, macOS)
using Handle = void *;

inline Handle open(const char *name) noexcept {
  return dlopen(name, RTLD_NOW | RTLD_LOCAL);
}

inline void *sym(Handle h, const char *name) noexcept { return dlsym(h, name); }

inline void close(Handle h) noexcept {
  if (h)
    dlclose(h);
}

inline std::string error() {
  const char *msg = dlerror();
  return msg ? std::string(msg) : std::string{};
}
#endif

/// Try a list of library names in order, return first successful handle.
inline Handle openFirst(const char *const names[]) noexcept {
  for (const char *const *name = names; *name; ++name) {
    Handle h = open(*name);
    if (h)
      return h;
  }
  return {};
}

/// Load a symbol and cast to a function pointer type.
template <typename Fn> Fn loadSym(Handle h, const char *name) noexcept {
  return reinterpret_cast<Fn>(sym(h, name));
}

} // namespace eonc::dynlib

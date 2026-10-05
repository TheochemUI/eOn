#pragma once
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
/*
** Thin loader for any rgpot engine plugin behind the generic C ABI
** (engine_c_abi.h): one symbol set, backend selection by library name,
** configuration as a Cap'n Proto message on the shared wire. New
** backends need no new loader code.
*/
#include "engine_c_abi.h"
#include <cstddef>
#include <string>
#include <vector>

struct GenericEngineOptions {
  /// Library basename to search for, e.g. "libuma_engine.so".
  std::string library;
  /// Optional env var naming an explicit .so path (checked first).
  std::string env_var;
  /// Optional explicit .so path from config.
  std::string engine_path;
  /// Cap'n Proto flat-array message (the engine's params struct on the
  /// shared Potentials.capnp wire) handed to rgpot_engine_create.
  std::vector<unsigned char> config;
  /// Backend tag for error messages ("uma", ...).
  std::string tag;
};

class GenericEngineLoader {
public:
  explicit GenericEngineLoader(const GenericEngineOptions &opt);
  ~GenericEngineLoader();
  GenericEngineLoader(const GenericEngineLoader &) = delete;
  GenericEngineLoader &operator=(const GenericEngineLoader &) = delete;

  [[nodiscard]] bool available() const noexcept { return m_pot != nullptr; }
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) const;

private:
  void *m_lib{nullptr};
  RgpotEnginePot *m_pot{nullptr};
  std::string m_tag;
  using create_fn = RgpotEnginePot *(*)(const void *, size_t, char *, size_t);
  using destroy_fn = void (*)(RgpotEnginePot *);
  using force_fn = int (*)(RgpotEnginePot *, long, const double *, const int *,
                           double *, double *, double *, const double *,
                           rgpot_engine_coord_transform, void *);
  create_fn m_create{nullptr};
  destroy_fn m_destroy{nullptr};
  force_fn m_force{nullptr};
};

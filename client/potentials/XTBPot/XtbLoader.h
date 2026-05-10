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

/// Runtime loader for the libxtb C library.
///
/// Mirrors the LammpsLoader / ARTnResource pattern: dlopen at first
/// access, cache function pointers, throw via require_loaded() when
/// the user actually requests an XTB potential. This lets a single
/// eOn binary support XTB iff libxtb is on the library search path,
/// without a build-time link against tblite + xtb's transitive
/// Fortran dependency graph.
///
/// The xtb opaque struct types (xtb_TEnvironment, xtb_TMolecule,
/// xtb_TCalculator, xtb_TResults) survive across the ABI as plain
/// pointer-to-struct. We re-typedef them as void* here so the loader
/// doesn't need xtb.h.

#include "../../DynLib.h"

#include <stdexcept>

namespace eonc {

class XtbLoader {
public:
  // Opaque pointer types. xtb.h defines them as
  // `typedef struct _xtb_TX* xtb_TX;`; for the loader they're plain
  // void* so we keep compile coupling to xtb.h zero. XTBPot still
  // includes xtb.h and casts at the call sites.
  using env_t = void *;
  using mol_t = void *;
  using calc_t = void *;
  using res_t = void *;

  // Function pointer types -- one per xtb.h entry point used in
  // XTBPot::force / XTBPot ctor / XTBPot dtor.
  using new_environment_fn = env_t (*)(void);
  using del_environment_fn = void (*)(env_t *);
  using check_environment_fn = int (*)(env_t);
  using get_error_fn = void (*)(env_t, char *, const int *);
  using release_output_fn = void (*)(env_t);
  using set_verbosity_fn = void (*)(env_t, int);

  using new_molecule_fn = mol_t (*)(env_t, const int *, const int *,
                                    const double *, const double *, const int *,
                                    const double *, const bool *);
  using del_molecule_fn = void (*)(mol_t *);
  using update_molecule_fn = void (*)(env_t, mol_t, const double *,
                                      const double *);

  using new_calculator_fn = calc_t (*)(void);
  using del_calculator_fn = void (*)(calc_t *);
  using load_gfnff_fn = void (*)(env_t, mol_t, calc_t, char *);
  using load_gfn0_fn = void (*)(env_t, mol_t, calc_t, char *);
  using load_gfn1_fn = void (*)(env_t, mol_t, calc_t, char *);
  using load_gfn2_fn = void (*)(env_t, mol_t, calc_t, char *);
  using set_accuracy_fn = void (*)(env_t, calc_t, double);
  using set_max_iter_fn = void (*)(env_t, calc_t, int);
  using set_electronic_temp_fn = void (*)(env_t, calc_t, double);

  using singlepoint_fn = void (*)(env_t, mol_t, calc_t, res_t);

  using new_results_fn = res_t (*)(void);
  using del_results_fn = void (*)(res_t *);
  using get_energy_fn = void (*)(env_t, res_t, double *);
  using get_gradient_fn = void (*)(env_t, res_t, double *);

  /// Thread-safe singleton accessor (Meyer's pattern).
  static XtbLoader &instance();

  // Loaded function pointers; null when libxtb is missing or the
  // installed copy lacks a symbol.
  new_environment_fn new_environment{nullptr};
  del_environment_fn del_environment{nullptr};
  check_environment_fn check_environment{nullptr};
  get_error_fn get_error{nullptr};
  release_output_fn release_output{nullptr};
  set_verbosity_fn set_verbosity{nullptr};
  new_molecule_fn new_molecule{nullptr};
  del_molecule_fn del_molecule{nullptr};
  update_molecule_fn update_molecule{nullptr};
  new_calculator_fn new_calculator{nullptr};
  del_calculator_fn del_calculator{nullptr};
  load_gfnff_fn load_gfnff{nullptr};
  load_gfn0_fn load_gfn0{nullptr};
  load_gfn1_fn load_gfn1{nullptr};
  load_gfn2_fn load_gfn2{nullptr};
  set_accuracy_fn set_accuracy{nullptr};
  set_max_iter_fn set_max_iter{nullptr};
  set_electronic_temp_fn set_electronic_temp{nullptr};
  singlepoint_fn singlepoint{nullptr};
  new_results_fn new_results{nullptr};
  del_results_fn del_results{nullptr};
  get_energy_fn get_energy{nullptr};
  get_gradient_fn get_gradient{nullptr};

  [[nodiscard]] bool is_loaded() const noexcept { return m_loaded; }

  /// Throws std::runtime_error if libxtb is not available.
  void require_loaded() const;

  XtbLoader(const XtbLoader &) = delete;
  XtbLoader &operator=(const XtbLoader &) = delete;

private:
  XtbLoader();
  ~XtbLoader();

  bool m_loaded{false};
  dynlib::Handle m_handle{};
};

inline XtbLoader &get_xtb_loader() { return XtbLoader::instance(); }

} // namespace eonc

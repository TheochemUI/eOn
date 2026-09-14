/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Pure C ABI for libmetatomic_pot.so — loadable via dlopen like Fortran pots /
** liblammps. Host eOn never includes torch headers.
*/
#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifdef EON_MTA_BUILD
#define EON_MTA_API __declspec(dllexport)
#else
#define EON_MTA_API __declspec(dllimport)
#endif
#else
#define EON_MTA_API __attribute__((visibility("default")))
#endif

/** Opaque pot handle. */
typedef struct EonMtaPot EonMtaPot;

/**
 * Flat config for create (mirrors Parameters::metatomic_options_t).
 * All string pointers may be NULL (treated as empty / defaults).
 */
typedef struct EonMtaConfig {
  const char *model_path;
  const char *device;
  const char *length_unit;
  const char *extensions_directory;
  int check_consistency;
  double uncertainty_threshold;
  const char *energy_output;
  const char *energy_uncertainty_output;
  const char *force_output;
  int non_conservative;
  int random_rotation;
  long n_symmetry_rotations;
  int deterministic;
  int deterministic_strict;
  const char *variant_base;
  const char *variant_energy;
  const char *variant_energy_uncertainty;
  const char *variant_force;
} EonMtaConfig;

/**
 * Create pot. On failure returns NULL and writes a message into errbuf
 * (if non-NULL and errlen > 0).
 */
EON_MTA_API EonMtaPot *eon_mta_pot_create(const EonMtaConfig *cfg, char *errbuf,
                                          size_t errlen);

EON_MTA_API void eon_mta_pot_destroy(EonMtaPot *pot);

/**
 * Evaluate energy + forces. Returns 0 on success, nonzero on failure.
 * variance may be NULL.
 */
EON_MTA_API int eon_mta_pot_force(EonMtaPot *pot, long nAtoms,
                                  const double *positions, const int *atomicNrs,
                                  double *forces, double *energy,
                                  double *variance, const double *box);

/** ABI version for loader checks. */
EON_MTA_API int eon_mta_abi_version(void);

#define EON_MTA_ABI_VERSION 1

#ifdef __cplusplus
}
#endif

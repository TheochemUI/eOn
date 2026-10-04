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
/**
 * @file eon_relax_engine.h
 * @brief Tier-2 C ABI: eOn as a dlopen relaxation engine.
 *
 * Reverse of include/eon/potentials/Rgpot/engine_c_abi.h: the host
 * supplies the surface (energy/forces/variance for a band), eOn owns
 * NudgedElasticBand / MinModeSaddleSearch. Configuration is a Cap'n
 * Proto flat-array whose root is RelaxEngineParams in
 * schema/eon_relax_engine.capnp (NebParams / SaddleParams arms). NULL config
 * selects Parameters NSDMI defaults and kind=NEB.
 *
 * Units: positions Angstrom, energy eV, forces eV/Angstrom.
 *
 * Plugin Matters use the Matter default (PBC on for a non-isolated
 * pot). fill and endpoint dirty write positions with PBC off so the
 * caller frame is kept (a 10 A box must not map x=-0.3 to 9.7).
 *
 * Layout discipline (dlpack.h is the model):
 *
 * - Every wire struct opens with an eon_relax_version_t and a uint64
 *   flags word. The writer stamps the version it was built against;
 *   the reader refuses a major it does not know
 *   (EON_RELAX_ABI_MISMATCH). Growth is append-only behind the
 *   version head, gated on minor.
 * - Fixed-width integer types only. Never long in a wire struct.
 * - The surface callback takes ONE request struct, so the callback
 *   contract can grow without a new function symbol.
 * - Return values are three-valued: 0 success (read the outcome),
 *   positive recoverable, negative a named fail. SURFACE_FATAL makes
 *   the instance unusable. Every other named fail leaves the handle
 *   reusable.
 * - Closed enums. Unknown kind is EON_RELAX_UNKNOWN_KIND; other
 *   unknown enumerants are EON_RELAX_INVALID_PARAMETER. Never a
 *   silent fallback.
 * - Optional symbols are discovered by dlsym; absence means the
 *   operation is unsupported.
 */
#pragma once

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifdef EON_RELAX_ENGINE_BUILD
#define EON_RELAX_API __declspec(dllexport)
#else
#define EON_RELAX_API __declspec(dllimport)
#endif
#else
#define EON_RELAX_API __attribute__((visibility("default")))
#endif

typedef struct EonRelaxEngine EonRelaxEngine;

#define EON_RELAX_ABI_MAJOR 2u
#define EON_RELAX_ABI_MINOR 1u

/** Version head carried by every wire struct (dlpack DLPackVersion). */
typedef struct {
  uint32_t major;
  uint32_t minor;
} eon_relax_version_t;

#define EON_RELAX_VERSION_INIT                                                 \
  { EON_RELAX_ABI_MAJOR, EON_RELAX_ABI_MINOR }

typedef enum {
  EON_RELAX_KIND_INVALID = 0,
  EON_RELAX_KIND_NEB = 1,
  EON_RELAX_KIND_SADDLE = 2
} eon_relax_kind_t;

#define EON_RELAX_KIND_IS_KNOWN(k)                                             \
  ((k) == EON_RELAX_KIND_NEB || (k) == EON_RELAX_KIND_SADDLE)

/** Call rc: 0 success (read outcome), positive recoverable, negative a
 *  named fail. SURFACE_FATAL makes the instance unusable. Every other
 *  named fail leaves the handle reusable. */
typedef enum {
  EON_RELAX_OK = 0,
  EON_RELAX_NULL_ENGINE = -1,
  EON_RELAX_NULL_BAND = -2,
  EON_RELAX_NULL_SURFACE = -3,
  EON_RELAX_NULL_OUTCOME = -4,
  EON_RELAX_NULL_POSITIONS = -5,
  EON_RELAX_NULL_ATOMIC_NRS = -6,
  EON_RELAX_NULL_BOXES = -7,
  EON_RELAX_NULL_MODE = -8,
  EON_RELAX_BAND_TOO_SHORT = -9,
  EON_RELAX_NATOMS = -10,
  EON_RELAX_UNKNOWN_KIND = -11,
  EON_RELAX_UNKNOWN_STATUS = -12,
  EON_RELAX_CAPNP_ROOT = -13,
  EON_RELAX_SURFACE_FATAL = -15,
  EON_RELAX_SADDLE_NIMAGES = -17,
  EON_RELAX_BAND_SIZE = -18,
  EON_RELAX_ALLOC = -19,
  EON_RELAX_UNAVAILABLE = -21,
  EON_RELAX_INVALID_PARAMETER = -22,
  EON_RELAX_ABI_MISMATCH = -23
} eon_relax_rc_t;

/** NudgedElasticBand::NEBStatus integers. Written to outcome.status for NEB. */
typedef enum {
  EON_RELAX_NEB_GOOD = 0,
  EON_RELAX_NEB_INIT = 1,
  EON_RELAX_NEB_BAD_MAX_ITERATIONS = 2,
  EON_RELAX_NEB_RUNNING = 3,
  EON_RELAX_NEB_MAX_UNCERTAINTY = 4
} eon_relax_neb_status_t;

#define EON_RELAX_IMAGE_NONE INT64_C(-1)

/**
 * Caller-owned band. The caller stamps version (EON_RELAX_VERSION_INIT)
 * and zeroes flags. positions is 3*n_atoms*n_images (image stride
 * 3*n_atoms). A successful NEB run or step writes interior images only;
 * endpoints stay in the caller frame. atomic_nrs is n_atoms (shared
 * composition), boxes is 9*n_images in Matter::getCell().data()
 * order (eOn Matrix3d is RowMajor). is_fixed may be NULL (all free) or
 * n_atoms (1 = fixed atom). mode may be NULL; required for
 * kind=saddle (3*n_atoms). image_ids may be NULL; dest run/step
 * emits dest path[] slots on the surface request (pointer identity)
 * and does not read this field. masses may be NULL (each atom mass
 * 1.0) or n_atoms (amu); readers honor masses only when
 * version.minor >= 1. For kind=NEB, n_images must equal
 * Parameters.neb_options.image_count + 2 (NSDMI 7). positions layout
 * matches Matter AtomMatrix (row-major xyz per atom).
 */
typedef struct {
  eon_relax_version_t version;
  uint64_t flags;
  int64_t n_images;
  int64_t n_atoms;
  double *positions;
  const int32_t *atomic_nrs;
  const double *boxes;
  const int32_t *is_fixed;
  const double *mode;
  const int64_t *image_ids;
  const double *masses;
} eon_relax_band_t;

/** Engine-written result. The engine stamps version and flags. */
typedef struct {
  eon_relax_version_t version;
  uint64_t flags;
  int32_t kind;
  int32_t status;
  int64_t iterations;
  int64_t climbing_image;
  double max_force;
  uint64_t version_hash;
  uint64_t surface_epoch;
} eon_relax_outcome_t;

/**
 * Host-surface request: one struct, engine-stamped version, so the
 * callback contract grows behind the version head instead of by new
 * symbols. n_images systems of one composition. positions/forces are
 * 3*n_atoms*n_images, boxes 9*n_images, energies n_images, variances
 * n_images (NULL when the engine does not consume variance),
 * image_ids n_images (NULL when unmapped). epoch_out may be NULL; a
 * written value keys Matter caches on the new surface generation.
 */
typedef struct {
  eon_relax_version_t version;
  uint64_t flags;
  int64_t n_images;
  int64_t n_atoms;
  const double *positions;
  const int32_t *atomic_nrs;
  const double *boxes;
  const int64_t *image_ids;
  double *energies;
  double *forces;
  double *variances;
  uint64_t *epoch_out;
} eon_relax_surface_request_t;

/** Return 0 on success, positive recoverable, negative fatal. */
typedef int (*eon_relax_surface_fn)(void *user,
                                    eon_relax_surface_request_t *req);

/** (major << 16) | minor. */
EON_RELAX_API int eon_relax_abi_version(void);
EON_RELAX_API int eon_relax_available(void);
EON_RELAX_API int eon_relax_abi_stamp(eon_relax_version_t *out);
/** Immortal dest identity string. Never NULL. */
EON_RELAX_API const char *eon_relax_version_hash_str(void);
/** 64-bit FNV-1a hash of eon_relax_version_hash_str(). Copied to
 *  outcome.version_hash. */
EON_RELAX_API uint64_t eon_relax_version_hash(void);

/**
 * Create. NULL config && config_len==0 selects Parameters NSDMI
 * defaults and kind=NEB. A non-empty buffer is a Cap'n Proto
 * flat-array whose root is RelaxEngineParams. kind tokens are
 * "neb" and "saddle"; any other token fails closed
 * (EON_RELAX_UNKNOWN_KIND). An unparseable buffer is
 * EON_RELAX_CAPNP_ROOT. config!=NULL && config_len==0 is refused.
 */
EON_RELAX_API EonRelaxEngine *eon_relax_create(const void *config,
                                               size_t config_len, char *errbuf,
                                               size_t errlen);

EON_RELAX_API int eon_relax_set_surface_epoch(EonRelaxEngine *eng,
                                              uint64_t epoch);
/** Last named create/run/step/reset/set_epoch rc on this handle. 0 if none. */
EON_RELAX_API int eon_relax_last_error(const EonRelaxEngine *eng);

/**
 * Run the bound engine to completion. Writes updated interior
 * coordinates into band->positions (caller-owned); endpoints stay in
 * the caller frame. Endpoint minimization stays with the host
 * (NudgedElasticBand::compute does not minimize endpoints).
 * surface may not be NULL. Kind is bound at create (NULL config =>
 * NEB). Returns eon_relax_rc_t; on 0, out->status is the dest NEB
 * or saddle table.
 */
EON_RELAX_API int eon_relax_run(EonRelaxEngine *eng, eon_relax_band_t *band,
                                eon_relax_surface_fn surface, void *user,
                                eon_relax_outcome_t *out);

/**
 * One band optimizer step (kind=NEB only). The first call on a fresh
 * or reset engine initializes the band state from the caller's
 * positions and captures the baseline force; every call syncs
 * band->positions in (the host may have moved images between steps),
 * assembles NEB forces on the host surface, takes ONE optimizer step,
 * and writes the stepped INTERIOR positions back (endpoints are the
 * caller's and are never rewritten, so the caller's frame survives
 * eOn's own cell wrapping). Climbing-image activation
 * follows the engine's own trigger rule; saddle POLICY (MMF bursts,
 * resampling, acquisition) stays with the host. outcome.status is
 * RUNNING until convergenceForce() reaches the configured tolerance
 * (GOOD), host variance exceeds the GP-surrogate threshold
 * (MAX_UNCERTAINTY), or neb.max_iterations is hit
 * (BAD_MAX_ITERATIONS). surface and user must not change across
 * steps without an eon_relax_reset.
 */
EON_RELAX_API int eon_relax_step(EonRelaxEngine *eng, eon_relax_band_t *band,
                                 eon_relax_surface_fn surface, void *user,
                                 eon_relax_outcome_t *out);

/**
 * Drop the stepper's band state and optimizer history; the next
 * eon_relax_step reinitializes from the caller's band. This is the
 * host's model-update hook: quasi-Newton memory taken on one surface
 * epoch must not survive onto the next.
 */
EON_RELAX_API int eon_relax_reset(EonRelaxEngine *eng);

EON_RELAX_API void eon_relax_destroy(EonRelaxEngine *eng);

/** NULL if kind or status is unknown. */
EON_RELAX_API const char *eon_relax_status_name(eon_relax_kind_t kind,
                                                int status);

#ifdef __cplusplus
}
#endif

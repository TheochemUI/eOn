/**
 * @file metatomic_c_abi.h
 * @brief C ABI between Metatomic frontend and libmetatomic_engine.so.
 *
 * Units: positions Angstrom, energy eV, forces eV/Angstrom (eOn/rgpot host).
 * Mirrors the NWChem engine split: thin consumer, optional heavy torch engine.
 */
#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifdef RGPOT_MTA_ENGINE_BUILD
#define RGPOT_MTA_API __declspec(dllexport)
#else
#define RGPOT_MTA_API __declspec(dllimport)
#endif
#else
#define RGPOT_MTA_API __attribute__((visibility("default")))
#endif

typedef struct RgpotMtaPot RgpotMtaPot;

/** Flat config (subset of MetatomicConfig). NULL strings => empty/default. */
typedef struct RgpotMtaConfig {
  const char *model_path;
  const char *device;
  const char *length_unit;
  const char *extensions_directory;
  int check_consistency;
  double uncertainty_threshold;
  const char *dtype_override;
  int random_rotation;
  long n_symmetry_rotations;
  int so3_probe_scatter;
  int torch_determinism_strict; /**< 0=Fast, 1=Strict */
} RgpotMtaConfig;

#define RGPOT_MTA_ABI_VERSION 1

RGPOT_MTA_API int rgpot_mta_abi_version(void);

RGPOT_MTA_API RgpotMtaPot *rgpot_mta_create(const RgpotMtaConfig *cfg,
                                            char *errbuf, size_t errlen);

RGPOT_MTA_API void rgpot_mta_destroy(RgpotMtaPot *pot);

/**
 * Energy + forces. Returns 0 on success.
 * variance may be NULL.
 */
RGPOT_MTA_API int rgpot_mta_force(RgpotMtaPot *pot, long nAtoms,
                                  const double *positions, const int *atomicNrs,
                                  double *forces, double *energy,
                                  double *variance, const double *box);

RGPOT_MTA_API int rgpot_mta_available(void);

#ifdef __cplusplus
}
#endif

/**
 * @file xtb_c_abi.h
 * @brief C ABI between XTBDlopen frontend and libxtb_engine.so.
 *
 * Engine wraps libxtb's ISO_C_BINDING API (`xtb.h`) behind a small stable
 * plugin surface (create / force / destroy). Host units: positions Å, energy
 * eV, forces eV/Å. Mirrors the metatomic engine split: thin consumer,
 * optional heavy plugin.
 */
#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifdef RGPOT_XTB_ENGINE_BUILD
#define RGPOT_XTB_API __declspec(dllexport)
#else
#define RGPOT_XTB_API __declspec(dllimport)
#endif
#else
#define RGPOT_XTB_API __attribute__((visibility("default")))
#endif

typedef struct RgpotXtbPot RgpotXtbPot;

/** GFN method codes (stable ABI). */
enum {
  RGPOT_XTB_METHOD_GFNFF = 0,
  RGPOT_XTB_METHOD_GFN0 = 1,
  RGPOT_XTB_METHOD_GFN1 = 2,
  RGPOT_XTB_METHOD_GFN2 = 3
};

typedef struct RgpotXtbConfig {
  int method; /**< RGPOT_XTB_METHOD_* */
  double accuracy;
  double electronic_temperature; /**< Kelvin */
  int max_iterations;
  double charge;
  int uhf;
} RgpotXtbConfig;

#define RGPOT_XTB_ABI_VERSION 1

RGPOT_XTB_API int rgpot_xtb_abi_version(void);

RGPOT_XTB_API RgpotXtbPot *rgpot_xtb_create(const RgpotXtbConfig *cfg,
                                            char *errbuf, size_t errlen);

RGPOT_XTB_API void rgpot_xtb_destroy(RgpotXtbPot *pot);

/**
 * Energy + forces. Returns 0 on success.
 * variance may be NULL.
 */
RGPOT_XTB_API int rgpot_xtb_force(RgpotXtbPot *pot, long nAtoms,
                                  const double *positions, const int *atomicNrs,
                                  double *forces, double *energy,
                                  double *variance, const double *box);

RGPOT_XTB_API int rgpot_xtb_available(void);

#ifdef __cplusplus
}
#endif

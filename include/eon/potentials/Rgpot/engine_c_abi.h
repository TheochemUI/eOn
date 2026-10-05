/**
 * @file engine_c_abi.h
 * @brief Generic C ABI for rgpot engine plugins.
 *
 * One symbol set for every heavy backend: an engine shared library
 * (libuma_engine.so, libmetatomic_engine.so, ...) exports these
 * functions, backend selection is which library the host dlopens, and
 * configuration travels as a Cap'n Proto flat-array message on the
 * shared Potentials.capnp wire (the engine names its root params
 * struct, e.g. UmaParams). Adding a backend therefore needs no new
 * host code and no new ABI header.
 *
 * Units: positions Angstrom, energy eV, forces eV/Angstrom.
 *
 * The optional coordinate transform lets the host inject a
 * pre-evaluation map over (positions, box) — for example the
 * molecular-box convention that recenters a gas-phase molecule in a
 * large fixed cell so compiled static-shape graphs see a constant
 * edge count. The engine applies it to its own scratch copy; the
 * caller's buffers are never mutated.
 *
 * Conventions:
 *
 * - Return values are three-valued: 0 success, positive recoverable
 *   (the evaluation failed at this input, the caller may back off or
 *   evaluate elsewhere), negative fatal (the instance is unusable).
 * - Capability discovery is by symbol presence: entry points marked
 *   optional are probed with dlsym. An absent symbol means the engine
 *   does not support the operation; a static (non-learning) engine
 *   simply omits the learning set.
 * - Engines with internal mutable state (learning surrogates) carry a
 *   monotone model-state counter. Any refit, hyperparameter change,
 *   or observation incorporation increments it; host-visible caches
 *   keyed on model output are valid only within one counter value.
 *   Acquisition (delegating a point to the oracle callback and
 *   returning its result verbatim) never increments the counter;
 *   incorporation happens only at declared boundaries (a batch call,
 *   a certify event, an explicit host hint), each of which may
 *   increment it and fire on_model_update.
 */
#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifdef RGPOT_ENGINE_BUILD
#define RGPOT_ENGINE_API __declspec(dllexport)
#else
#define RGPOT_ENGINE_API __declspec(dllimport)
#endif
#else
#define RGPOT_ENGINE_API __attribute__((visibility("default")))
#endif

typedef struct RgpotEnginePot RgpotEnginePot;

/**
 * Optional host-side transform applied by the engine to its scratch
 * copy of (positions, box) before evaluation. positions holds
 * 3*nAtoms doubles, box 9 doubles (row-major). Return 0 on success;
 * nonzero aborts the force call with that code.
 */
typedef int (*rgpot_engine_coord_transform)(void *user, long nAtoms,
                                            double *positions, double *box);

#define RGPOT_ENGINE_ABI_VERSION 1

RGPOT_ENGINE_API int rgpot_engine_abi_version(void);

RGPOT_ENGINE_API int rgpot_engine_available(void);

/**
 * Create a potential from a Cap'n Proto flat-array message whose root
 * struct is the engine's params arm in Potentials.capnp (UmaParams for
 * the UMA engine). config points at config_len bytes (word-aligned
 * capnp framing as produced by capnp::messageToFlatArray). Returns
 * NULL on failure with a message in errbuf.
 */
RGPOT_ENGINE_API RgpotEnginePot *rgpot_engine_create(const void *config,
                                                     size_t config_len,
                                                     char *errbuf,
                                                     size_t errlen);

/**
 * Gate settings for a learning engine, passed as the config blob.
 *
 * INTERIM. The config arm is specified as a Cap'n Proto message and
 * becomes GprParams once Potentials.capnp carries a nested
 * PotentialSpec (rgpot-b7zw). Until then a learning engine has no way
 * to be configured at all, which leaves certifyBelowForce at zero and
 * so leaves a driver's convergence claim unbacked by ground truth --
 * the one guarantee such an engine exists to provide. A fixed-layout
 * struct with a magic word closes that hole without pretending to be
 * capnp: an engine dispatches on the magic, so the capnp message can
 * replace it later without either side guessing which it received.
 */
#define RGPOT_ENGINE_GPR_GATE_MAGIC 0x47505247u /* "GRPG" */
#define RGPOT_ENGINE_GPR_GATE_VERSION 1u

typedef struct {
  unsigned magic;   /**< RGPOT_ENGINE_GPR_GATE_MAGIC */
  unsigned version; /**< RGPOT_ENGINE_GPR_GATE_VERSION */
  /** Answer from the model below this multiple of its likelihood noise. */
  double sigmaLoFactor;
  /** Above this multiple the geometry is not learned from. */
  double sigmaHiFactor;
  /** Predicted force magnitude under this (eV/Angstrom) forces a
   *  ground-truth call, so a driver's convergence claim is checked
   *  rather than trusted. Non-positive disables the certificate. */
  double certifyBelowForce;
  /** Refit cadence in training-set appends; zero refits at every
   *  incorporation boundary. */
  unsigned long long refitTrainingRevisions;
} RgpotEngineGprGate;

RGPOT_ENGINE_API void rgpot_engine_destroy(RgpotEnginePot *pot);

/**
 * Energy + forces. transform (with transform_user) may be NULL;
 * variance may be NULL. Returns 0 on success.
 */
RGPOT_ENGINE_API int rgpot_engine_force(
    RgpotEnginePot *pot, long nAtoms, const double *positions,
    const int *atomicNrs, double *forces, double *energy, double *variance,
    const double *box, rgpot_engine_coord_transform transform,
    void *transform_user);

/* ------------------------------------------------------------------ */
/* Optional symbols: probe with dlsym; absence means unsupported.     */
/* ------------------------------------------------------------------ */

/**
 * Strided-uniform batched evaluation: nSystems systems sharing one
 * composition (atomicNrs, nAtoms). positions and forces hold
 * 3*nAtoms*nSystems doubles with system stride 3*nAtoms, energies
 * nSystems doubles, variances likewise (may be NULL), boxes
 * 9*nSystems doubles with system stride 9. The transform is applied
 * per system to the engine's scratch copy. For learning engines this
 * call is also an incorporation boundary. A grouped variant (batches
 * of uniform groups) is a separate symbol, never a flag here.
 */
RGPOT_ENGINE_API int rgpot_engine_force_batch(
    RgpotEnginePot *pot, long nSystems, long nAtoms, const double *positions,
    const int *atomicNrs, double *forces, double *energies, double *variances,
    const double *boxes, rgpot_engine_coord_transform transform,
    void *transform_user);

/**
 * Monotone model-state counter for engines with internal mutable
 * state. Increments on every refit / hyperparameter change /
 * observation incorporation; host caches keyed on model output are
 * valid only while the value is unchanged. Static engines omit the
 * symbol (equivalent to a counter frozen at 0).
 */
RGPOT_ENGINE_API unsigned long long rgpot_engine_model_state(
    const RgpotEnginePot *pot);

/** Why a learning engine delegated a point to the oracle. */
typedef enum {
  RGPOT_ACQUIRE_VARIANCE = 1, /**< predictive variance above sigma_lo */
  RGPOT_ACQUIRE_CERTIFY = 2,  /**< verifying a surrogate convergence claim */
  RGPOT_ACQUIRE_BATCH = 3,    /**< batched acquisition sweep */
  RGPOT_ACQUIRE_REJECT = 4    /**< geometry above sigma_hi, evaluated but
                                   not incorporated */
} rgpot_engine_acquire_reason;

/**
 * Host callbacks for learning engines. oracle evaluates the ground
 * truth (same conventions and three-valued return as
 * rgpot_engine_force; variance is absent because the oracle is
 * exact). on_acquire fires when the engine delegates a point;
 * on_model_update fires after each incorporation with the new
 * model-state value and marks the boundary at which history-based
 * host state (quasi-Newton memory, cached factorizations) must be
 * reset. Either notification pointer may be NULL.
 */
typedef struct {
  void *ctx;
  int (*oracle)(void *ctx, long nAtoms, const double *positions,
                const int *atomicNrs, double *forces, double *energy,
                const double *box);
  void (*on_acquire)(void *ctx, rgpot_engine_acquire_reason reason);
  void (*on_model_update)(void *ctx, unsigned long long model_state);
} RgpotEngineCallbacks;

/**
 * Register host callbacks. The struct is copied; ctx must outlive the
 * potential or a subsequent call replacing it. Learning engines
 * require an oracle before the first evaluation and fail create-time
 * validation without one.
 */
RGPOT_ENGINE_API int rgpot_engine_set_callbacks(
    RgpotEnginePot *pot, const RgpotEngineCallbacks *callbacks);

/**
 * Runtime hint/tunable message: a Cap'n Proto flat-array whose root is
 * the engine's hint arm on the shared Potentials.capnp wire (never a
 * string-keyed channel). Also the host-initiated incorporation
 * boundary for learning engines. The engine validates the resulting
 * configuration as a whole before it takes effect; on failure the
 * previous configuration stays active and errbuf holds the reason.
 */
RGPOT_ENGINE_API int rgpot_engine_hint(RgpotEnginePot *pot, const void *msg,
                                       size_t msg_len, char *errbuf,
                                       size_t errlen);

/**
 * Neighbor-list declaration. An engine exporting this symbol expects
 * the host to provision every declared list before evaluation; an
 * engine omitting it self-provisions (the UMA molecular-box engine
 * builds its own complete graph). strict zero permits over-inclusion.
 */
typedef struct {
  double cutoff;  /**< Angstrom */
  int full_list;  /**< 1 full, 0 half */
  int strict;     /**< 0: over-inclusion tolerated */
} RgpotEngineNeighborOptions;

/**
 * Writes up to capacity declarations into opts and returns the total
 * count (call with capacity 0 to size). Declarations may change after
 * on_model_update; hosts re-query on that boundary.
 */
RGPOT_ENGINE_API long rgpot_engine_neighbor_options(
    const RgpotEnginePot *pot, RgpotEngineNeighborOptions *opts,
    long capacity);

#ifdef __cplusplus
}
#endif

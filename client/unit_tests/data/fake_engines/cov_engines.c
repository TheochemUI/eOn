/* Stand-in engines for the runtime loaders. Behavior is selected with
   environment variables so one library can take the success path, the
   version mismatch, a failed create, and a failed force. */
#include <stdlib.h>
#include <string.h>

static int env_int(const char *key, int fallback) {
  const char *v = getenv(key);
  if (v == NULL || v[0] == '\0') {
    return fallback;
  }
  return atoi(v);
}

static int fail_create(const char *key, char *err, size_t errlen) {
  if (!env_int(key, 0)) {
    return 0;
  }
  if (err != NULL && errlen > 0) {
    strncpy(err, "stand-in create failed", errlen - 1);
    err[errlen - 1] = '\0';
  }
  return 1;
}

static int fill_force(const char *key, long n, double *f, double *u,
                      double *var) {
  long i;
  if (env_int(key, 0) != 0) {
    return 1;
  }
  if (u != NULL) {
    *u = 0.5;
  }
  if (var != NULL) {
    *var = 0.25;
  }
  if (f != NULL) {
    for (i = 0; i < 3 * n; ++i) {
      f[i] = 0.1;
    }
  }
  return 0;
}

static int live = 1;

int rgpot_xtb_abi_version(void) { return env_int("EON_FAKE_XTB_ABI", 1); }

void *rgpot_xtb_create(const void *cfg, char *err, size_t errlen) {
  (void)cfg;
  if (fail_create("EON_FAKE_XTB_CREATE_FAIL", err, errlen)) {
    return NULL;
  }
  return &live;
}

void rgpot_xtb_destroy(void *pot) { (void)pot; }

int rgpot_xtb_force(void *pot, long n, const double *r, const int *z, double *f,
                    double *u, double *var, const double *box) {
  (void)pot;
  (void)r;
  (void)z;
  (void)box;
  return fill_force("EON_FAKE_XTB_FORCE_RC", n, f, u, var);
}

int rgpot_mta_abi_version(void) { return env_int("EON_FAKE_MTA_ABI", 1); }

void *rgpot_mta_create(const void *cfg, char *err, size_t errlen) {
  (void)cfg;
  if (fail_create("EON_FAKE_MTA_CREATE_FAIL", err, errlen)) {
    return NULL;
  }
  return &live;
}

void rgpot_mta_destroy(void *pot) { (void)pot; }

int rgpot_mta_force(void *pot, long n, const double *r, const int *z, double *f,
                    double *u, double *var, const double *box) {
  (void)pot;
  (void)r;
  (void)z;
  (void)box;
  return fill_force("EON_FAKE_MTA_FORCE_RC", n, f, u, var);
}

int rgpot_engine_abi_version(void) { return env_int("EON_FAKE_ENGINE_ABI", 1); }

void *rgpot_engine_create(const void *cfg, size_t len, char *err,
                          size_t errlen) {
  (void)cfg;
  (void)len;
  if (fail_create("EON_FAKE_ENGINE_CREATE_FAIL", err, errlen)) {
    return NULL;
  }
  return &live;
}

void rgpot_engine_destroy(void *pot) { (void)pot; }

int rgpot_engine_force(void *pot, long n, const double *r, const int *z,
                       double *f, double *u, double *var, const double *box,
                       void *transform, void *user) {
  (void)pot;
  (void)r;
  (void)z;
  (void)box;
  (void)transform;
  (void)user;
  return fill_force("EON_FAKE_ENGINE_FORCE_RC", n, f, u, var);
}

int eon_mta_abi_version(void) { return env_int("EON_FAKE_EON_MTA_ABI", 1); }

void *eon_mta_pot_create(const void *cfg, char *err, size_t errlen) {
  (void)cfg;
  if (fail_create("EON_FAKE_EON_MTA_CREATE_FAIL", err, errlen)) {
    return NULL;
  }
  return &live;
}

void eon_mta_pot_destroy(void *pot) { (void)pot; }

int eon_mta_pot_force(void *pot, long n, const double *r, const int *z,
                      double *f, double *u, double *var, const double *box) {
  (void)pot;
  (void)r;
  (void)z;
  (void)box;
  return fill_force("EON_FAKE_EON_MTA_FORCE_RC", n, f, u, var);
}

int eon_covplug_marker(void) { return 7; }

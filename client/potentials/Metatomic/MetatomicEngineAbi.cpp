/*
** libmetatomic_engine.so — C ABI for RGPOT-metatomic (wraps
*MetatomicPotential).
** Fat eOn still uses MetatomicPotential directly via potential=Metatomic.
*/
#define RGPOT_MTA_ENGINE_BUILD
#include "eon/Parameters.h"
#include "eon/potentials/Metatomic/MetatomicPotential.h"
#include "eon/potentials/Rgpot/metatomic_c_abi.h"
#include <cstdio>
#include <cstring>
#include <exception>
#include <memory>
#include <string>
#include <vector>

struct RgpotMtaPot {
  std::unique_ptr<MetatomicPotential> impl;
};

static const char *nz(const char *s) { return s ? s : ""; }
static void set_err(char *b, size_t n, const char *m) {
  if (b && n)
    std::snprintf(b, n, "%s", m ? m : "unknown");
}

extern "C" {

int rgpot_mta_abi_version(void) { return RGPOT_MTA_ABI_VERSION; }
int rgpot_mta_available(void) { return 1; }

RgpotMtaPot *rgpot_mta_create(const RgpotMtaConfig *cfg, char *errbuf,
                              size_t errlen) {
  if (!cfg || !cfg->model_path || !cfg->model_path[0]) {
    set_err(errbuf, errlen, "model_path required");
    return nullptr;
  }
  try {
    Parameters p;
    p.potential_options.potential = PotType::METATOMIC;
    auto &o = p.metatomic_options;
    o.model_path = cfg->model_path;
    o.device = nz(cfg->device)[0] ? cfg->device : "cpu";
    o.length_unit = nz(cfg->length_unit)[0] ? cfg->length_unit : "angstrom";
    o.extensions_directory = nz(cfg->extensions_directory);
    o.check_consistency = cfg->check_consistency != 0;
    o.uncertainty_threshold = cfg->uncertainty_threshold;
    o.deterministic = !cfg->torch_determinism_strict; // map: strict => det
    o.deterministic_strict = cfg->torch_determinism_strict != 0;
    auto pot = std::make_unique<RgpotMtaPot>();
    pot->impl = std::make_unique<MetatomicPotential>(p);
    return pot.release();
  } catch (const std::exception &e) {
    set_err(errbuf, errlen, e.what());
    return nullptr;
  } catch (...) {
    set_err(errbuf, errlen, "create failed");
    return nullptr;
  }
}

void rgpot_mta_destroy(RgpotMtaPot *pot) { delete pot; }

int rgpot_mta_force(RgpotMtaPot *pot, long nAtoms, const double *positions,
                    const int *atomicNrs, double *forces, double *energy,
                    double *variance, const double *box) {
  if (!pot || !pot->impl || !positions || !atomicNrs || !forces || !energy ||
      !box)
    return 1;
  try {
    pot->impl->force(nAtoms, positions, atomicNrs, forces, energy, variance,
                     box);
    return 0;
  } catch (...) {
    return 2;
  }
}

} // extern "C"

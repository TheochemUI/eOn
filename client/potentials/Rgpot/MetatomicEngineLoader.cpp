#include "MetatomicEngineLoader.h"
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>
#ifndef _WIN32
#include <dlfcn.h>
#endif

namespace {
void *open_lib(const char *path) {
#ifndef _WIN32
  return dlopen(path, RTLD_NOW | RTLD_GLOBAL);
#else
  return nullptr;
#endif
}
void close_lib(void *h) {
#ifndef _WIN32
  if (h)
    dlclose(h);
#endif
}
void *load_sym(void *h, const char *n) {
#ifndef _WIN32
  return dlsym(h, n);
#else
  return nullptr;
#endif
}
} // namespace

MetatomicEngineLoader::MetatomicEngineLoader(const MetatomicEngineOptions &opt) {
  std::vector<std::string> paths;
  if (!opt.engine_path.empty())
    paths.push_back(opt.engine_path);
  if (const char *e = std::getenv("RGPOT_METATOMIC_ENGINE"))
    if (e && *e)
      paths.emplace_back(e);
  if (const char *e = std::getenv("METATOMIC_ENGINE"))
    if (e && *e)
      paths.emplace_back(e);
  paths.emplace_back("libmetatomic_engine.so");
  auto add_dirs = [&](const char *env) {
    if (!env)
      return;
    std::string s(env);
    size_t i = 0;
    while (i < s.size()) {
      auto j = s.find(':', i);
      if (j == std::string::npos)
        j = s.size();
      if (j > i) {
        std::string d = s.substr(i, j - i);
        if (!d.empty() && d.back() != '/')
          d += '/';
        paths.push_back(d + "libmetatomic_engine.so");
      }
      i = j + 1;
    }
  };
  add_dirs(std::getenv("EON_POTENTIALS_PATH"));
  add_dirs(std::getenv("RGPOT_ENGINE_PATH"));

  std::string last_dlerr;
  for (const auto &p : paths) {
    m_lib = open_lib(p.c_str());
    if (m_lib)
      break;
#ifndef _WIN32
    if (const char *e = dlerror())
      last_dlerr = e;
#endif
  }
  if (!m_lib) {
    std::string msg =
        "RGPOT(metatomic): libmetatomic_engine.so not found "
        "(set RGPOT_METATOMIC_ENGINE or EON_POTENTIALS_PATH)";
    if (!last_dlerr.empty())
      msg += std::string("; last dlerror: ") + last_dlerr;
    throw std::runtime_error(msg);
  }

  auto abi = reinterpret_cast<int (*)()>(load_sym(m_lib, "rgpot_mta_abi_version"));
  m_create = reinterpret_cast<create_fn>(load_sym(m_lib, "rgpot_mta_create"));
  m_destroy = reinterpret_cast<destroy_fn>(load_sym(m_lib, "rgpot_mta_destroy"));
  m_force = reinterpret_cast<force_fn>(load_sym(m_lib, "rgpot_mta_force"));
  if (!abi || !m_create || !m_destroy || !m_force ||
      abi() != RGPOT_MTA_ABI_VERSION) {
    close_lib(m_lib);
    m_lib = nullptr;
    throw std::runtime_error("RGPOT(metatomic): engine C ABI missing/mismatch");
  }
  char err[1024]{};
  RgpotMtaConfig cfg{};
  cfg.model_path = opt.model_path.c_str();
  cfg.device = opt.device.c_str();
  cfg.length_unit = opt.length_unit.c_str();
  cfg.extensions_directory = opt.extensions_directory.c_str();
  cfg.check_consistency = opt.check_consistency ? 1 : 0;
  cfg.uncertainty_threshold = opt.uncertainty_threshold;
  cfg.torch_determinism_strict = opt.torch_determinism_strict ? 1 : 0;
  m_pot = m_create(&cfg, err, sizeof err);
  if (!m_pot) {
    close_lib(m_lib);
    m_lib = nullptr;
    throw std::runtime_error(std::string("RGPOT(metatomic): create failed: ") +
                             err);
  }
}

MetatomicEngineLoader::~MetatomicEngineLoader() {
  // Destroy the pot while the engine is still mapped. Do not dlclose
  // libmetatomic_engine.so: torch/metatomic static teardown after dlclose
  // routinely SEGV on process exit.
  if (m_pot && m_destroy)
    m_destroy(m_pot);
  m_pot = nullptr;
  m_lib = nullptr;
}

void MetatomicEngineLoader::force(long N, const double *R, const int *atomicNrs,
                                  double *F, double *U, double *variance,
                                  const double *box) const {
  if (!m_pot || !m_force)
    throw std::runtime_error("RGPOT(metatomic): engine not available");
  double var = 0.0;
  const int rc = m_force(m_pot, N, R, atomicNrs, F, U, variance ? variance : &var,
                         box);
  if (rc != 0)
    throw std::runtime_error("RGPOT(metatomic): force failed");
}

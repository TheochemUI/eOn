#include "XTBEngineLoader.h"
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
  return dlopen(path, RTLD_NOW | RTLD_LOCAL);
#else
  (void)path;
  return nullptr;
#endif
}
void close_lib(void *h) {
#ifndef _WIN32
  if (h)
    dlclose(h);
#else
  (void)h;
#endif
}
void *load_sym(void *h, const char *n) {
#ifndef _WIN32
  return dlsym(h, n);
#else
  (void)h;
  (void)n;
  return nullptr;
#endif
}
} // namespace

XTBEngineLoader::XTBEngineLoader(const XTBEngineOptions &opt) {
  std::vector<std::string> paths;
  if (!opt.engine_path.empty())
    paths.push_back(opt.engine_path);
  if (const char *e = std::getenv("RGPOT_XTB_ENGINE"))
    if (e && *e)
      paths.emplace_back(e);
  if (const char *e = std::getenv("XTB_ENGINE"))
    if (e && *e)
      paths.emplace_back(e);
  paths.emplace_back("libxtb_engine.so");
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
        paths.push_back(d + "libxtb_engine.so");
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
        "RGPOT(xtb): libxtb_engine.so not found "
        "(set RGPOT_XTB_ENGINE or [RgpotPot] engine_path)";
    if (!last_dlerr.empty())
      msg += std::string("; last dlerror: ") + last_dlerr;
    throw std::runtime_error(msg);
  }

  auto abi =
      reinterpret_cast<int (*)()>(load_sym(m_lib, "rgpot_xtb_abi_version"));
  m_create = reinterpret_cast<create_fn>(load_sym(m_lib, "rgpot_xtb_create"));
  m_destroy = reinterpret_cast<destroy_fn>(load_sym(m_lib, "rgpot_xtb_destroy"));
  m_force = reinterpret_cast<force_fn>(load_sym(m_lib, "rgpot_xtb_force"));
  if (!abi || !m_create || !m_destroy || !m_force ||
      abi() != RGPOT_XTB_ABI_VERSION) {
    close_lib(m_lib);
    m_lib = nullptr;
    throw std::runtime_error("RGPOT(xtb): engine C ABI missing/mismatch");
  }
  char err[1024]{};
  RgpotXtbConfig cfg{};
  cfg.method = opt.method;
  cfg.accuracy = opt.accuracy;
  cfg.electronic_temperature = opt.electronic_temperature;
  cfg.max_iterations = opt.max_iterations;
  cfg.charge = opt.charge;
  cfg.uhf = opt.uhf;
  m_pot = m_create(&cfg, err, sizeof err);
  if (!m_pot) {
    close_lib(m_lib);
    m_lib = nullptr;
    throw std::runtime_error(std::string("RGPOT(xtb): create failed: ") + err);
  }
}

XTBEngineLoader::~XTBEngineLoader() {
  if (m_pot && m_destroy)
    m_destroy(m_pot);
  m_pot = nullptr;
  // Safe to close: xtb has no torch-style static teardown hazard.
  if (m_lib)
    close_lib(m_lib);
  m_lib = nullptr;
}

void XTBEngineLoader::force(long N, const double *R, const int *atomicNrs,
                            double *F, double *U, double *variance,
                            const double *box) const {
  if (!m_pot || !m_force)
    throw std::runtime_error("RGPOT(xtb): engine not available");
  double var = 0.0;
  const int rc = m_force(m_pot, N, R, atomicNrs, F, U,
                         variance ? variance : &var, box);
  if (rc != 0)
    throw std::runtime_error("RGPOT(xtb): force failed");
}

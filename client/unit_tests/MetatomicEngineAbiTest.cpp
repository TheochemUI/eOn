#include "catch2/catch_amalgamated.hpp"
#include "eon/DynLib.h"

#include <cstdlib>
#include <string>
#include <vector>

// Structural + optional live check: libmetatomic_engine exports rgpot_mta_*
// ABI. Uses DynLib (LoadLibrary on Windows, dlopen elsewhere) - not raw
// dlfcn.h - so the test builds and runs on MSVC.
TEST_CASE("metatomic engine C ABI symbols present when engine is loadable",
          "[metatomic][rgpot][engine]") {
  const char *env = std::getenv("RGPOT_METATOMIC_ENGINE");
  std::vector<std::string> candidates;
  if (env && *env) {
    candidates.emplace_back(env);
  }
  // Platform default sonames / import names.
  candidates.emplace_back("libmetatomic_engine.so");
  candidates.emplace_back("libmetatomic_engine.dylib");
  candidates.emplace_back("metatomic_engine.dll");
  candidates.emplace_back("libmetatomic_engine.dll");

  // Build-tree relative path from test env (colon/semicolon separated).
  const char *pp = std::getenv("EON_POTENTIALS_PATH");
  if (pp && *pp) {
    std::string d(pp);
#ifdef _WIN32
    const char sep = ';';
    const char slash = '\\';
#else
    const char sep = ':';
    const char slash = '/';
#endif
    auto pos = d.find(sep);
    if (pos != std::string::npos)
      d = d.substr(0, pos);
    if (!d.empty() && d.back() != '/' && d.back() != '\\')
      d += slash;
    candidates.push_back(d + "libmetatomic_engine.so");
    candidates.push_back(d + "libmetatomic_engine.dylib");
    candidates.push_back(d + "metatomic_engine.dll");
    candidates.push_back(d + "libmetatomic_engine.dll");
  }

  eonc::dynlib::Handle h{};
  for (const auto &path : candidates) {
    h = eonc::dynlib::open(path.c_str());
    if (h)
      break;
  }
  if (!h) {
    WARN("libmetatomic_engine not on path; skip live ABI check");
    return;
  }

  REQUIRE(eonc::dynlib::sym(h, "rgpot_mta_abi_version") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "rgpot_mta_create") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "rgpot_mta_destroy") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "rgpot_mta_force") != nullptr);
  auto abi = eonc::dynlib::loadSym<int (*)()>(h, "rgpot_mta_abi_version");
  REQUIRE(abi != nullptr);
  REQUIRE(abi() == 1);
  eonc::dynlib::close(h);
}

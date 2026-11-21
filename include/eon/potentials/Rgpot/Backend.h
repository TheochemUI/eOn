#pragma once

#include <cctype>
#include <cstdlib>
#include <string>

namespace eonc {

inline std::string canonicalRgpotBackend(std::string backend) {
  for (auto &c : backend)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (backend == "nwchem" || backend == "nwchempot")
    return "nwchemc";
  if (backend == "cpmd" || backend == "cpmdpot")
    return "cpmdc";
  if (backend == "mta" || backend == "metatomicpot")
    return "metatomic";
  if (backend == "omol" || backend == "umapot")
    return "uma";
  if (backend == "gfn" || backend == "gfnxtb" || backend == "xtbpot")
    return "xtb";
  return backend;
}

inline std::string configuredRgpotBackend(const std::string &backend) {
  const auto *overrideName = std::getenv("RGPOT_BACKEND");
  return canonicalRgpotBackend(overrideName ? overrideName : backend);
}

} // namespace eonc

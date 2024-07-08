#include "Parser.hpp"

#include <filesystem>
#include <spdlog/spdlog.h>

namespace fs = std::filesystem;

namespace eonc {

toml::table loadTOML(const std::string &fname) {
  fs::path fpath(fname);

  if (fpath.extension() != ".toml") {
    throw std::runtime_error(
        fmt::format("Filename must have a .toml extension\n Got: {}", fname));
  }

  if (!fs::exists(fpath)) {
    throw std::runtime_error(
        fmt::format("Input filename does not exist\n Got: {}", fname));
  }

  try {
    toml::parse_result result = toml::parse_file(fname);
    return result;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format("Error parsing the configuration file: {}", e.what()));
  }
}

} // namespace eonc

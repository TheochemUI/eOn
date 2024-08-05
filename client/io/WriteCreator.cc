#include "client/Parser.hpp"

#include "client/io/WriteCreator.hpp"
#include "client/io/writers/con/ConWriter.hpp"
#include "client/io/writers/con/ConvelWriter.hpp"
#include "client/io/writers/tibble/TibbleWriter.hpp"
#include "client/io/writers/xyz/XYZWriter.hpp"

namespace eonc::io {
std::unique_ptr<WriteBase> mkWriter(const toml::table &config) {
  config_section(config, "Main");
  auto wtype = get_enum_toml<WriteType>(config["Main"]["write"]);
  switch (wtype) {
  case WriteType::CON: {
    return (std::make_unique<ConWriter>());
    break;
  }
  case WriteType::CONVEL: {
    return (std::make_unique<ConvelWriter>());
    break;
  }
  case WriteType::XYZ: {
    return (std::make_unique<XYZWriter>());
    break;
  }
  case WriteType::TIBBLE: {
    return (std::make_unique<TibbleWriter>());
    break;
  }
  default:
    SPDLOG_ERROR("No known writer could be constructed from {}",
                 magic_enum::enum_name(wtype));
    throw std::runtime_error("Terminating");
    break;
  }
}
} // namespace eonc::io

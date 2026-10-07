#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"
#include "eon/potentials/SocketNWChem/SocketNWChemPot.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

TEST_CASE("SocketNWChem listens and writes a template", "[socketnwchem]") {
  eonc::Parameters params;
  eonc::ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::SocketNWChem;
  auto &opt = eonc::ParametersLoadAccess::socket_nwchem_options(params);
  opt.unix_socket_mode = true;
  opt.unix_socket_path = "eon_sock_t";
  opt.mem_in_gb = 1;
  opt.nwchem_settings = "theory.nw";

  SocketNWChemPot pot(params);
  const auto path =
      std::filesystem::temp_directory_path() / "eon-nwchem-template.nw";
  pot.write_nwchem_template(path.string(), 2, {"H", "O"});
  std::ifstream in(path);
  std::stringstream text;
  text << in.rdbuf();
  REQUIRE(text.str().find("memory 1 gb") != std::string::npos);
  REQUIRE(text.str().find("geometry units bohr") != std::string::npos);
  std::filesystem::remove(path);
}

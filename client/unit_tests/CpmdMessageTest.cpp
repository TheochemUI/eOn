#include "eon/potentials/Rgpot/CpmdMessage.h"
#include "eon/potentials/Rgpot/RGPotEngine.h"
#include "catch2/catch_amalgamated.hpp"

#include <capnp/message.h>
#include <capnp/serialize.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

using namespace Catch::Matchers;

namespace tests {
namespace {

void clear_env(const char *name) {
#ifdef _WIN32
  _putenv_s(name, "");
#else
  unsetenv(name);
#endif
}

std::filesystem::path write_message(const std::filesystem::path &path,
                                    ::capnp::MallocMessageBuilder &src) {
  auto words = capnp::messageToFlatArray(src);
  std::ofstream out(path, std::ios::binary);
  out.write(reinterpret_cast<const char *>(words.begin()),
            static_cast<std::streamsize>(words.size() * sizeof(capnp::word)));
  REQUIRE(out);
  return path;
}

} // namespace

TEST_CASE("NWChem DFT theory is matched without regard to case",
          "[params][nwchem][RGPOT]") {
  const std::string mixed = nwchemDftInputBlock("Dft", "B3LYP", 1, "");
  REQUIRE(mixed.find("xc B3LYP") != std::string::npos);
  REQUIRE(mixed.find("mult 1") != std::string::npos);
  REQUIRE(nwchemDftInputBlock("DFT", "pbe", 3, "").find("xc pbe") !=
          std::string::npos);
  REQUIRE(nwchemDftInputBlock("B3LYP", "rhf", 1, "").find("xc B3LYP") !=
          std::string::npos);
  REQUIRE(nwchemDftInputBlock("dft", "b3lyp", 1, "keep") == "keep");
  REQUIRE(nwchemDftInputBlock("scf", "rhf", 1, "").empty());
}

TEST_CASE("input_block appends and leaves params_path sections",
          "[params][cpmd][RGPOT]") {
  clear_env("RGPOT_CPMD_INPUT_BLOCK");
  clear_env("RGPOT_PARAMS_PATH");
  namespace fs = std::filesystem;
  const auto dir = fs::temp_directory_path() /
                   ("eon_cpmd_layers_" + std::to_string(std::rand()));
  fs::create_directories(dir);

  ::capnp::MallocMessageBuilder src;
  auto root = src.initRoot<::CPMDParams>();
  root.setFunctional("PBE");
  root.setCutOffRy(55.5);
  root.setCharge(4);
  auto sections = root.initInputSections(2);
  sections[0].setRaw("SECTION_A");
  sections[1].setRaw("SECTION_B");
  auto prior = root.initInputBlocks(1);
  prior.set(0, "BLOCK_FROM_FILE");
  const auto path = write_message(dir / "two-sections.bin", src);

  RGPotEngineOptions opt;
  opt.backend = "cpmdc";
  opt.params_path = path.string();
  opt.functional = "LDA";
  opt.cutoff_ry = 10.0;
  opt.charge = 1;
  opt.scratch_dir = "/scratch/run";
  opt.input_block = "FROM_INI";

  ::capnp::MallocMessageBuilder msg;
  auto built = eon::fillCpmdParams(msg, opt);
  REQUIRE(built.getInputSections().size() == 2);
  REQUIRE(std::string(built.getInputSections()[0].getRaw().cStr()) ==
          "SECTION_A");
  REQUIRE(std::string(built.getInputSections()[1].getRaw().cStr()) ==
          "SECTION_B");
  REQUIRE(built.getInputBlocks().size() == 2);
  REQUIRE(std::string(built.getInputBlocks()[0].cStr()) == "BLOCK_FROM_FILE");
  REQUIRE(std::string(built.getInputBlocks()[1].cStr()) == "FROM_INI");
  REQUIRE(std::string(built.getFunctional().cStr()) == "PBE");
  REQUIRE(built.getCutOffRy() == Catch::Approx(55.5));
  REQUIRE(built.getCharge() == 4);
  REQUIRE(std::string(built.getScratchDir().cStr()) == "/scratch/run");

  fs::remove_all(dir);
}

TEST_CASE("an empty input_block keeps the sections from params_path",
          "[params][cpmd][RGPOT]") {
  clear_env("RGPOT_CPMD_INPUT_BLOCK");
  namespace fs = std::filesystem;
  const auto dir = fs::temp_directory_path() /
                   ("eon_cpmd_keep_" + std::to_string(std::rand()));
  fs::create_directories(dir);

  ::capnp::MallocMessageBuilder src;
  auto root = src.initRoot<::CPMDParams>();
  auto sections = root.initInputSections(2);
  sections[0].setRaw("SECTION_A");
  sections[1].setRaw("SECTION_B");
  const auto path = write_message(dir / "sections-only.bin", src);

  RGPotEngineOptions opt;
  opt.params_path = path.string();
  opt.input_block = "";
  opt.functional = "LDA";

  ::capnp::MallocMessageBuilder msg;
  auto built = eon::fillCpmdParams(msg, opt);
  REQUIRE(built.getInputSections().size() == 2);
  REQUIRE(built.getInputBlocks().size() == 0);
  REQUIRE(std::string(built.getFunctional().cStr()) == "BLYP");

  fs::remove_all(dir);
}

TEST_CASE("scalars fill the message when params_path is empty",
          "[params][cpmd][RGPOT]") {
  clear_env("RGPOT_CPMD_INPUT_BLOCK");
  RGPotEngineOptions opt;
  opt.functional = "PBE";
  opt.cutoff_ry = 40.0;
  opt.charge = 2;
  opt.multiplicity = 3;
  opt.input_block = "ONLY_BLOCK";

  ::capnp::MallocMessageBuilder msg;
  auto built = eon::fillCpmdParams(msg, opt);
  REQUIRE(std::string(built.getFunctional().cStr()) == "PBE");
  REQUIRE(built.getCutOffRy() == Catch::Approx(40.0));
  REQUIRE(built.getCharge() == 2);
  REQUIRE(built.getMultiplicity() == 3);
  REQUIRE(built.getInputSections().size() == 0);
  REQUIRE(built.getInputBlocks().size() == 1);
  REQUIRE(std::string(built.getInputBlocks()[0].cStr()) == "ONLY_BLOCK");
}

} // namespace tests

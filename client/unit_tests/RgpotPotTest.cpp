#include "catch2/catch_amalgamated.hpp"
#include "eon/MatrixHelpers.hpp"
#include "eon/Matter.h"
#include "eonc_test_aliases.hpp"
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

using namespace Catch::Matchers;

namespace tests {

namespace {

bool env_nonempty(const char *name) {
  const char *v = std::getenv(name);
  return v != nullptr && v[0] != '\0';
}

void clear_env(const char *name) {
#ifdef _WIN32
  _putenv_s(name, "");
#else
  unsetenv(name);
#endif
}

} // namespace

TEST_CASE("RgpotPot in-process nwchemc force (no potserv)",
          "[PotTest][RGPOT][nwchemc]") {
#ifndef WITH_RGPOT
  SKIP("built without WITH_RGPOT");
#else
  // Engine is a runtime dlopen dep; packaging CI may build WITH_RGPOT without
  // a local libnwchemc — skip rather than fail the default suite.
  if (!(env_nonempty("NWCHEMC_LIBRARY") ||
        env_nonempty("RGPOT_NWCHEMC_ENGINE") ||
        env_nonempty("RGPOT_NWCHEM_ENGINE"))) {
    SKIP("no NWChem embed library in environment "
         "(set NWCHEMC_LIBRARY / RGPOT_NWCHEMC_ENGINE)");
  }

  Parameters params{};
  ParametersLoadAccess::potential_options(params).potential = PotType::RGPOT;
  ParametersLoadAccess::rgpot_options(params).backend = "nwchemc";
  ParametersLoadAccess::rgpot_options(params).basis = "sto-3g";
  ParametersLoadAccess::rgpot_options(params).theory = "scf";
  ParametersLoadAccess::rgpot_options(params).scf_type = "rhf";
  ParametersLoadAccess::rgpot_options(params).charge = 0;
  ParametersLoadAccess::rgpot_options(params).multiplicity = 1;

  auto pot = eonc::helpers::sharePotential(eonc::helpers::makePotential(
      params.potential_options().potential, params));
  REQUIRE(pot != nullptr);
  REQUIRE(pot->getType() == PotType::RGPOT);

  auto matter = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(matter->con2matter("pos.con")));
  REQUIRE(matter->numberOfAtoms() == 9);

  double energy = 0.0;
  AtomMatrix forces = MatrixXd::Zero(matter->numberOfAtoms(), 3);
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy, nullptr,
             matter->getCell().data());

  REQUIRE(std::isfinite(energy));
  REQUIRE(std::abs(energy) > 1e-6);

  // Two successive forces must both succeed (warm multi-call path)
  double energy2 = 0.0;
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy2, nullptr,
             matter->getCell().data());
  REQUIRE(std::isfinite(energy2));
  REQUIRE(std::abs(energy - energy2) < 1e-4);
#endif
}

TEST_CASE("RgpotPot in-process cpmdc force (no potserv)",
          "[PotTest][RGPOT][cpmdc]") {
#ifndef WITH_RGPOT
  SKIP("built without WITH_RGPOT");
#else
  // Engine is a runtime dlopen dep; packaging CI may build WITH_RGPOT without
  // a local libcpmdc — skip rather than fail the default suite.
  if (!(env_nonempty("CPMDC_LIBRARY") || env_nonempty("RGPOT_CPMDC_ENGINE") ||
        env_nonempty("RGPOT_CPMD_ENGINE"))) {
    SKIP("no CPMD embed library in environment "
         "(set CPMDC_LIBRARY / RGPOT_CPMDC_ENGINE)");
  }

  Parameters params{};
  ParametersLoadAccess::potential_options(params).potential = PotType::RGPOT;
  ParametersLoadAccess::rgpot_options(params).backend = "cpmdc";
  ParametersLoadAccess::rgpot_options(params).functional = "BLYP";
  ParametersLoadAccess::rgpot_options(params).cutoff_ry = 70.0;
  ParametersLoadAccess::rgpot_options(params).charge = 0;
  ParametersLoadAccess::rgpot_options(params).multiplicity = 1;

  auto pot = eonc::helpers::sharePotential(eonc::helpers::makePotential(
      params.potential_options().potential, params));
  REQUIRE(pot != nullptr);
  REQUIRE(pot->getType() == PotType::RGPOT);

  auto matter = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(matter->con2matter("pos.con")));
  REQUIRE(matter->numberOfAtoms() == 9);

  double energy = 0.0;
  AtomMatrix forces = MatrixXd::Zero(matter->numberOfAtoms(), 3);
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy, nullptr,
             matter->getCell().data());

  REQUIRE(std::isfinite(energy));
  REQUIRE(std::abs(energy) > 1e-6);

  // Two successive forces must both succeed (warm multi-call path)
  double energy2 = 0.0;
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy2, nullptr,
             matter->getCell().data());
  REQUIRE(std::isfinite(energy2));
  REQUIRE(std::abs(energy - energy2) < 1e-4);
#endif
}

TEST_CASE("RgpotPot reads params_path and ranks_per_image from INI",
          "[params][ini][RGPOT]") {
  Parameters p;
  REQUIRE(p.rgpot_options().params_path.empty());
  REQUIRE(p.rgpot_options().ranks_per_image == 0);
  REQUIRE(p.load_ini_text("[Potential]\npotential = rgpot\n\n"
                          "[RgpotPot]\nbackend = cpmdc\n"
                          "params_path = /data/si3n4.bin\n"
                          "ranks_per_image = 6\n") == 0);
  REQUIRE(p.potential_options().potential == PotType::RGPOT);
  REQUIRE(p.rgpot_options().params_path == "/data/si3n4.bin");
  REQUIRE(p.rgpot_options().ranks_per_image == 6);
}

TEST_CASE("cpmd section supplies cutOffRy and overrides RgpotPot",
          "[params][ini][RGPOT]") {
  Parameters p;
  REQUIRE(p.load_ini_text("[Potential]\npotential = rgpot\n\n"
                          "[RgpotPot]\nbackend = cpmdc\n"
                          "cutOffRy = 10\ncharge = 1\n"
                          "input_block = FROM_SHARED\n"
                          "params_path = /data/si3n4.bin\n\n"
                          "[cpmd]\nfunctional = PBE\n"
                          "cutOffRy = 55.5\ncharge = 4\n"
                          "input_block = DEMO_BLOCK\n") == 0);
  REQUIRE(p.rgpot_options().functional == "PBE");
  REQUIRE(p.rgpot_options().cutoff_ry == Catch::Approx(55.5));
  REQUIRE(p.rgpot_options().charge == 4);
  REQUIRE(p.rgpot_options().input_block == "DEMO_BLOCK");
  REQUIRE(p.rgpot_options().params_path == "/data/si3n4.bin");
}

TEST_CASE("cpmd section still loads the older cutoff keys",
          "[params][ini][RGPOT]") {
  Parameters legacy;
  REQUIRE(legacy.load_ini_text("[Potential]\npotential = rgpot\n\n"
                               "[RgpotPot]\nbackend = cpmdc\n\n"
                               "[cpmd]\ncutoff_ry = 40\n") == 0);
  REQUIRE(legacy.rgpot_options().cutoff_ry == Catch::Approx(40.0));

  Parameters older;
  REQUIRE(older.load_ini_text("[Potential]\npotential = rgpot\n\n"
                              "[RgpotPot]\nbackend = cpmdc\n\n"
                              "[cpmd]\ncpmd_cut_off_ry = 33\n") == 0);
  REQUIRE(older.rgpot_options().cutoff_ry == Catch::Approx(33.0));
}

TEST_CASE("nwchemc ignores the cpmd section", "[params][ini][RGPOT]") {
  Parameters p;
  REQUIRE(p.load_ini_text("[Potential]\npotential = rgpot\n\n"
                          "[RgpotPot]\nbackend = nwchemc\n\n"
                          "[cpmd]\ncutOffRy = 12.5\n"
                          "input_block = SHOULD_NOT_APPLY\n") == 0);
  REQUIRE(p.rgpot_options().cutoff_ry == Catch::Approx(70.0));
  REQUIRE(p.rgpot_options().input_block.empty());
}

TEST_CASE("RgpotPot cutOffRy wins over cutoff_ry and cpmd_cut_off_ry",
          "[params][ini][RGPOT]") {
  Parameters p;
  REQUIRE(p.load_ini_text("[Potential]\npotential = rgpot\n\n"
                          "[RgpotPot]\nbackend = nwchemc\n"
                          "cutOffRy = 22\n"
                          "cutoff_ry = 11\n"
                          "cpmd_cut_off_ry = 9\n") == 0);
  REQUIRE(p.rgpot_options().cutoff_ry == Catch::Approx(22.0));
}

TEST_CASE("nwchem_basis fills basis when basis is absent",
          "[params][ini][RGPOT]") {
  Parameters filled;
  REQUIRE(filled.load_ini_text("[Potential]\npotential = rgpot\n\n"
                               "[RgpotPot]\nbackend = nwchemc\n"
                               "nwchem_basis = 6-31g\n") == 0);
  REQUIRE(filled.rgpot_options().basis == "6-31g");

  Parameters both;
  REQUIRE(both.load_ini_text("[Potential]\npotential = rgpot\n\n"
                             "[RgpotPot]\nbackend = nwchemc\n"
                             "basis = sto-3g\n"
                             "nwchem_basis = 6-31g\n") == 0);
  REQUIRE(both.rgpot_options().basis == "sto-3g");
}

TEST_CASE("xtb_charge follows charge for backend nwchemc",
          "[params][ini][RGPOT]") {
  Parameters p;
  REQUIRE(p.load_ini_text("[Potential]\npotential = rgpot\n\n"
                          "[RgpotPot]\nbackend = nwchemc\n"
                          "charge = 4\n") == 0);
  REQUIRE(p.rgpot_options().xtb_charge == Catch::Approx(4.0));

  Parameters ignored;
  REQUIRE(ignored.load_ini_text("[Potential]\npotential = rgpot\n\n"
                                "[RgpotPot]\nbackend = nwchemc\n"
                                "charge = 4\n\n"
                                "[XTBPot]\ncharge = -2\n") == 0);
  REQUIRE(ignored.rgpot_options().xtb_charge == Catch::Approx(4.0));
}

TEST_CASE("XTBPot charge overlays xtb_charge for an xtb backend",
          "[params][ini][RGPOT]") {
  Parameters fromCharge;
  REQUIRE(fromCharge.load_ini_text("[Potential]\npotential = rgpot\n\n"
                                   "[RgpotPot]\nbackend = xtb\n"
                                   "charge = 3\n") == 0);
  REQUIRE(fromCharge.rgpot_options().xtb_charge == Catch::Approx(3.0));

  Parameters explicitKey;
  REQUIRE(explicitKey.load_ini_text("[Potential]\npotential = rgpot\n\n"
                                    "[RgpotPot]\nbackend = xtb\n"
                                    "charge = 3\n"
                                    "xtb_charge = 1.5\n") == 0);
  REQUIRE(explicitKey.rgpot_options().xtb_charge == Catch::Approx(1.5));

  Parameters overlay;
  REQUIRE(overlay.load_ini_text("[Potential]\npotential = rgpot\n\n"
                                "[RgpotPot]\nbackend = xtb\n"
                                "charge = 3\n"
                                "xtb_charge = 1.5\n\n"
                                "[XTBPot]\ncharge = -2\n") == 0);
  REQUIRE(overlay.rgpot_options().xtb_charge == Catch::Approx(-2.0));
}

TEST_CASE("RgpotPot cpmdc refuses a params_path it cannot read",
          "[PotTest][RGPOT][cpmdc]") {
#ifndef WITH_RGPOT
  SKIP("built without WITH_RGPOT");
#else
  // The file is read before any engine is loaded, so no libcpmdc is needed.
  clear_env("RGPOT_PARAMS_PATH");
  namespace fs = std::filesystem;
  const auto dir = fs::temp_directory_path() /
                   ("eon_rgpot_params_" + std::to_string(std::rand()));
  fs::create_directories(dir);

  Parameters params{};
  ParametersLoadAccess::potential_options(params).potential = PotType::RGPOT;
  ParametersLoadAccess::rgpot_options(params).backend = "cpmdc";

  SECTION("a missing file") {
    ParametersLoadAccess::rgpot_options(params).params_path =
        (dir / "absent.bin").string();
    REQUIRE_THROWS_WITH(eonc::helpers::makePotential(PotType::RGPOT, params),
                        ContainsSubstring("cannot open params_path"));
  }
  SECTION("a file that is not a Cap'n Proto message") {
    const auto path = dir / "seven-bytes.bin";
    {
      std::ofstream out(path, std::ios::binary);
      out << "1234567";
    }
    ParametersLoadAccess::rgpot_options(params).params_path = path.string();
    REQUIRE_THROWS_WITH(eonc::helpers::makePotential(PotType::RGPOT, params),
                        ContainsSubstring("not a capnp flat message"));
  }
  fs::remove_all(dir);
#endif
}

} // namespace tests

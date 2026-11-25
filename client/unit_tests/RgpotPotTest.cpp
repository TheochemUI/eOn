#include "catch2/catch_amalgamated.hpp"
#include "eon/CalculatorGroupUse.h"
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

TEST_CASE("Calculator-group use reports POP efficiencies", "[RGPOT][groups]") {
  eonc::CalculatorGroupUse use;
  // Four groups inside 100 s of driver wall: 80, 60, 60, 40 s busy.
  use.busy = {80.0, 60.0, 60.0, 40.0};
  use.systems = {10.0, 8.0, 8.0, 6.0};
  use.wall = 100.0;
  use.singleWall = 12.0;
  use.batches = 6;
  use.singles = 2;
  REQUIRE(use.groups() == 4);
  REQUIRE(use.meanBusy() == Catch::Approx(60.0));
  REQUIRE(use.loadBalance() == Catch::Approx(0.75));
  REQUIRE(use.communicationEfficiency() == Catch::Approx(0.8));
  REQUIRE(use.parallelEfficiency() == Catch::Approx(0.6));
  REQUIRE(use.idleFraction(0) == Catch::Approx(0.2));
  REQUIRE(use.idleFraction(3) == Catch::Approx(0.6));
  const std::string t = use.table();
  REQUIRE_THAT(t, ContainsSubstring("4 groups, 6 batches and 2 single"));
  REQUIRE_THAT(t, ContainsSubstring("load balance 0.750"));
  REQUIRE_THAT(t, ContainsSubstring("parallel efficiency 0.600"));

  // No grouped work: every ratio is the neutral 1.
  eonc::CalculatorGroupUse idle;
  REQUIRE(idle.loadBalance() == 1.0);
  REQUIRE(idle.communicationEfficiency() == 1.0);
  REQUIRE(idle.idleFraction(0) == 0.0);
}

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

namespace {

// The fixture table is synthetic and its functions are written in its header:
// F(rho) = -1.5 (1 - e^-rho), rho(r) = e^-1.2 r t(r), phi(r) = 2 e^-1.8 r t(r),
// t(r) = (1 - r/rc)^3 for r < rc = 4. The reference below never reads the
// file, so the engine and the formulas agree only if the file was read.
struct EamTerm {
  double value;
  double slope;
};

EamTerm eamTaper(double r) {
  constexpr double rc = 4.0;
  if (r >= rc)
    return {0.0, 0.0};
  const double s = 1.0 - r / rc;
  return {s * s * s, -3.0 * s * s / rc};
}

EamTerm eamDensity(double r) {
  const auto t = eamTaper(r);
  const double e = std::exp(-1.2 * r);
  return {e * t.value, e * (t.slope - 1.2 * t.value)};
}

EamTerm eamPair(double r) {
  const auto t = eamTaper(r);
  const double e = 2.0 * std::exp(-1.8 * r);
  return {e * t.value, e * (t.slope - 1.8 * t.value)};
}

EamTerm eamEmbedding(double rho) {
  return {-1.5 * (1.0 - std::exp(-rho)), -1.5 * std::exp(-rho)};
}

double eamReference(const AtomMatrix &pos, AtomMatrix &forces) {
  const int n = static_cast<int>(pos.rows());
  Eigen::VectorXd rho = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (i != j)
        rho(i) += eamDensity((pos.row(i) - pos.row(j)).norm()).value;
  double energy = 0.0;
  for (int i = 0; i < n; ++i)
    energy += eamEmbedding(rho(i)).value;
  forces = AtomMatrix::Zero(n, 3);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const Eigen::RowVector3d d = pos.row(i) - pos.row(j);
      const double r = d.norm();
      const double slope = eamPair(r).slope + (eamEmbedding(rho(i)).slope +
                                               eamEmbedding(rho(j)).slope) *
                                                  eamDensity(r).slope;
      energy += eamPair(r).value;
      forces.row(i) -= slope * d / r;
      forces.row(j) += slope * d / r;
    }
  }
  return energy;
}

std::string eamIni(const std::string &table, const std::string &format) {
  return "[Potential]\npotential = rgpot\n\n[RgpotPot]\nbackend = eam\n"
         "table_path = " +
         table + "\ntable_format = " + format + "\nneighbor_skin = 0.5\n";
}

std::string eamFixturePath() {
  const char *dir = std::getenv("EON_TEST_SYSTEMS_DIR");
  REQUIRE(dir != nullptr);
  return (std::filesystem::path(dir) / "rgpot_eam_table" /
          "synthetic_one_element.setfl")
      .string();
}

} // namespace

TEST_CASE("RgpotPot eam reads table_path, table_format and neighbor_skin "
          "from the INI",
          "[params][ini][RGPOT][eam]") {
  Parameters defaults;
  REQUIRE(defaults.rgpot_options().table_path.empty());
  REQUIRE(defaults.rgpot_options().table_format == "setfl");
  REQUIRE(defaults.rgpot_options().neighbor_skin == Catch::Approx(1.0));

  Parameters p;
  REQUIRE(p.load_ini_text(eamIni("/data/any.eam.alloy", "eam/fs")) == 0);
  REQUIRE(p.potential_options().potential == PotType::RGPOT);
  REQUIRE(p.rgpot_options().backend == "eam");
  REQUIRE(p.rgpot_options().table_path == "/data/any.eam.alloy");
  REQUIRE(p.rgpot_options().table_format == "eam/fs");
  REQUIRE(p.rgpot_options().neighbor_skin == Catch::Approx(0.5));
}

TEST_CASE("RgpotPot in-process eam evaluates the table named in the INI",
          "[PotTest][RGPOT][eam]") {
#ifndef WITH_RGPOT
  SKIP("built without WITH_RGPOT");
#else
  clear_env("RGPOT_BACKEND");
  Parameters params;
  REQUIRE(params.load_ini_text(eamIni(eamFixturePath(), "setfl")) == 0);

  auto pot = eonc::helpers::sharePotential(eonc::helpers::makePotential(
      params.potential_options().potential, params));
  REQUIRE(pot != nullptr);
  REQUIRE(pot->getType() == PotType::RGPOT);

  constexpr int n = 3;
  AtomMatrix pos(n, 3);
  pos << 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.4, 1.3, 0.2;
  Eigen::Matrix3d cell = Eigen::Matrix3d::Identity() * 20.0;
  const int nrs[n] = {18, 18, 18};

  double energy = 0.0;
  AtomMatrix forces = AtomMatrix::Zero(n, 3);
  pot->force(n, pos.data(), nrs, forces.data(), &energy, nullptr, cell.data());

  AtomMatrix expected;
  const double expectedEnergy = eamReference(pos, expected);
  REQUIRE(std::isfinite(energy));
  REQUIRE(energy != 0.0);
  REQUIRE(energy == Catch::Approx(expectedEnergy).epsilon(1e-4));
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < 3; ++k)
      REQUIRE(forces(i, k) ==
              Catch::Approx(expected(i, k)).epsilon(1e-3).margin(1e-5));
  REQUIRE(forces.colwise().sum().norm() < 1e-10);

  // The forces are the gradient of the energy the engine reports.
  const double h = 1e-5;
  for (const auto &[atom, axis] :
       {std::pair{0, 0}, std::pair{1, 0}, std::pair{2, 1}}) {
    double ep = 0.0, em = 0.0;
    AtomMatrix scratch = AtomMatrix::Zero(n, 3);
    AtomMatrix shifted = pos;
    shifted(atom, axis) += h;
    pot->force(n, shifted.data(), nrs, scratch.data(), &ep, nullptr,
               cell.data());
    shifted(atom, axis) -= 2.0 * h;
    pot->force(n, shifted.data(), nrs, scratch.data(), &em, nullptr,
               cell.data());
    REQUIRE(forces(atom, axis) ==
            Catch::Approx(-(ep - em) / (2.0 * h)).epsilon(1e-4).margin(1e-6));
  }
#endif
}

TEST_CASE("RgpotPot eam refuses a table it cannot read or understand",
          "[PotTest][RGPOT][eam]") {
#ifndef WITH_RGPOT
  SKIP("built without WITH_RGPOT");
#else
  clear_env("RGPOT_BACKEND");
  Parameters absent;
  REQUIRE(absent.load_ini_text(eamIni("/nonexistent/table.setfl", "setfl")) ==
          0);
  REQUIRE_THROWS_WITH(eonc::helpers::makePotential(PotType::RGPOT, absent),
                      ContainsSubstring("cannot open"));

  Parameters unknown;
  REQUIRE(unknown.load_ini_text(eamIni(eamFixturePath(), "no-such-layout")) ==
          0);
  REQUIRE_THROWS_WITH(eonc::helpers::makePotential(PotType::RGPOT, unknown),
                      ContainsSubstring("unknown table format"));
#endif
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

/*
** Integration: RgpotPot against potserv + fake engines (env-driven).
*/
#include "../MatrixHelpers.hpp"
#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include "catch2/catch_amalgamated.hpp"

#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>

using Catch::Matchers::WithinAbs;

TEST_CASE("RgpotPot force via potserv (requires WITH_RGPOT_CLIENT + server)",
          "[rgpot][integration]") {
#ifndef WITH_RGPOT_CLIENT
  SKIP("built without with_rgpot");
#else
  const char *host = std::getenv("RGPOT_POTSERV_HOST");
  const char *port_s = std::getenv("RGPOT_POTSERV_PORT");
  if (!host || !port_s) {
    SKIP("set RGPOT_POTSERV_HOST and RGPOT_POTSERV_PORT to potserv");
  }

  Parameters params;
  params.potential_options.potential = PotType::RGPOT;
  params.rgpot_options.host = host;
  params.rgpot_options.port = std::atoi(port_s);
  const char *backend = std::getenv("RGPOT_BACKEND");
  params.rgpot_options.backend = backend ? backend : "NWChem";
  // Optional: RGPOT_NWCHEM_BASIS / THEORY / SCF_TYPE (e.g. 6-31g* / dft / b3lyp)
  if (const char *b = std::getenv("RGPOT_NWCHEM_BASIS"))
    params.rgpot_options.nwchem_basis = b;
  else
    params.rgpot_options.nwchem_basis = "sto-3g";
  if (const char *t = std::getenv("RGPOT_NWCHEM_THEORY"))
    params.rgpot_options.nwchem_theory = t;
  else
    params.rgpot_options.nwchem_theory = "scf";
  if (const char *s = std::getenv("RGPOT_NWCHEM_SCF_TYPE"))
    params.rgpot_options.nwchem_scf_type = s;
  else
    params.rgpot_options.nwchem_scf_type = "rhf";
  params.rgpot_options.cpmd_functional = "BLYP";
  params.rgpot_options.cpmd_task = "gradient";

  auto pot = eonc::helpers::makePotential(params);
  REQUIRE(pot != nullptr);
  REQUIRE(pot->getType() == PotType::RGPOT);

  // Water geometry (Å) — matches rgpot NWChem B3LYP demos when env sets method
  const long N = 3;
  double R[9] = {0.0, 0.0, 0.11779, 0.0, 0.75545, -0.47116,
                 0.0, -0.75545, -0.47116};
  int Z[3] = {8, 1, 1};
  double F[9] = {};
  double U = 0.0;
  double var = 0.0;
  double box[9] = {100, 0, 0, 0, 100, 0, 0, 0, 100};

  pot->force(N, R, Z, F, &U, &var, box);
  REQUIRE(std::isfinite(U));
  REQUIRE(U != 0.0);
  bool any_f = false;
  for (int i = 0; i < 9; ++i) {
    REQUIRE(std::isfinite(F[i]));
    if (std::abs(F[i]) > 1e-12)
      any_f = true;
  }
  REQUIRE(any_f);

  const bool want_b3lyp =
      (std::getenv("RGPOT_NWCHEM_SCF_TYPE") &&
       std::string(std::getenv("RGPOT_NWCHEM_SCF_TYPE")).find("b3lyp") == 0) ||
      (std::getenv("RGPOT_NWCHEM_THEORY") &&
       std::string(std::getenv("RGPOT_NWCHEM_THEORY")).find("b3lyp") == 0);
  if (want_b3lyp) {
    // Water B3LYP/6-31G* ~ -2079 eV (not LDA ~-2063)
    REQUIRE(U < -2070.0);
  }

  // Second call (configure once). Real DFT may differ slightly between SCFs.
  double U2 = 0.0;
  pot->force(N, R, Z, F, &U2, &var, box);
  REQUIRE(std::isfinite(U2));
  const bool dft_like =
      want_b3lyp ||
      (std::getenv("RGPOT_NWCHEM_THEORY") != nullptr &&
       std::string(std::getenv("RGPOT_NWCHEM_THEORY")) == "dft");
  const double etol = dft_like ? 1e-3 : 1e-9;
  REQUIRE_THAT(U2, WithinAbs(U, etol));
  if (want_b3lyp)
    REQUIRE(U2 < -2070.0);
#endif
}

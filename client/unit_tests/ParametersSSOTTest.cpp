#include "eon/ParametersSSOT.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace eonc::config {
JobType job_from_ssot(std::string_view);
PotType pot_from_ssot(std::string_view);
OptType opt_from_ssot(std::string_view);
} // namespace eonc::config

namespace tests {
static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("Parameters constructor applies Cap'n Proto SSoT defaults",
          "[params][ssot]") {
  Parameters p;
  // From schema/eon_params.capnp MainOptions / PotentialOptions
  REQUIRE(p.main_options().temperature == Catch::Approx(300.0));
  REQUIRE(p.main_options().finiteDifference == Catch::Approx(0.01));
  REQUIRE(p.main_options().job == JobType::Process_Search);
  REQUIRE(p.potential_options().potential == PotType::LJ);
  REQUIRE(p.optimizer_options().method == OptType::CG);
  REQUIRE(p.optimizer_options().max_iterations == 1000);
  REQUIRE(p.process_search_options().minimize_first == true);
  REQUIRE(p.structure_comparison_options().distance_difference ==
          Catch::Approx(0.1));
  REQUIRE(p.rgpot_options().backend == "nwchemc");
  REQUIRE(p.rgpot_options().cutoff_ry == Catch::Approx(70.0));
}

TEST_CASE("schema names map onto job, potential, and optimizer types",
          "[params][ssot][names]") {
  using eonc::config::job_from_ssot;
  using eonc::config::opt_from_ssot;
  using eonc::config::pot_from_ssot;
  REQUIRE(job_from_ssot("process_search") == JobType::Process_Search);
  REQUIRE(job_from_ssot("minimization") == JobType::Minimization);
  REQUIRE(job_from_ssot("saddle_search") == JobType::Saddle_Search);
  REQUIRE(job_from_ssot("basin_hopping") == JobType::Basin_Hopping);
  REQUIRE(job_from_ssot("parallel_replica") == JobType::Parallel_Replica);
  REQUIRE(job_from_ssot("unbiased_parallel_replica") ==
          JobType::Parallel_Replica);
  REQUIRE(job_from_ssot("nudged_elastic_band") == JobType::Nudged_Elastic_Band);
  REQUIRE(job_from_ssot("dynamics") == JobType::Dynamics);
  REQUIRE(job_from_ssot("molecular_dynamics") == JobType::Dynamics);
  REQUIRE(job_from_ssot("hessian") == JobType::Hessian);
  REQUIRE(job_from_ssot("point") == JobType::Point);
  REQUIRE(job_from_ssot("prefactor") == JobType::Prefactor);
  REQUIRE(job_from_ssot("monte_carlo") == JobType::Monte_Carlo);
  REQUIRE(job_from_ssot("structure_comparison") ==
          JobType::Structure_Comparison);
  REQUIRE(job_from_ssot("gp_surrogate") == JobType::GP_Surrogate);
  REQUIRE(job_from_ssot("safe_hyperdynamics") == JobType::Safe_Hyperdynamics);
  REQUIRE(job_from_ssot("tad") == JobType::TAD);
  REQUIRE(job_from_ssot("replica_exchange") == JobType::Replica_Exchange);
  REQUIRE(job_from_ssot("finite_difference") == JobType::Finite_Difference);
  REQUIRE(job_from_ssot("finite_differences") == JobType::Finite_Difference);
  REQUIRE(job_from_ssot("global_optimization") == JobType::Global_Optimization);
  REQUIRE(job_from_ssot("instanton") == JobType::Instanton);
  REQUIRE(job_from_ssot("no-such-job") == JobType::Process_Search);
  REQUIRE(pot_from_ssot("lj") == PotType::LJ);
  REQUIRE(pot_from_ssot("eam_al") == PotType::EAM_AL);
  REQUIRE(pot_from_ssot("emt") == PotType::EMT);
  REQUIRE(pot_from_ssot("lammps") == PotType::LAMMPS);
  REQUIRE(pot_from_ssot("morse_pt") == PotType::MORSE_PT);
  REQUIRE(pot_from_ssot("metatomic") == PotType::METATOMIC);
  REQUIRE(pot_from_ssot("xtb") == PotType::XTB);
  REQUIRE(pot_from_ssot("rgpot") == PotType::RGPOT);
  REQUIRE(pot_from_ssot("ext_pot") == PotType::EXT_POT);
  REQUIRE(pot_from_ssot("no-such-pot") == PotType::LJ);
  REQUIRE(opt_from_ssot("cg") == OptType::CG);
  REQUIRE(opt_from_ssot("lbfgs") == OptType::LBFGS);
  REQUIRE(opt_from_ssot("qm") == OptType::QM);
  REQUIRE(opt_from_ssot("quickmin") == OptType::QM);
  REQUIRE(opt_from_ssot("sd") == OptType::SD);
  REQUIRE(opt_from_ssot("fire") == OptType::FIRE);
  REQUIRE(opt_from_ssot("xtsci") == OptType::XTSCI);
  REQUIRE(opt_from_ssot("no-such-opt") == OptType::CG);
}

TEST_CASE("ssot_has_field knows covered catalog keys", "[params][ssot]") {
  REQUIRE(eonc::config::ssot_has_field("Main", "temperature"));
  REQUIRE(eonc::config::ssot_has_field("Potential", "potential"));
  REQUIRE(eonc::config::ssot_has_field("RgpotPot", "cutoff_ry"));
  REQUIRE(eonc::config::ssot_has_field("RgpotPot", "cutOffRy"));
  REQUIRE(eonc::config::ssot_has_field("RgpotPot", "params_path"));
  REQUIRE_FALSE(eonc::config::ssot_has_field("Main", "not_a_real_field"));
  REQUIRE_FALSE(eonc::config::ssot_has_field("Dimer", "rotation_angle"));
}

TEST_CASE("Parameters::load INI still overrides SSoT defaults",
          "[params][ssot][ini]") {
  namespace fs = std::filesystem;
  auto dir = fs::temp_directory_path() / "eon_ssot_ini_XXXXXX";
  // unique path
  dir = fs::temp_directory_path() /
        ("eon_ssot_ini_" + std::to_string(std::rand()));
  fs::create_directories(dir);
  auto ini = dir / "config.ini";
  {
    std::ofstream out(ini);
    out << "[Main]\njob = minimization\ntemperature = 450.0\n"
           "random_seed = 7\n"
           "[Potential]\npotential = eam_al\n"
           "[Optimizer]\nopt_method = lbfgs\nmax_iterations = 42\n";
  }
  Parameters p;
  REQUIRE(p.load(ini.string()) == 0);
  REQUIRE(p.main_options().job == JobType::Minimization);
  REQUIRE(p.main_options().temperature == Catch::Approx(450.0));
  REQUIRE(p.main_options().randomSeed == 7);
  REQUIRE(p.potential_options().potential == PotType::EAM_AL);
  REQUIRE(p.optimizer_options().method == OptType::LBFGS);
  REQUIRE(p.optimizer_options().max_iterations == 42);
  REQUIRE(p.last_load_source() == ini.string());
  REQUIRE(p.last_load_error() == 0);
  Parameters copied = p;
  REQUIRE(copied.last_load_source() == p.last_load_source());
  REQUIRE(copied.last_load_error() == 0);
  REQUIRE(copied.main_options().temperature == Catch::Approx(450.0));
  Parameters assigned;
  assigned = p;
  REQUIRE(assigned.last_load_source() == p.last_load_source());
  fs::remove_all(dir);
}

TEST_CASE("Parameters::load accepts AMS DFTB and FORCEFIELD engines",
          "[params][ams]") {
  namespace fs = std::filesystem;
  auto dir = fs::temp_directory_path() /
             ("eon_ams_engine_" + std::to_string(std::rand()));
  fs::create_directories(dir);
  const auto ini = dir / "config.ini";
  const auto load_body = [&](const char *body) {
    {
      std::ofstream out(ini);
      out << body;
    }
    Parameters loaded;
    const int rc = loaded.load(ini.string());
    return std::make_pair(rc, std::move(loaded));
  };

  {
    const auto [rc, loaded] = load_body("[Potential]\n"
                                        "potential = ams\n"
                                        "[AMS]\n"
                                        "engine = DFTB\n"
                                        "resources = DFTB\n");
    REQUIRE(rc == 0);
    REQUIRE(loaded.potential_options().potential == PotType::AMS);
    REQUIRE(loaded.ams_options().engine == "DFTB");
    REQUIRE(loaded.ams_options().resources == "DFTB");
    REQUIRE(loaded.ams_options().forcefield.empty());
    REQUIRE(loaded.ams_options().model.empty());
    REQUIRE(loaded.ams_options().xc.empty());
    REQUIRE(loaded.last_load_error() == 0);
  }
  {
    const auto [rc, loaded] = load_body("[Potential]\n"
                                        "potential = ams\n"
                                        "[AMS]\n"
                                        "engine = forcefield\n");
    REQUIRE(rc == 0);
    REQUIRE(loaded.ams_options().engine == "forcefield");
    REQUIRE(loaded.last_load_error() == 0);
  }
  {
    const auto [rc, loaded] = load_body("[Potential]\n"
                                        "potential = ams_io\n"
                                        "[AMS_IO]\n"
                                        "engine = FORCEFIELD\n");
    REQUIRE(rc == 0);
    REQUIRE(loaded.potential_options().potential == PotType::AMS_IO);
    REQUIRE(loaded.ams_options().engine == "FORCEFIELD");
    REQUIRE(loaded.last_load_error() == 0);
  }
  {
    const auto [rc, loaded] = load_body("[Potential]\n"
                                        "potential = ams\n"
                                        "[AMS]\n"
                                        "engine = DFTB\n");
    REQUIRE(rc != 0);
    REQUIRE(loaded.last_load_error() != 0);
  }
  {
    const auto rejected = load_body("[Potential]\n"
                                    "potential = ams\n"
                                    "[AMS]\n"
                                    "engine = MOPAC\n");
    REQUIRE(rejected.first != 0);
  }
  {
    const auto [rc, loaded] = load_body("[Potential]\n"
                                        "potential = ams\n"
                                        "[AMS]\n"
                                        "engine = MOPAC\n"
                                        "model = PM3\n");
    REQUIRE(rc == 0);
    REQUIRE(loaded.ams_options().model == "PM3");
  }
  {
    const auto rejected = load_body("[Potential]\n"
                                    "potential = ams\n"
                                    "[AMS]\n"
                                    "engine = MOPAC\n"
                                    "forcefield = ff\n"
                                    "model = PM3\n"
                                    "xc = LDA\n");
    REQUIRE(rejected.first != 0);
  }
  fs::remove_all(dir);
}

TEST_CASE("Parameters load-state Impl records a missing file",
          "[params][pimpl]") {
  Parameters p;
  REQUIRE(p.last_load_source().empty());
  REQUIRE(p.last_load_error() == 0);
  REQUIRE(p.load("no-such-eon-config.ini") != 0);
  REQUIRE(p.last_load_source() == "no-such-eon-config.ini");
  REQUIRE(p.last_load_error() != 0);
  Parameters moved = std::move(p);
  REQUIRE(moved.last_load_source() == "no-such-eon-config.ini");
  REQUIRE(p.last_load_source().empty());
}

TEST_CASE("Parameters::load rejects INI parse errors", "[params][ini]") {
  Parameters memory;
  const std::string missing_equals =
      "[Main]\ntemperature 450\njob = minimization\n";
  REQUIRE(memory.load_ini_text(missing_equals) != 0);
  REQUIRE(memory.last_load_source() == "<ini>");
  REQUIRE(memory.last_load_error() != 0);
  REQUIRE(memory.main_options().job == JobType::Process_Search);
  REQUIRE(memory.main_options().temperature == Catch::Approx(300.0));

  namespace fs = std::filesystem;
  const auto dir = fs::temp_directory_path() /
                   ("eon_ini_parse_" + std::to_string(std::rand()));
  fs::create_directories(dir);
  const auto path = dir / "config.ini";
  {
    std::ofstream out(path);
    out << "[Main]\njob minimization\n";
  }
  Parameters from_path;
  REQUIRE(from_path.load(path.string()) != 0);
  REQUIRE(from_path.last_load_error() != 0);
  REQUIRE(from_path.main_options().job == JobType::Process_Search);

  {
    std::ofstream out(path);
    // Longer than INI_MAX_LINE (65536).
    out << "[Main]\ntemperature = 450" << std::string(70000, '0') << '\n';
  }
  Parameters overlong;
  REQUIRE(overlong.load(path.string()) != 0);
  REQUIRE(overlong.last_load_error() != 0);
  REQUIRE(overlong.main_options().temperature == Catch::Approx(300.0));

  FILE *handle = std::fopen(path.string().c_str(), "rb");
  REQUIRE(handle != nullptr);
  Parameters from_file;
  REQUIRE(from_file.load(handle) != 0);
  REQUIRE(from_file.last_load_source() == "<FILE*>");
  REQUIRE(from_file.last_load_error() != 0);
  REQUIRE(from_file.main_options().temperature == Catch::Approx(300.0));
  std::fclose(handle);
  fs::remove_all(dir);
}

TEST_CASE("Parameters INI reads instanton, OH-TST, and socket keys",
          "[params][ini]") {
  Parameters p;
  REQUIRE(p.load_ini_text("[Main]\n"
                          "job = instanton\n"
                          "temperature = 250\n"
                          "[Potential]\n"
                          "potential = lj\n"
                          "lammps_threads = 2\n"
                          "lammps_logging = true\n"
                          "[Instanton]\n"
                          "mode = rate\n"
                          "springs = eco\n"
                          "initial_hessians = finite_difference\n"
                          "friction = explicit\n"
                          "temperatures = 200, 300\n"
                          "discretization = 1, 2\n"
                          "friction_eta_beads = 0.1, 0.2\n"
                          "beads = 16\n"
                          "[OH_TST]\n"
                          "thermostat = gle\n"
                          "equil_steps = 12\n"
                          "sample_steps = 24\n"
                          "gle_a_file = drift.txt\n"
                          "[SocketNWChemPot]\n"
                          "unix_socket_mode = true\n"
                          "unix_socket_path = eon_sock\n"
                          "mem_in_gb = 2\n"
                          "[Hyperdynamics]\n"
                          "bias_potential = bond_boost\n") == 0);
  REQUIRE(p.instanton_options().beads == 16);
  REQUIRE(p.instanton_options().springs == "eco");
  REQUIRE(p.instanton_options().temperatures.size() == 2);
  REQUIRE(p.oh_tst_options().equil_steps == 12);
  REQUIRE(p.oh_tst_options().gle_a_file == "drift.txt");
  REQUIRE(p.socket_nwchem_options().unix_socket_mode);
  REQUIRE(p.socket_nwchem_options().unix_socket_path == "eon_sock");
  REQUIRE(p.potential_options().LAMMPSThreads == 2);
}

TEST_CASE("Parameters INI reads the optional potential sections",
          "[params][ini]") {
  Parameters xtb;
  REQUIRE(xtb.load_ini_text("[Potential]\npotential = xtb\n"
                            "[XTBPot]\nparamset = GFN1\naccuracy = 0.2\n"
                            "charge = 1\n") == 0);
  REQUIRE(xtb.xtb_options().paramset == "GFN1");
  REQUIRE(xtb.xtb_options().acc == Catch::Approx(0.2));

  Parameters zbl;
  REQUIRE(zbl.load_ini_text("[Potential]\npotential = zbl\n"
                            "[ZBLPot]\ncut_inner = 0.5\ncut_global = 2.5\n") ==
          0);
  REQUIRE(zbl.zbl_options().cut_inner == Catch::Approx(0.5));
  REQUIRE(zbl.zbl_options().cut_global == Catch::Approx(2.5));

  Parameters crossed;
  REQUIRE_THROWS_AS(
      crossed.load_ini_text("[Potential]\npotential = zbl\n"
                            "[ZBLPot]\ncut_inner = 3\ncut_global = 1\n"),
      std::runtime_error);

  Parameters dispersion;
  REQUIRE(dispersion.load_ini_text("[Potential]\npotential = dftd3\n"
                                   "[D3Pot]\nfunctional = blyp\natm = false\n"
                                   "damping = zero\n") == 0);
  REQUIRE(dispersion.dftd_options().functional == "blyp");
  REQUIRE_FALSE(dispersion.dftd_options().atm);
  REQUIRE(dispersion.dftd_options().d3_damping == "zero");

  Parameters expr;
  REQUIRE(expr.load_ini_text("[Potential]\npotential = expr\n"
                             "[ExprPot]\nexpression = lj\nterms = pair\n") ==
          0);
  REQUIRE(expr.expr_options().expression == "lj");
  REQUIRE(expr.expr_options().terms == "pair");

  Parameters mopac;
  REQUIRE(mopac.load_ini_text(
              "[Potential]\npotential = mopac\n"
              "[MOPACPot]\ncharge = -1\nspin = 1\nmodel = 2\n"
              "engine_path = libmopac.so\n") == 0);
  REQUIRE(mopac.mopac_options().charge == -1);
  REQUIRE(mopac.mopac_options().engine_path == "libmopac.so");

  Parameters ams;
  REQUIRE(ams.load_ini_text("[Potential]\npotential = ams\n"
                            "[AMS]\nengine = dftb\nforcefield = uff\n"
                            "[AMS_ENV]\namshome = /opt/ams\n") == 0);
  REQUIRE(ams.ams_options().engine == "dftb");
  REQUIRE(ams.ams_options().env.amshome == "/opt/ams");

  Parameters sections;
  REQUIRE(sections.load_ini_text(
              "[Potential]\npotential = lj\n"
              "[CG]\ncg_line_search = true\ncg_max_iter_line_search = 7\n"
              "[FIRE]\ntime_step = 0.5\ntime_step_max = 2\n"
              "[Xtsci]\nmethod = lbfgs\nqn_step = full\nprecon = none\n"
              "accept = ratio\nhighs = true\nmanifold = off\n"
              "[ASE_ORCA]\norca_path = orca\nnproc = 2\ncharge = -1\n"
              "[ASE_NWCHEM]\nnwchem_path = nwchem\nbasis = 6-31g\n"
              "memory = 400\n"
              "[Metatomic]\nmodel_path = model.pt\ndevice = cpu\n"
              "check_consistency = true\nn_symmetry_rotations = 3\n"
              "variant_base = base\n"
              "[Serve]\nhost = 127.0.0.1\nport = 4321\nreplicas = 2\n"
              "endpoints = lj:1\n"
              "[Saddle Search]\ndisplace_atom_list = 1, 2\n"
              "confine_positive = true\nbowl_active_atoms = 4\n"
              "confine_positive_min_active = 2\n") == 0);
  REQUIRE(sections.optimizer_options().cg.line_search);
  REQUIRE(sections.optimizer_options().cg.line_search_max_iter == 7);
  REQUIRE(sections.optimizer_options().time_step_input == Catch::Approx(0.5));
  REQUIRE(sections.optimizer_options().xtsci.method == "lbfgs");
  REQUIRE(sections.optimizer_options().xtsci.highs);
  REQUIRE(sections.ase_orca_options().charge == -1);
  REQUIRE(sections.ase_nwchem_options().basis == "6-31g");
  REQUIRE(sections.metatomic_options().model_path == "model.pt");
  REQUIRE(sections.metatomic_options().n_symmetry_rotations == 3);
  REQUIRE(sections.serve_options().host == "127.0.0.1");
  REQUIRE(sections.serve_options().port == 4321);
  REQUIRE(sections.saddle_search_options().displace_atom_list.size() == 2);
  REQUIRE(sections.saddle_search_options().confine_positive.enabled);
  REQUIRE(sections.saddle_search_options().confine_positive.min_active == 2);
}

} // namespace tests
